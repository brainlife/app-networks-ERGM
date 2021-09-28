#### 1 INPUT
setwd("/Users/vito.dichio/ownCloud/Vit-Fab/Coding/BrainLife")
## 1.1 INPUT-EXT
library("rjson")
config <- fromJSON(file = "config.json")
#for(i in 1:length(sys.argv)){config[i]=sys.argv[i]}

my_formula = as.character(config[1])
nsim_gof <- as.integer(config[2])
namefile = as.character(config[3])
num_ecov = as.integer(config[4])
num_ncov = as.integer(config[5])
unfiltered = as.logical(config[6])

## 1.2 INPUT-INT

# 1.2.1 Network
#ECO Neural Network Filter----
ECO <- function(data, k=3, directed = FALSE)
{
  #corrlation matrix
  if(nrow(data)==ncol(data)){cormat<-data
  }else{cormat<-cor(data)}
  C<-cormat
  n<-ncol(C)
  S<-C
  if(directed)
  {
    numcon<-k*n
    ind<-which(C!=0)
  }else{C<-upper.tri(C,diag=TRUE)
  numcon<-k/2*n
  ind<-which(upper.tri(C,diag=TRUE)!=0)}
  S<-ifelse(C==1,S,0)
  if(numcon>length(ind))
  {
    stop("Input matrix is too sparse")
  }
  sorind<-matrix(0,nrow=length(ind),ncol=2)
  G<-S
  S<-abs(S)
  x<-S[ind]
  y<-ind
  h<-cbind(ind,S[ind])
  sorind<-h[order(-h[,2]),]
  C[sorind[(numcon+1):nrow(sorind),1]]<-0
  if(directed)
  {W<-C}else{W<-C+t(C)
  diag(W)<-1}
  J<-G+t(G)
  diag(J)<-1
  W<-ifelse(W!=0,J,0)
  W<-as.data.frame(W)
  colnames(W)<-colnames(data)
  row.names(W)<-colnames(data)
  W<-as.matrix(W)
  return(W)
}
#----

data <- read.table(namefile, sep="")
if (unfiltered == TRUE){
  RB = ECO(as.matrix(data), 9, FALSE)
  data = ifelse(RB > 0, 1, 0)
}

library(network)
bnet = network(data, directed = FALSE)

# 1.2.2 Edge-Covariates
if (num_ecov != 0){
ecov = vector("list", length = num_ecov)
for (i in 1:num_ecov){
ecov[[i]] = as.matrix( read.table(paste("ecov_",i,".txt",sep=""),header=FALSE) )
if (length(ecov[[i]] != network.size(bnet))) {
  stop("Edgecov is not of size [#nodes]x[#nodes]!!!")
} 
}
}

# 1.2.3 Node-Covariates
if (num_ncov != 0){
attrs <- read.table("attrs.txt",header=TRUE,stringsAsFactors=FALSE)
for (i in 1:num_ncov){
if (length(attrs[[1]] != network.size(bnet))) {
  stop("Attribute is not of size [#nodes]!!!")
  }  
bnet%v%names(attrs)[i] <- attrs[[i]]
}
}

#### 2 ERGM 
library(ergm)
my_log <- file("output/log-computation.txt")
sink(my_log, append = TRUE, type = "output")

## 2.2 ESTIMATION 
summary(bnet ~ edges + triangle + degree(1:10))
bfit <- eval(parse(text = paste("ergm(bnet ~", my_formula, ")")))
pdf("output/mcmc-diagnostic.pdf")
mcmc.diagnostics(bfit)
dev.off()

## 2.3 GOF 
bfit.gof <- gof(bfit,control=control.gof.formula(nsim = nsim_gof)) #output simulate 
bfit.gof

pdf("output/gof.pdf")
plot(bfit.gof)
dev.off()

#### 3 OUTPUT 

closeAllConnections()

my_log_est <- file("output/estimation.txt")
sink(my_log_est, append = TRUE, type = "output")
print("Estimated values and covariance matrix:")
estimate = list(bfit$coefficients,bfit$covar)
print(estimate)
print("Summary of the fitting procedure:")
summary(bfit)


closeAllConnections()



