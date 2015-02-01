# TODO: Fit a super learned model
# 
# Author: solomon
###############################################################################
setwd('~/Dropbox/creditClaimingProjects/het/HetPackage')
set.seed(23432)

source('HetEffects.R')

n <- 500
p <- 2
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X) <- paste("X", 1:p, sep="")
#X <- data.frame(X)
T1 <- c("C", "T")[rbinom(n, 1, .5) +1]
T2 <- c("C", "T")[rbinom(n, 1, .5) +1]
Y <- as.numeric(as.factor(T1)) * X[, 1] + 
    as.numeric(as.factor(T2))  
data = data.frame(Y, T1, T2, X)
formula = Y ~ T1*T2*X
treatments = c('T1', 'T2')
weights = NULL
method = 'method.NNLS'
id = NULL 
verbose = FALSE
na.action = NULL
control = list()
cvControl = list()
family=gaussian()
R = 2

# build Library and run Super Learner
SL.glmnet_0.5 <- function(..., alpha = 0.5){ SL.glmnet(..., alpha = 0.5)}


SL.library <- c("SL.glm", "SL.mean")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights

SL.library <- c("SL.glm", "SL.randomForest")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights

SL.library <- c("SL.glm", "SL.glmnet_0.5")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights

SL.library <- c("SL.glm", "SL.glmnet_0.5", "SL.randomForest")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights

SL.library <- c("SL.glm", "SL.glmnet_0.5", "SL.randomForest", "SL.gam")
hetmod = HetEffects(formula, data, treatments, SL.library = SL.library)
hetmod[[1]]$weights





##alright, beginning to 
SL.glmnet_0.5<- function(..., alpha = 0.5){ SL.glmnet(..., alpha = 0.5)}
ee<- SuperLearner(Y, X, SL.library = c("SL.glm", "SL.randomForest","SL.krls",  "SL.glmnet_0.5"))



n <- 500
p <- 50
t_p<- 3
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X) <- paste("X", 1:p, sep="")
X <- data.frame(X)
Treat<- matrix(NA, nrow=n, ncol=t_p)
for(z in 1:n){
	for(k in 1:3){
		Treat[z,k]<- rbinom(1, size = 1, p = 0.5)
		}
	}
	
Y <- pnorm(X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + 0.5*Treat[,1] + 0.25*Treat[,2] + -2*Treat[,3] + 5*X[,1]*Treat[,3]  + rnorm(n))

Y_bin<- c()
for(z in 1:length(Y)){
	Y_bin[z]<- rbinom(1, size =1, prob = Y[z])
	}

Y_bin[which(Y_bin==0)]<- -1	

mkstand<- function(x){
	x<- x - mean(x, na.rm=T)
	if(sd(x, na.rm=T)>0){
		x<- x/sd(x, na.rm=T)
		}
	return(x)
	}
	


SDsToRescaleX<- apply(X, 2, sd)
SDsToRescaleX[SDsToRescaleX==0]<- 1
Xstd<- apply(X, 2, mkstand)
colnames(Treat)<- paste('Treat', 1:3)
dd<- FindIt(Y_bin, Xstd, Treat, type='multiple', scale.c = SDsToRescaleX, search.lambdas=TRUE, fit.glmnet=TRUE)

ert<- cbind(1, Xstd, Treat)%*%dd$coefs/2 + 0.5

  treatments <- attr(Xmat, "treatments")
  Treat_columns = c()
  for(t in treatments){
    Treat_columns <- c(Treat_columns, grep(paste('^', treatments[i],sep=''),  colnames(Xmat))) 
  }
  Treat_columns = Treat_columns[!duplicated(Treat_columns)]
  attr(Xmat, "Treat_columns") <- Treat_columns


