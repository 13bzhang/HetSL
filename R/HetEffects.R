# TODO: Fit a super learned model
# 
# Author: solomon
###############################################################################


HetEffects <- function(formula, 
    data, 
    treatments, 
    subset, 
    weights = NULL, 
    family = gaussian(),
    SL.library, 
    method = 'method.NNLS', 
    id = NULL, 
    verbose = FALSE, 
    control = list(), 
    cvControl = list(), 
    cores = 1,
    R = 2,
    na.action = NULL, ...){
  # TODO: allow people to pass covariate matrix + weights?
  # TODO: Consider setting tuning params on the fly
    # glmnet::alpha
    # knn::k
  
  require('SuperLearner')
  
  # vector for treatment column names. Check to make sure all is well:
  if (!is.character(treatments)) stop("'treatments' must be a character scalar/vector")
  if (!is.vector(treatments)) stop("'treatments' must be a character scalar/vector")
  for (i in 1:length(treatments)){
    if (length(grep(treatments[i], formula, fixed=TRUE))==0) {
      stop("All treatments must be in the formula")
    }
  }
  
  # Parse formula, data into model matrix, Y
  # TODO: Decide on what to do about NAs
  X_frame <- model.frame(formula, data=data)
  toInclude <- complete.cases(X_frame)
  Xmat <- model.matrix(formula, data=data)
  Xmat <- Xmat[toInclude,]
  Y <- X_frame[[1]]
  Y <- Y[toInclude]
  Xmat <- data.frame(Xmat)
  
  # Take all covariates, all treatments, interactions
  # rename treatment columns T_whatever so that we can grep
  # them in FindIt.
  
  attr(Xmat, "treatments") <- treatments
  
  # check to see if variables are numeric, if so put them in quintiles
  quintilize = function(x){
    qs = quantile(x, probs = c(0, .25, .5, .75, 1))
    if(length(unique(qs)) < 5) warning("Binning quantitative variable to < 5 bins")
    qs[1] = qs[1] - 0.01 # ensures cut doesn't set values at lowest quintile to NA   
    cut(x, unique(qs))
  }
  
  col_numeric <- apply(Xmat, 2, function(x) length(unique(x)) > 15 )
  Xmat_buckets = Xmat
  Xmat_buckets[col_numeric] = apply(Xmat_buckets[col_numeric], 2, 
      quintilize) 
#  anyNA = function(x) any(is.na(x))
#  apply(Xmat_buckets, 2, anyNA)
  
  # bootstrap
  # Generate matrix for prediction:
  veclist = apply(Xmat_buckets[,-1], 2, unique)
  veclist = lapply(veclist, na.omit)
  veclist = lapply(veclist, as.vector)
  if (prod(unlist(lapply(veclist, length))) > 1e10) {
    stop("Too many factor combinations to estimate effects")
  }
  newX = expand.grid(veclist, stringsAsFactors=FALSE)
  
  # Set columns in newX to numeric that are numeric in Xmat_buckets
  numeric_cols = unlist(lapply(Xmat_buckets, is.numeric))
  newX[,numeric_cols[-1]] = lapply(newX[,numeric_cols[-1]], as.numeric)
#  Xmat_buckets contrasts are different perhaps? 
  
  # bootstrap SuperLearner 
  # get effects for every combination of covariates
  # install.packages('doMC')
#  require('doMC')  # for easy parallelization on a nice multicore machine
#  registerDoMC(cores = cores)
#  bootres = foreach(i=1:R, .combine=rbind, .errorhandling='remove') %dopar% {	

  # Put everything in mod matrix first:
  Xmatb <- data.frame(model.matrix(~ -1 + ., Xmat_buckets[,-1]))
  newXb <- data.frame(model.matrix(~ -1 +., newX))
#  names(newXb) == names(Xmatb)
  
  res = list()
  for( i in 1:R){
    d = sample(x = 1:length(Xmatb[,1]), 
               size = length(Xmatb[,1]), 
               replace=TRUE)
    bsl = SuperLearner(Y = Y[d], 
            X = Xmatb[d,], # don't pass intercept
            newX = newXb, 
            family = family, 
            SL.library = SL.library, 
            method = method, 
            id = id, 
            verbose = verbose, 
            control = control, 
            cvControl = cvControl, 
            obsWeights = weights) 
     res[[i]] = list(weights = bsl$coef, predictions = bsl$library.predict)
#     return( list(weights = bsl$coef, predictions = bsl$library.predict) )
  }
  return(list(boostrap_samples = res, effectX = newXb))
}

