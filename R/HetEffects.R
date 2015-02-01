# Fits a super learned ensamble to estimate heterogenous treatment effects.
# Key decisions and TODOs:
#   - Do we want users to be able to specify newX themselves?
#     - Do we want to force users to do this themselves,
#       and just provide helper functions?
#   - How to handle continuous variables?
#     - If estimating MCATEs, do we estimate at the mean?
#     - If not, do we fix at the mean?
#   - Still need to process and format results, display for the user.
#     We can base this on previous code.
#   - Test, test, test.
#   -
#   -
###############################################################################

#' Estimating Heterogenous Treatments and Effects using Super Learning
#'
#' The 'HetEffects' command uses the method outlined in Grimmer, Messing
#' and Westwood (2014) to estimate heterogenous treatments and heterogenous
#' treatment effects.  This method uses an ensemble of constituent methods,
#' weighted based on cross-validated model performance.
#' An ensemble of constituent models are fit based on user-supplied formulas,
#' weighted according to cross-validated model performance, then used to
#' predict the expected marginal conditional averages.  Then the treatment
#' and control are differenced, conditional on key user-specified covariates
#' to produce Marginal Conditional Average Treatment Effects (MCATEs).
#'
#'
#' @usage
#' HetEffects(
#'    formula,
#'    data,
#'    treatments,
#'    subset,
#'    weights = NULL,
#'    family = gaussian(),
#'    SL.library,
#'    method = 'method.NNLS',
#'    id = NULL,
#'    verbose = TRUE,
#'    control = list(),
#'    cvControl = list(),
#'    bootstrap = TRUE,
#'    cores = 1,
#'    R = 2,
#'    na.action = NULL, ...)
#'
#' @param formula The formula used to estimate hetergenous effects
#' @param data The data to be used.  Data MUST be in XYZ format.
#' @param treatments Specify which of columns in the data are treatments.
#' @param subset The subset of variables for which to return MCATEs.
#'        Setting to NULL will return all variables in the formula.
#' @param weights observation weights to be used to fit models.
#'        Verify SuperLearner can utilize weights in each model.
#' @param Family link function for models. Currently either gaussian() or binomial().
#' @param SL.library library of algorithms
#' @param method Loss function and model to estimate ensemble weights.
#'        Currently method.NNLS, method.NNloglik, or custom
#'        (see SuperLearner::create.method.template()).
#' @param id cluster id
#' @param verbose View progress of estimation? Defaults to TRUE.
#' @param control Optional controls for the SuperLearner package.
#' @param cvControl List for CV control, see SuperLearner documentation.
#' @param bootstrap Bootstrap estimates?  Defaults to FALSE.
#' @param cores Number of cores to use if using doMC for parallelization.
#' @param R Number of bootstrap replicate samples to estimate if bootstrap == TRUE.
#' @param na.action How to handle NA values.
#' @param \dots additional parameters to be passed to the plot
#'
#' @return if bootstrap is FALSE, returns a list of model weights and
#'         predictions for each model, and the corresponding covariate
#'         combinations, effectX.  If bootstrap is TRUE,
#'         returns a list of the above (of length 1:R) for each boostrap
#'         replicate, and the corresponding covariate
#'         combinations, effectX
#' @export
#' @author Solomon Messing and Sean Westwood and Justin Grimmer

HetEffects <- function(formula,
    data,
    treatments,
    subset,
    weights = NULL,
    family = gaussian(),
    SL.library,
    method = 'method.NNLS',
    id = NULL,
    verbose = TRUE,
    control = list(),
    cvControl = list(),
    bootstrap = TRUE,
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
  X_frame <- model.frame(formula, data=data)
  toInclude <- complete.cases(X_frame)
  Xmat <- model.matrix(formula, data=data)
  Xmat <- Xmat[toInclude,]
  Y <- X_frame[[1]]
  Y <- Y[toInclude]
  Xmat <- data.frame(Xmat)

  # TODO: Take all covariates, all treatments, interactions
  # rename treatment columns T_whatever so that we can grep
  # them in FindIt.

  attr(Xmat, "treatments") <- treatments


  # Check to see if variables are numeric, if so put them in quintiles

  # TODO: find a better way to do this...
  col_numeric <- apply(Xmat, 2, function(x) length(unique(x)) > 15 )

  # SM: Do we really want to do this?  We might be better off forcing
  # users to do this themselves.  What we have now could yield unexpected
  # behavior.

  # TODO: Instead perhaps we should move this to a new function,
  # and allow users to customize this.  We should then allow
  # users to pass the results as an argument to HetEffects, and/or
  # call the function here with sensible defaults.

  # TODO: make this available to users outside of the
  # HetEffects() scope.
  quintilize = function(x){
    qs = quantile(x, probs = c(0, .25, .5, .75, 1))
    if(length(unique(qs)) < 5) warning("Binning quantitative variable to < 5 bins")
    qs[1] = qs[1] - 0.0001 # ensures cut doesn't set values at lowest quintile to NA
    cut(x, unique(qs))
  }

  # Generate matrix for prediction:
  Xmat_buckets = Xmat
  Xmat_buckets[col_numeric] = apply(Xmat_buckets[col_numeric], 2,
      quintilize)

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

  # TODO: Decide on what to do about NAs
  #  anyNA = function(x) any(is.na(x))
  #  apply(Xmat_buckets, 2, anyNA)


  # Put everything in mod matrix first:
  Xmatb <- data.frame(model.matrix(~ -1 + ., Xmat_buckets[,-1]))
  newXb <- data.frame(model.matrix(~ -1 +., newX))
#  names(newXb) == names(Xmatb)

  if(!bootstrap){

    bsl = SuperLearner(Y = Y,
            X = Xmatb, # don't pass intercept
            newX = newXb,
            family = family,
            SL.library = SL.library,
            method = method,
            id = id,
            verbose = verbose,
            control = control,
            cvControl = cvControl,
            obsWeights = weights)

    res = list(
      weights = bsl$coef,
      predictions = bsl$library.predict,
      effectX = newXb
    )

  } else {
    # bootstrap SuperLearner
    # get effects for every combination of covariates
    # install.packages('doMC')
    # require('doMC')  # for easy parallelization on a nice multicore machine
    # registerDoMC(cores = cores)
    # bootres = foreach(i=1:R, .combine=rbind, .errorhandling='remove') %dopar% {
    res = list()
    for( i in 1:R){
      d = sample(x = 1:length(Xmatb[,1]),
                 size = length(Xmatb[,1]),
                 replace=TRUE)

      bsl = SuperLearner(Y = Y[d],
              X = Xmatb[d,],
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

    }
  }
  return(list(boostrap_samples = res, effectX = newXb))
}

# TODO: process and format results, display for the user.
