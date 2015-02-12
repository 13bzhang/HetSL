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
#'        Data MUST NOT contain NAs.
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

  # Make sure there are no numeric variables
  col_numeric <- apply(Xmat, 2, function(x) length(unique(x)) > 20 )
  if (any(col_numeric)) stop(
    "Columns in design matrix must be discrete (20 unique values or fewer).")

  # Generate matrix for prediction:

  # First make sure there are not too many factor combinations
  veclist = apply(Xmat[,-1], 2, unique)
  veclist = lapply(veclist, na.omit)
  veclist = lapply(veclist, as.vector)
  if (prod(unlist(lapply(veclist, length))) > 1e10) {
    stop("Too many factor combinations to estimate effects")
  }

  newX = expand.grid(veclist, stringsAsFactors=FALSE)

  # TODO: Check to make sure there are no NAs
  #  anyNA = function(x) any(is.na(x))
  #  apply(Xmat, 2, anyNA)


  # Put everything in mod matrix first:
  Xmatb <- data.frame(model.matrix(~ -1 + ., Xmat[,-1]))
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


#' Split by Quantile
#'
#' The split_by_quantile function splits a quantitative variable into
#' n - 1 quantile chunks.
#'
#' @usage
#' split_by_quantile(x, n)
#'
#' @param x The variable
#' @param n The n-th quantiles to split the variable, results in n - 1 levels.
#'
#' @return A factor with n - 1 levels.
#' @export
#' @author Solomon Messing
split_by_quantile = function(x, n){
  qs = quantile(x, probs = seq(0, 1, length.out = n))
  if(length(unique(qs)) < n) warning("Binning quantitative variable to < n bins")
  qs[1] = qs[1] - 0.0000001 # ensures cut doesn't set values at lowest quintile to NA
  cut(x, unique(qs))
}

# Example usage:
# Xmat_buckets = Xmat
# Xmat_buckets[col_numeric] = apply(Xmat_buckets[col_numeric], 2,
#     quintilize)






# TODO: process and format results, display for the user.
