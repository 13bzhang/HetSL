SL.krls<- function(Y, X, newX, whichkernel='gaussian', ...){
	require('KRLS')
	fit.krls<- krls(X = X, y = Y, whichkernel= 'gaussian')
	pred<- predict(fit.krls, newdata = newX)$fit
	
	fit<- list(object = fit.krls)
	out<- list(pred = pred, fit = fit)
	class(out$fit)<- c('SL.krls')
	return(out)
	}
	
	
	predict.SL.krls<- function(object, newdata,...){
		require('KRLS')
		pred<- predict(object, newdata)$fit
		return(pred)
		}
		