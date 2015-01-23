SL.findit<- function(Y, X, newX, ...){
	require(FindIt)
	treatments <- attr(Xmat, "treatments")
 	Treat_columns = c()
	for(t in treatments){
    	Treat_columns <- c(Treat_columns, grep(paste('^', treatments[i],sep=''),  colnames(Xmat))) 
	  }
  	Treat_columns = Treat_columns[!duplicated(Treat_columns)]
  	attr(Xmat, "Treat_columns") <- Treat_columns
	
	mkstand<- function(x){
	x<- x - mean(x, na.rm=T)
	if(sd(x, na.rm=T)>0){
		x<- x/sd(x, na.rm=T)
		}
	return(x)
	}
	
	
	Y_find<- Y
	
	Y_find[which(Y==0)]<- -1	

	Xst<- 
	
	fit.findit<- 