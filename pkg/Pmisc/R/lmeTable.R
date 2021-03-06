#' @export
lmeTable = function(x) {
  fromLme = x
		thetable = summary(fromLme)$tTable
		toRep = dim(thetable)[2]-1
		if(length(fromLme$modelStruct$reStruct)==1) {
			newnames = c(rownames(thetable), c('$\\sigma$', '$\\tau$'))
			thetable = rbind(thetable, 
				c(sqrt(as.matrix(fromLme$modelStruct$reStruct[[1]]))* fromLme$sigma, rep(NA, toRep)),
				c(fromLme$sigma, rep(NA, toRep)) 
			)
		} else {
		newnames = c(rownames(thetable), c(names(fromLme$modelStruct$reStruct), '$\\tau$'))
		thetable = rbind(thetable, 
		cbind(sqrt(diag(as.matrix(fromLme$modelStruct$reStruct[[1]])))* fromLme$sigma, matrix(NA, length(fromLme$modelStruct$reStruct), toRep)),
		c(fromLme$sigma, rep(NA, toRep)) )
		}
		rownames(thetable) = newnames
		colnames(thetable)[colnames(thetable)=='Value'] = 'MLE'

		# if model is ar1
		if('corStruct' %in% names(summary(fromLme$modelStruct))){
			rangeNugget = c(
					summary(fromLme$modelStruct)$corStruct[[1]],
					summary(fromLme$modelStruct)$corStruct[[2]])
			rangeNugget = exp(rangeNugget)
			thetable = rbind(thetable, range=NA, sigmav = NA)
			tauRow = grep("tau\\$$", rownames(thetable))
			sigmaRow = grep("sigma\\$$", rownames(thetable))
			thetable['range','MLE'] = rangeNugget[1]
			tausq = rangeNugget[2] * thetable[tauRow, 'MLE']^2
  		thetable['sigmav','MLE'] = sqrt( 
					thetable[tauRow, 'MLE']^2 - tausq)	
			thetable[tauRow, 'MLE'] =  sqrt(tausq)
			rownames(thetable) = gsub("sigmav", "$\\\\sigma_V$",
					rownames(thetable))
			rownames(thetable)[sigmaRow] = '$\\sigma_U$'
		}
		thetable
}
		