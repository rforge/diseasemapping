setGeneric('glgm', 
		function(
				formula, data, grid, 
				covariates, 
				...) {
			standardGeneric("glgm")
		}
)

 # sort out formula
# null formula
setMethod("glgm", 
		signature("NULL"), 
		gm.nullFormula
		)


setMethod("glgm", 
		signature("numeric"),  
		gm.numericFormula
		)

# change character to formula
setMethod("glgm", 
		signature("character"),  
	gm.characterFormula
		)
		

# numeric cells, create raster from data bounding box

setMethod("glgm", 
		signature("formula", "ANY", "numeric"),
		gm.gridNumeric
		)



# extrat covariates for data, convert covariates to a stack
setMethod("glgm", 
		signature("formula", "Raster", "Raster"),
		gm.dataRaster
		)


		setMethod("glgm", 
				signature("formula", "Spatial", "Raster", "NULL"),
				gm.dataSpatial
		)

setMethod("glgm", 
		signature("formula", "Spatial", "Raster", "list"),
		gm.dataSpatial
		)

setMethod("glgm", 
				signature("formula", "Spatial", "Raster", "Raster"),
				gm.dataSpatial
		)
		
setMethod("glgm", 
		
				signature("formula", "Spatial", "Raster", "data.frame"),
				function(formula, data, grid, covariates, ...) {

		dataDF = data@data

		callGeneric(
			formula=formula, data=dataDF,
			grid = grid,
			covariates=covariates, 
			...
	)
}
	)
		

#################
#### the real work
##################


setMethod("glgm", 
		signature("formula", "data.frame", "Raster", "data.frame"), 
		function(formula, data,  grid, 
				covariates=data.frame(), 
				shape=1, priorCI=NULL, 
				mesh=FALSE,...) {

			
			
			if(!any(names(grid)=='space'))
				warning("grid must have a layer called space with inla cell ID's")


			if(!all(all.vars(formula)%in% names(data)))
				warning("some covariates seem to be missing: formula ", paste(all.vars(formula), collapse=" "), ", data: ", paste(names(data), collapse=" "))
			
			cells = trim(grid[['space']])
			firstCell = values(cells)[1]
			cellDim = dim(cells)[1:2]
			# first cell = 2 * buffer^2 + ncolSmall * buffer + buffer
			# buffer = -(nrowSmall+1) + sqrt (  (nrowSmall+1)^2 + 8 firstCell / 4
			buffer = (-(cellDim[1]+1) + sqrt(  (cellDim[1]+1)^2 + 8* (firstCell-1) ))/4

			
			
	# data, cells, and covariates must have varilable called 'space'		
	# values of cells must be index numbers, and cells shouldnt include the buffer		
	thedots = list(...)
			
	# priors for spatial standard deviation and nugget std dev.
	sdNames = unique(c("sd",grep("^sd", names(priorCI), value=TRUE)))
	# if model is Gaussian, look for prior for sdNugget
	if(!any(names(thedots)=="family")) {
		thedots$family =  "gaussian"
	}
	if(thedots$family=="gaussian") {
		sdNames = unique(c(sdNames, "sdNugget"))
	}
	
	precPrior=list()
	for(Dsd in sdNames) {
			
		if(any(names(priorCI)==Dsd)) {
		obj1 = sort(priorCI[[Dsd]]^-2)
		cifun = function(pars) {
				theci = 	pgamma(obj1, shape=pars[1], 
						rate=pars[2],log.p=T)
				
				(log(0.025) - theci[1])^2 +
				(2*(log(0.975) - theci[2]))^2		

			}
	
			
			
		precPrior2=optim(c(.5,.5/mean(obj1)), cifun, 
				lower=c(0.000001,0.0000001),method="L-BFGS-B")
		names(precPrior2$par) = c("shape","rate")
		precPrior[[Dsd]] = precPrior2$par 
				
 		#pgamma(obj1, shape= precPrior["shape"], rate=precPrior["rate"],log.p=F)
		#pgamma(obj1, shape= precPrior["shape"], rate=precPrior["rate"],log.p=T)
		#log(c(0.025, 0.975))
		#precPrior2
		#pgamma(obj1, shape=precPrior["shape"], rate=precPrior["rate"],log.p=T)
		#log(c(0.025, 0.975)) 
 		#1/sqrt(qgamma(c(0.975,0.025), shape=precPrior["shape"], rate=precPrior["rate"]))
		#priorCI$sd
		
		} else {
			precPrior[[Dsd]] = c(shape=0.01, rate=0.01)
		}
	}
		
	if("range" %in% names(priorCI)) {
		if(priorCI$range[1] < xres(cells)/4) {
			priorCI$range[1] = xres(cells)/4
			warning("lower bound of range CI too small, setting it to 1/4 cell size")
			
		}
		
		# rang parameter, in terms of cells, not km.
		obj1=sort(priorCI$range/xres(cells))
		
		cifun = function(pars) {
			theci = 		pgamma(obj1, shape=pars[1], rate=pars[2], log.p=T)
			
			(theci[1] - log(0.025))^2 +
					(theci[2] - log(0.925))^2 
		}
			
		ratePrior2=optim(c(2,2/mean(obj1)), cifun, 
				lower=c(0.001,0.001),method="L-BFGS-B")
		ratePrior = ratePrior2$par
		names(ratePrior ) = c("shape","rate")
if(FALSE) {
		ratePrior
		pgamma(obj1, shape=ratePrior["shape"], rate=ratePrior["rate"])
		 qgamma(c(0.975,0.025), shape=ratePrior["shape"], rate=ratePrior["rate"])
		xres(cells)*qgamma(c(0.025,0.975), shape=ratePrior["shape"], rate=ratePrior["rate"])
	}
	} else {
		ratePrior = c(shape=0.01, rate=0.01)
	}
	spaceFormula = paste(".~.+ f(space, model='matern2d', ",
				"nrow=", nrow(cells)+2*buffer, 
				", ncol=", ncol(cells)+2*buffer,
				", nu=", shape, 
				", hyper = list(",
				 "range=list( param=c(",
				      paste(ratePrior, collapse=","),
				"), prior='loggamma'),",
				"prec=list( param=c(",
				paste(precPrior$sd, collapse=","),
				"),prior='loggamma')",
				" ) )" 
			)
	
	formula = update.formula(formula,	as.formula(spaceFormula))


	
	# sort out factors
	thevars = rownames(attributes(terms(formula))$factors)
	thevars = grep("^factor\\(", thevars, value=TRUE)
	varsInData = apply(data, 2, is.factor)
	varsInData = names(data)[varsInData]
	thevars = c(varsInData, thevars)
	
	if(length(thevars)){
		thevars = gsub("^factor\\(|\\)", "", thevars)
		# loop through factors
		for(D in thevars){
			# biggest category is baseline
			thetable = table(data[,D])
			thebase = names(sort(thetable,decreasing=TRUE))[1]
			newLevels = unique(c(thebase, levels(factor(data[,D]))))
			data[,D] = factor(data[,D], levels=newLevels)
			covariates[,D] = factor(covariates[,D],
					levels=levels(data[,D]))
			
		}
	}
	theFactors = thevars

	
	# create linear combinations object for prediction.
	# create formula, strip out left variable and f(...) terms
 	formulaForLincombs = unlist(strsplit(as.character(formula), "~"))
	formulaForLincombs = formulaForLincombs[length(formulaForLincombs)]
	formulaForLincombs =
			gsub("\\+?[[:space:]]*f\\([[:print:]]*\\)[[:space:]]?($|\\+)", "+", formulaForLincombs)
	# strip out offsets
formulaForLincombs =
		gsub("\\+?[[:space:]]*offset\\([[:print:]]*\\)[[:space:]]?($|\\+)", "+", formulaForLincombs)

	# convert multiple + to a single +
	formulaForLincombs = gsub(
			"\\+[[:space:]]?\\+([[:space:]]?\\+)?", "+",
			formulaForLincombs)
	# strip out trailing +
formulaForLincombs = gsub("\\+[[:space:]]?$", "", formulaForLincombs)

	# if we have covariates
	if(nchar(formulaForLincombs)) {

		formulaForLincombs=as.formula(paste("~", formulaForLincombs))
	
	
		# variables in the model but not in prediction rasters
		thevars = rownames(attributes(terms(formulaForLincombs))$factors)
		thevars = gsub("^factor\\(|\\)", "", thevars)
		varsInPredict = thevars[thevars %in% names(covariates)]
		cantPredict = thevars[! thevars %in% names(covariates)]
		theFactors2 = grep("factor\\([[:print:]]*\\)", cantPredict)
		if(length(theFactors2)) {
			temp = cantPredict
			cantPredict = cantPredict[-theFactors2]
			theFactorsInFormula = temp[theFactors2]
		}
		if(length(cantPredict))
			covariates[,cantPredict]= 0

		if(FALSE){
		# check for factors	
 
		facData = apply(data, 2, is.factor)
		facData = names(facData[which(facData)])
		facCovariates = apply(covariates, 2, is.factor)
		facCovariates = names(facCovariates[which(facCovariates)])
		theFactors = unique(c(facData, facCovariates))
		theFactors = match(theFactors, names(covariates))
		
		# convert raster data into factors using same levels as in points data
		for(D in theFactors) {
			# convert columns in lincombs to factors
			# get rid of levels not present in the observed data

		covariates[,D] = factor(covariates[,D],
					levels=levels(data[,D]))
		}
	}
		covariates = covariates[,c("space", varsInPredict),drop=FALSE]
		lincombMat = model.matrix(update.formula(
						formulaForLincombs, ~.+space),
						covariates, na.action=NULL)
	
	} else { # no covariates
		lincombMat = cbind("(Intercept)"=rep(1, ncell(cells)), space=values(cells))
	# for some reason can't give name (Intercept) in the data.frame call.
	}

	
	
	thelincombs=list()	
	lincombMat[lincombMat==0] = NA

	for(D in 1:nrow(lincombMat)) {
		thisrow = lincombMat[D, ]
		thisrow = as.list(thisrow[!is.na(thisrow)  ])
		
		thelincombs[[D]] = 
			do.call(inla.make.lincomb, thisrow)$lc
		spaceCol = which(unlist(lapply(thelincombs[[D]],names))=='space')

		if(length(spaceCol))
			thelincombs[[D]][[spaceCol]]$space = list(
				weight=1, 
				idx=thelincombs[[D]][[spaceCol]]$space$weight)
}


names(thelincombs) = paste("c", lincombMat[,"space"],sep="")


	# get rid of observations with NA's in covariates
	allVars = all.vars(formula)
	theNA = apply(data[,allVars], 1, function(qq) any(is.na(qq)))
 
	data = data[!theNA,]
	if(any(names(thedots)=='Ntrials'))
		thedots$Ntrials = thedots$Ntrials[!theNA]

	forInla = thedots
	forInla$lincomb = c(thelincombs, forInla$lincomb)
	forInla$data = data
	forInla$formula = formula
	
	

		
	# if model is gaussian, add prior for nugget
	if(!is.null(precPrior$sdNugget)) {
		forInla$control.family$hyper$prec =
				list(prior="loggamma",
						param=precPrior$sdNugget
				) 
	}
	
	
	# get rid of some elements of forInla that aren't required
	forInla = forInla[grep("^buffer$", names(forInla), invert=TRUE)]

#	return(forInla)
 
	inlaResult = do.call(inla, forInla) 
 
	
	if(all(names(inlaResult)=="logfile"))
		return(c(forInla, inlares=inlaResult))
	

 	# parameter priors for result


	params = list(
			range = list(userPriorCI=priorCI$range, 
					priorCI = 
							xres(cells)*
							qgamma(c(0.025,0.975), 
									shape=ratePrior["shape"], 
									rate=ratePrior["rate"]),
					priorCIcells = 
							qgamma(c(0.975,0.025), 
									shape=ratePrior["shape"], 
									rate=ratePrior["rate"]),
					params.intern = ratePrior))
	rangeLim = 	qgamma(c(0.001,0.999), 
			shape=ratePrior["shape"], 
			rate=ratePrior["rate"])
	rangeLim = xres(cells) *rangeLim
	rangeSeq = seq(min(rangeLim), max(rangeLim), len=1000)
	rangeSeqCells = rangeSeq/xres(cells)
	params$range$prior=cbind(
			x=rangeSeq,
			y=dgamma(rangeSeqCells, shape=ratePrior["shape"], 
					rate=ratePrior["rate"])  / xres(cells)
	)
	if(FALSE) {
		plot(params$range$prior, type='l')
		
		sum(params$range$prior[,"y"])*range(diff(params$range$prior[,"x"]))
	}	
	

	
	for(Dsd in names(precPrior)) {
		params[[Dsd]] = list(userPriorCI=priorCI[[Dsd]], 
			priorCI = 1/sqrt(
				qgamma(c(0.975,0.025), 
						shape=precPrior[[Dsd]]["shape"], 
						rate=precPrior[[Dsd]]["rate"])),
					params.intern=precPrior[[Dsd]])
	
	precLim = 	qgamma(c(0.999,0.001), 
			shape=precPrior[[Dsd]]["shape"], 
			rate=precPrior[[Dsd]]["rate"])
	sdLim = 1/sqrt(precLim)
	sdSeq = seq(min(sdLim), max(sdLim), len=1000)
	precSeq = sdSeq^(-2)
	params[[Dsd]]$prior=cbind(
			x=sdSeq,
			y=dgamma(precSeq, shape=precPrior[[Dsd]]["shape"], 
					rate=precPrior[[Dsd]]["rate"]) *2* (precSeq)^(3/2) 
	)

	}

	
	# random into raster
# E exp(random)

if("summary.random" %in% names(inlaResult)) {

temp=unlist(
		lapply(inlaResult$marginals.random$space, function(qq) {
					sum(
							exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
					)
				})
)
inlaResult$summary.random[['space']][,"exp"] = temp



		forRast = 	as.matrix(inlaResult$summary.random[["space"]][values(cells),])
		resRasterRandom = 
				brick(extent(cells), nrows=nrow(cells),
						ncols=ncol(cells), crs=projection(cells),
						nl=dim(forRast)[2])
		names(resRasterRandom) = 
				paste("random.", colnames(forRast),sep="")
		
		values(resRasterRandom) = as.vector(forRast)
		
	} else {
		return(list(inla=inlaResult, parameters=params))
	}

	inlaResult$marginals.random$space = inlaResult$marginals.random$space[values(cells)]
	
	# E exp(lincombs)
	temp=unlist(
			lapply(inlaResult$marginals.lincomb.derived, function(qq) {
						sum(
								exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
						)
					})
	)
	inlaResult$summary.lincomb.derived[,"exp"] = temp
	
	# E inv logit(lincombs)
	if(length(grep("binomial",inlaResult$.args$family))) {
		temp=unlist(
				lapply(inlaResult$marginals.lincomb.derived, function(qq) {
							eqqx = exp(qq[,"x"])
							sum(
								eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
							)
						})
		)
		inlaResult$summary.lincomb.derived[,"invlogit"] = temp		
	}
 

	# lincombs into raster
theSpaceName = grep("^c[[:digit:]]+$", names(inlaResult$marginals.lincomb.derived), value=TRUE)
theSpace = as.integer(gsub("^c", "", theSpaceName))

	linc = inlaResult$summary.lincomb.derived[theSpaceName,]
	linc$space = theSpace

	missingCells = values(cells)[! values(cells) %in% theSpace]

	if(length(missingCells)) {
		toadd = matrix(NA, length(missingCells), dim(linc)[2], 
					dimnames=list(
							paste("c", missingCells, sep=""), 
							colnames(linc)
			)
			)
		toadd[,"space"] = missingCells
	
		linc = rbind(linc, toadd)

		missingNames = rownames(toadd)
		missingMarginals = vector("list", length(missingNames))
		names(missingMarginals) = missingNames
		
		inlaResult$marginals.lincomb.derived = c(	
				inlaResult$marginals.lincomb.derived,
				missingMarginals)
		
	}
	linc = as.matrix(linc[paste("c", values(cells), sep=""),])

	resRasterFitted = 
			brick(extent(cells), nrows=nrow(cells),
					ncols=ncol(cells), crs=projection(cells),
					nl=ncol(linc))
	names(resRasterFitted) = 
			paste("predict.", colnames(linc),sep="")
	
	values(resRasterFitted) = as.vector(linc)
	
	
	# Add in empty lists for the marginals of missing cells
	
	inlaResult$marginals.predict = 
			inlaResult$marginals.lincomb.derived[
					paste("c", values(cells), sep="")
					]
			
	
	
	# posterior distributions




params$range$posterior=inlaResult$marginals.hyperpar[["Range for space"]]
params$range$posterior[,"x"] =  xres(cells) * params$range$posterior[,"x"]
params$range$posterior[,"y"] = params$range$posterior[,"y"] / xres(cells)


# sum(c(0,diff(params$range$posterior[,"x"])) * params$range$posterior[,"y"])
# sum(c(0,diff(params$range$prior[,"x"])) * params$range$prior[,"y"])


params$summary = inlaResult$summary.fixed

params$summary = cbind(params$summary, 
		meanExp = unlist(
				lapply(inlaResult$marginals.fixed,
						function(qq) {
							sum(
									exp(qq[,"x"])*c(0,diff(qq[,"x"]))*qq[,"y"]	
							)
						}
				))
)

if(length(grep("binomial",inlaResult$.args$family))) {
	params$summary = cbind(params$summary, 
			meanInvLogit = unlist(
					lapply(inlaResult$marginals.fixed, function(qq) {
								eqqx = exp(qq[,"x"])
								sum(
										eqqx/(1+eqqx)*c(0,diff(qq[,"x"]))*qq[,"y"]	
								)
							}
							)
	))
}


thecols = paste(c("0.975", "0.5","0.025"), "quant", sep="")

thesd = c(
		sdNugget= grep("^Precision[[:print:]]*Gaussian observations$", 
				names(inlaResult$marginals.hyperpar), value=TRUE),
		sd = grep("^Precision[[:print:]]*space$", 
				names(inlaResult$marginals.hyperpar), value=TRUE)
)

params$summary = rbind(params$summary,
		matrix(NA, nrow=length(thesd)+1, ncol=ncol(params$summary),
				dimnames = list(c("range", names(thesd)), NULL))
)


# convert precisions to standard deviations
for(Dsd in names(thesd)) {
	
	params[[Dsd]]$posterior=
			inlaResult$marginals.hyperpar[[thesd[Dsd]]]
	params[[Dsd]]$posterior[,"y"] = params[[Dsd]]$posterior[,"y"] * 2*  
			params[[Dsd]]$posterior[,"x"]^(3/2) 
	params[[Dsd]]$posterior[,"x"] = 1/sqrt(params[[Dsd]]$posterior[,"x"])  
	params[[Dsd]]$posterior = params[[Dsd]]$posterior[
			seq(dim(params[[Dsd]]$posterior)[1],1),]		

	params$summary[Dsd, thecols] = 
				1/sqrt(inlaResult$summary.hyperpar[
								thesd[Dsd],rev(thecols)])

		
	params$summary[Dsd,"mean"] =sum(
		1/sqrt(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"])*
			c(0,diff(inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"x"]))*
				inlaResult$marginals.hyperpar[[thesd[Dsd]]][,"y"]
	)
}

# put range in summary, in units of distance, not numbers of cells
thecolsFull =c("mean","sd",thecols,"mode") 
params$summary["range",thecolsFull]=				
		xres(cells)*
		inlaResult$summary.hyperpar[
				"Range for space",
				thecolsFull
		]
dimnames(params$summary) = lapply(dimnames(params$summary),
		function(qq) {
			qq=gsub("_", "\\\\textunderscore~", qq)
			qq=gsub("\\$", "\\\\textdollar~", qq)
			qq=gsub("<", "\\\\textless~", qq)
			qq=gsub(">", "\\\\textgreater~", qq)
			qq
		}
)
params$summary = as.data.frame(params$summary)

for(Dvar in names(covariates)) {
	theLevels =levels(covariates[[Dvar]])[[1]]
	if(!is.null(nrow(theLevels))){
	for(D in 1:nrow(theLevels)) {
		rownames(params$summary) = gsub(
			paste("(factor)?(\\()?", Dvar, "(\\))?:?", 
					theLevels[D,1],"$",sep=""),
			paste(Dvar, ":",theLevels[D,2],sep=""), 
					rownames(params$summary))
	}
}
}



resRaster=stack(resRasterRandom, resRasterFitted, cells)

#	if(!is.null(smallBbox))
#		resRaster = crop(resRaster, extent(smallBbox))
		
	result=list(inla=inlaResult,
					raster=resRaster,
					parameters=params
			)

	result
	
}
)