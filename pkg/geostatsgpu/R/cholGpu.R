
cholGpu = function(x,D, 
	control = list(
		workgroupSize = NULL, 
		localWorkgroupSize = NULL, 
		localStorage = NULL, 
		sizeOfDouble = NULL), 
	verbose=FALSE) {

	file <- system.file("CL", "cholGpu.cl", package = "geostatsgpu")
	kernel <- readChar(file, file.info(file)$size)

	if(missing(D)) {
		D = vclVector(0, nrow(x),
			type='double',
			ctx_id = x@.context_index)
	}


	Scontrol = c(
		'workgroupSize', 'localWorkgroupSize', 
		'localStorage',
		'sizeOfDouble')

	missingControl = setdiff(Scontrol, 
		names(control))

	if(length(missingControl) | 
		any(unlist(lapply(control, is.null)))) {
	
		localInfo = gpuNlocal(
			kernel, 
			'cholOffDiag',
			x@.device_index)

		for(Dcontrol in c(
			'localWorkgroupSize', 
			'sizeOfDouble')) {

			if(!length(control[[Dcontrol]])) {
				control[[Dcontrol]] = localInfo[[Dcontrol]]
			}
		}

		if(!length(control$workgroupSize)) {
			control$workgroupSize = localInfo$maxWorkgroupSize
		}

		if(!length(control$localStorage)) {
			control$localStorage = pmax(
				control$localWorkgroupSize,
				floor(
				0.95*localInfo$localMemory / 
				control$sizeofDouble) ) 
		}
	}

	if(verbose){
		message(paste('global work items', 
			control$workgroupSize, 
			'local work items', 
			control$localWorkgroupSize, 
			'local storage', 
			control$localStorage, 
			"N", ncol(x)))
	}


	diagTimesRowOfA = gpuR::vclVector(0, nrow(x),
            type='double',
            ctx_id = x@.context_index)

	diagWorking = gpuR::vclVector(0, 
		1 + round(
			control$workgroupSize / 
			control$localWorkgroupSize),
            type='double',
            ctx_id = x@.context_index)

	fromC = cpp_cholGpu(
		x, 
		D,
		diagWorking,
		diagTimesRowOfA,
		control$workgroupSize, 
		control$localWorkgroupSize, 
		control$localStorage,
		1,1,# NEED TO FIX colGroupwise, Ncrossprod
		TRUE, # verbose
		kernel
		)

	res = list(L = x, D=D, logDet = fromC, 
		control = control)

	res
}