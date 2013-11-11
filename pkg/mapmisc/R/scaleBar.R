
scaleBar = function(crs, pos="bottomright",scale.cex=1,...) {
	
	if(is.character(crs))
		crs = CRS(crs)
	if(class(crs) != "CRS")
		crs = CRS(proj4string(crs))
	

#	dash = "\u2517\u2501\u2501\u2501\u2537\u2501\u2501\u2501\u251B"
	
	oldcex = par("cex")
	forLegend = list(...)
	if(length(forLegend$cex)){
		par(cex=forLegend$cex)
	} 
	forLegend$cex=1
	
	
	xpoints = t(bbox(extent(par("usr"))))
	
	dashTemplate = " 2000 km "
	xcentre = apply(xpoints, 2, mean)
	xpoints = rbind(centre=xcentre, 
			dashright = xcentre + c(strwidth(dashTemplate),0)
	)

	
	xpoints = SpatialPoints(xpoints, proj4string=crs)


	
	xll = spTransform(xpoints, CRS("+init=epsg:4326"))
	

	up = matrix(coordinates(xll)["centre",]+c(0,0.1),ncol=2,
			dimnames=list("up",NULL))
	
	
	xll=rbind(xll, SpatialPoints(up, 
					proj4string=CRS(proj4string(xll))))
	
	
	xpoints2 = spTransform(xll[c("up","centre")], CRS(proj4string(xpoints)))
	thediff=apply(coordinates(xpoints2), 2,diff)
	north=atan(thediff[1]/thediff[2])
	

	dashdist = spDists(xll[c("centre","dashright"),], 
			longlat=TRUE)[1,2]*1000 
	bardist = 	dashdist*scale.cex
		
	theb = log10(bardist)
	candidates = 10^c(floor(theb), ceiling(theb))
	candidates = c(candidates[1]*c(1,2,5), candidates[2])
	segdist=candidates[order(abs(candidates - bardist))[1]]
	
	segscale = ( strwidth(dashTemplate)/strwidth("m") ) *segdist / dashdist
	


	
	if(segdist >900) {
		lunits="km"
		segdist = segdist / 1000
	} else {
		lunits="m"
	}

	eps = 0.175
	
	dimIn = par("pin")
	dimUser = par("usr")
	dimUser = c(dimUser[2]-dimUser[1], dimUser[4]-dimUser[3])
	InPerUnit = dimIn/dimUser
	

	theN = c(0, 0+1i, eps+1i, (1-eps)+1.5*eps*1i,
			(1-eps)+1i, 1+1i, 1+0i,
			(1-eps)+0i, eps+(1-1.5*eps)*1i,
			eps) - 0.5 - 0.5*1i
	theHat = c(-0.25+1i, 0.5+1.6i, 1.25+1i)
	theHat = c(theHat, rev(theHat) + 1.5*eps*1i)- 0.5-0.5*1i
	theN =  strwidth("N")*Re(theN) + 
			1i*strwidth("N")*InPerUnit[1]/InPerUnit[2]*
			Im(theN)
	theHat =  strwidth("N") * Re(theHat)+ 
			1i*strwidth("N")*InPerUnit[1]/InPerUnit[2]*
			Im(theHat)
	
	theN = theN * exp(-1i*north)
	theHat = theHat * exp(-1i*north)
	
	
	
	thelabel = paste(segdist, lunits,sep="")
	

	
	defaults = list(col='black', 
			xjust=0.7, inset=0.001, bg="white",
			x=pos, text.width=strwidth("I"), pt.cex=1)

	for(D in names(defaults)) {
		if(is.null(forLegend[[D]]))
			forLegend[[D]] = defaults[[D]]			
	}

	defaults = list(text.col = forLegend$col,
			title.adj = 0.5*(segscale*strwidth("m")/
						(segscale*strwidth("m") + 
							2*strwidth("m")+
							forLegend$text.width)))
	
	for(D in names(defaults)) {
		if(is.null(forLegend[[D]]))
			forLegend[[D]] = defaults[[D]]			
	}
	
	forLegend$lty = 1
	forLegend$pch = NA
	forLegend$seg.len = segscale
	forLegend$title=thelabel
	forLegend$legend = NA
	forLegend$lwd=3

	
	if(forLegend$seg.len*strwidth("m") < strwidth(forLegend$title)) {
		forLegend$title=NA
	}

		
	thelegend = do.call(legend, forLegend)
			
	if(is.na(forLegend$title))
		text(thelegend$text$x - (2/3)*strwidth("m")*forLegend$seg.len,
				thelegend$rect$top , 
				label=thelabel, pos=1, cex=0.75, offset=1.25)
		
	thecentre =  thelegend$text$x + 1i*thelegend$text$y
	 polygon(forLegend$pt.cex*theN +thecentre, 
			 col=forLegend$text.col,border=NA)
	 polygon(forLegend$pt.cex*theHat + thecentre, 
			 col=forLegend$text.col,border=NA)
	 
	par(cex=oldcex)
	 
	 return(invisible(c(thelegend,forLegend))	)
}