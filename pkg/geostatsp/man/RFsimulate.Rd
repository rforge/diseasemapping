\name{RFsimulate}
\alias{RFsimulate}
\alias{modelRandomFields}
\alias{RFsimulate}
\alias{RFsimulate-methods}
\alias{RFsimulate,ANY,Raster-method}
\alias{RFsimulate,numeric,SpatialGrid-method}
\alias{RFsimulate,numeric,SpatialPixels-method}
\alias{RFsimulate,numeric,SpatialPoints-method}
\alias{RFsimulate,numeric,GridTopology-method}
\alias{RFsimulate,RMmodel,GridTopology-method}
\alias{RFsimulate,RMmodel,SpatialPoints-method}
\alias{RFsimulate,matrix,Raster-method} 
\alias{RFsimulate,matrix,Spatial-method} 
\alias{RFsimulate,data.frame,ANY-method} 


\title{Simulation of Random Fields}

\description{
  This function simulates conditional and unconditional Gaussian random fields, calling the function in the RandomFields package of the same name.
}

\usage{
\S4method{RFsimulate}{ANY,Raster}(model, x,	data=NULL,
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{numeric,SpatialGrid}(model, x,data=NULL,
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{numeric,SpatialPixels}(model, x, data=NULL, 
err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{numeric,SpatialPoints}(model, x, data=NULL, 
 	err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{numeric,GridTopology}(model, x, data=NULL, 
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{RMmodel,GridTopology}(model, x, data=NULL, 
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{RMmodel,SpatialPoints}(model, x, data=NULL, 
 err.model=NULL, n = 1, ...)
\S4method{RFsimulate}{matrix,Raster}(model, x, 	data=NULL, 
 err.model=NULL, n = nrow(model), ...)
\S4method{RFsimulate}{matrix,Spatial}(model, x,	data=NULL, 
	err.model=NULL, n = nrow(model), ...)
\S4method{RFsimulate}{data.frame,ANY}(model, x,	data=NULL, 
	err.model=NULL, n = nrow(model), ...)
modelRandomFields(param, includeNugget=FALSE)
}

\arguments{
 \item{model}{object of class \code{\link[RandomFields]{RMmodel}}, a vector of named model
 parameters, or a matrix where each column is a model parameter
 }
 \item{x}{Object of type \code{\link[sp:GridTopology-class]{GridTopology}} or Raster
  or \command{\link[sp]{SpatialPoints}} or \command{\link[sp:SpatialPixels-class]{SpatialPixels}}.}
 \item{data}{For conditional simulation and random imputing only.
 If \code{data} is missing, unconditional
 simulation is performed.Object of class
 \command{\link[sp:SpatialPointsDataFrame-class]{SpatialPointsDataFrame}};
 coordinates and response values of
 measurements in case that conditional simulation is to
 be performed
 }
 \item{err.model}{For conditional simulation and random imputing only.\cr
 Usually \code{err.model=RMnugget(var=var)}, or not given at all
 (error-free measurements).
 }
 \item{n}{number of realizations to generate.
 }
 \item{...}{for advanced use:
 further options and control parameters for the simulation
 that are passed to and processed by \command{\link[RandomFields]{RFoptions}}
 }
 \item{param}{A vector of named parameters}
 \item{includeNugget}{If \code{FALSE}, the \code{nugget} parameter is ignored. }
}

\details{
If \code{model} is a matrix, a different set of parameters is used for each simulation.  If 
\code{data} has the same number of columns as \code{model} has rows, 
a different column \code{i} is used with parameters in row \code{i}.   
}

\value{
An object of the same class as \code{x}, with the exception of \code{x} being a \code{GridTopology}
where a Raster is returned.
 }


\author{Patrick E. Brown \email{patrick.brown@utoronto.ca}
}

\seealso{
 \command{\link[RandomFields]{RFsimulate}},
 \command{\link[RandomFields]{RFfit}},
 \command{\link[RandomFields]{RFgetModelInfo}},
 \command{\link[RandomFields]{RFgui}},
 \command{\link[RandomFields]{RMmodel}},
 \command{\link[RandomFields]{RFoptions}},
 \command{\link[RandomFields]{RFsimulateAdvanced}},
 \command{\link[RandomFields]{RFsimulate.more.examples}}
}

\examples{
library('geostatsp')
model <- c(var=5, range=1,shape=0.5)


myraster = raster(nrows=20,ncols=30,xmn=0,ymn=0,xmx=6,ymx=4, 
	crs=CRS("+proj=utm +zone=17 +datum=NAD27 +units=m +no_defs"))

set.seed(0) 

simu <- RFsimulate(model, x=myraster, n=3)

plot(simu[['sim2']])

# conditional simulation
firstSample = RFsimulate(
	model, 
	x=SpatialPoints(myraster)[seq(1,ncell(myraster), len=100), ],
	n=3
)

secondSample = RFsimulate(
	model = cbind(var=5:3, range=1:3, shape=seq(0.5, 1.5, len=3)),
	x= myraster,
	data=firstSample,n=4
)

plot(secondSample)

# convert the model to RandomFields format and plot
if(requireNamespace('RandomFields', quietly=TRUE)) {
	RandomFields::plot(modelRandomFields(model))
}


}

\keyword{spatial}

