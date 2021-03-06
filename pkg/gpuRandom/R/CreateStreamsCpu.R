#' @title Create streams
#'
#' @description Create streams in R.
#' 
#' @param x A vector indicating the number of work-items.
#' @return A R-streams.
#' @useDynLib gpuRandom
#' @export




CreateStreamsCpu = function(Nstreams) {

  result = matrix(as.integer(0), nrow=Nstreams, ncol=18)

  colnames(result) = c("current.g1.1", "current.g1.2", "current.g1.3", "current.g2.1", "current.g2.2", "current.g2.3",
                       "initial.g1.1", "initial.g1.2", "initial.g1.3", "initial.g2.1", "initial.g2.2", "initial.g2.3",
                       "substream.g1.1", "substream.g1.2", "substream.g1.3", "substream.g2.1", "substream.g2.2", "substream.g2.3")


  cpp_mrg31k3pCreateStreams(result)


}


