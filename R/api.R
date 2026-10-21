#' Build a current TaruWind response
#'
#' This is the current six-mode solver, not the retired spectral/shallow-layer
#' interface. Supply the explicit metric-grid arguments documented below.
#' @param ... Arguments passed to [taruwind_modal_prepare()].
#' @return A prepared TaruWind response.
#' @export
#' @examples
#' h = matrix(0, 5, 5)
#' op = taruwind_build(h, h + .03, 90, 0, 0, 180, 180, 900, .03)
#' taruwind_predict(op, 3, 4)
taruwind_build <- function(...) taruwind_modal_prepare(...)

#' Predict wind speed from a prepared response
#' @param object A response from [taruwind_build()].
#' @param u,v Target-interpolated component vectors, or matching matrices with
#'   four columns for the four source nodes. Components must use the grid frame.
#' @param weights Four non-negative interpolation weights summing to one,
#'   required for matrix input. The same weights apply to every row.
#' @return A tibble with a speed column in metres per second.
#' @export
taruwind_predict <- function(object,u,v,weights=NULL) {
  if(!inherits(object,"taruwind_modal_response")) stop("Current prepared TaruWind response required")
  if(is.matrix(u)||is.matrix(v)) {
    if(!is.matrix(u)||!is.matrix(v)||!is.numeric(u)||!is.numeric(v)||
       !identical(dim(u),dim(v))||ncol(u)!=4L) stop("Matching numeric four-column u/v matrices required")
    if(!is.numeric(weights)||length(weights)!=4L||any(!is.finite(weights))||
       any(weights<0)||abs(sum(weights)-1)>1e-10) stop("Four non-negative weights summing to one required")
    return(tibble::tibble(speed=wind_apply_native(u,v,weights,
      as.vector(t(object$matrix)),object$surface$factor,2L)))
  }
  if(!is.null(weights)) stop("weights are only used with four-node matrix input")
  taruwind_modal_predict(object,u,v)
}

#' Resolve portable batch concurrency
#' @param workers "auto", or an integer from 1 to 16.
#' @param n_sites Number of independent requests.
#' @param cores Detected logical cores; injectable for reproducible policy tests.
#' @return An integer worker count. Auto uses half the cores, rounded down,
#'   bounded below by one and above by 16 and the request count.
#' @export
taruwind_workers <- function(workers="auto",n_sites=1L,cores=parallel::detectCores()) {
  if(!is.numeric(n_sites)||length(n_sites)!=1L||!is.finite(n_sites)||n_sites<1||n_sites!=floor(n_sites)) stop("Positive integer n_sites required")
  if(identical(workers,"auto")) {
    if(!is.numeric(cores)||length(cores)!=1L||!is.finite(cores)||cores<1) cores<-1
    workers<-min(16,max(1,floor(cores/2)))
  }
  if(!is.numeric(workers)||length(workers)!=1L||!is.finite(workers)||
     workers<1||workers>16||workers!=floor(workers)) stop("workers must be auto or an integer from 1 to 16")
  as.integer(min(workers,n_sites))
}

#' Prepare independent metric-grid sites with portable workers
#' @param requests Non-empty list of named argument lists for [taruwind_build()].
#' @param workers "auto" by default; see [taruwind_workers()].
#' @param max_total_cells Limit on padded cells concurrently in flight.
#'   This is an allocation guard, not a measurement of available RAM.
#' @return A list of prepared responses in input order.
#' @export
taruwind_build_many <- function(requests,workers="auto",max_total_cells=4e6) {
  if(!is.list(requests)||!length(requests)||!all(vapply(requests,is.list,TRUE))) stop("Non-empty list of request lists required")
  if(!is.numeric(max_total_cells)||length(max_total_cells)!=1L||!is.finite(max_total_cells)||max_total_cells<=0) stop("Positive max_total_cells required")
  cells<-vapply(requests,function(q) {
    if(!is.matrix(q$dem)) stop("Each request needs a named dem matrix")
    length(q$dem)*(if(isFALSE(q$reflection)) 1 else 9)
  },0.)
  w<-taruwind_workers(workers,length(requests))
  if(max(cells)>max_total_cells) stop("One request exceeds max_total_cells")
  safe<-max(1,floor(max_total_cells/max(cells)))
  if(w>safe) {warning("Worker count reduced by the padded-cell memory guard");w<-as.integer(safe)}
  if(w==1L) return(lapply(requests,function(q) do.call(taruwind_build,q)))
  cl<-parallel::makePSOCKcluster(w);on.exit(parallel::stopCluster(cl),add=TRUE)
  init<-function(paths) {.libPaths(paths);Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1);NULL}
  # Do not serialize this package's namespace before setting worker libraries:
  # that could load an older globally installed taruwind on the worker.
  environment(init)<-baseenv()
  parallel::clusterCall(cl,init,.libPaths())
  task<-function(q) do.call(taruwind::taruwind_build,q)
  environment(task)<-baseenv()
  parallel::parLapply(cl,requests,task)
}
