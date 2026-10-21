#' Prepare the current reduced-mode terrain and surface-memory response
#'
#' Low-level metric-grid interface to the manuscript formulation. Grid rows run
#' south to north. Source roughness is an explicit effective input, not inferred
#' from target land cover. This interface does not download or reproject data.
#' Heights other than 10 m are supported numerically, not validated universally.
#' @param dem,roughness Finite numeric matrices on one square-cell metric grid.
#'   Roughness is aerodynamic length in metres, not land-cover class codes.
#' @param dx Cell spacing in metres.
#' @param x0,y0 Coordinates of the southwest cell centre.
#' @param x,y Target coordinates, in the same metric coordinate system.
#' @param depth Positive neutral reference depth, in metres; supplied explicitly.
#' @param source_z0 Effective source-profile roughness, in metres.
#' @param height Evaluation height above the local surface, in metres.
#' @param profile_z0 Terrain-profile reference roughness, in metres.
#' @param min_scale Smallest vertical mode scale, in metres.
#' @param reflection Reflect terrain by one domain on each side.
#' @param max_cells Maximum cells after reflection, checked before allocation.
#' @return A prepared response with a 2-by-2 terrain matrix, 72 directional
#'   surface factors and solver diagnostics. Construction is single-threaded.
#' @export
taruwind_modal_prepare <- function(dem, roughness, dx, x0, y0, x, y,
    depth, source_z0, height=10, profile_z0=.03, min_scale=10,
    reflection=TRUE, max_cells=2e6) {
  if (!is.matrix(dem) || !is.numeric(dem) || min(dim(dem))<3 ||
      any(!is.finite(dem))) stop("Finite numeric DEM matrix, at least 3 by 3, required")
  if (!is.matrix(roughness) || !is.numeric(roughness) ||
      !identical(dim(dem),dim(roughness))) stop("Roughness must match DEM dimensions")
  inputs <- list(dx,x0,y0,x,y,depth,source_z0,height,profile_z0,min_scale)
  if (any(!vapply(inputs,function(z) is.numeric(z)&&length(z)==1L&&is.finite(z),TRUE)))
    stop("Finite numeric scalar geometry and profile inputs required")
  if(dx<=0 || min_scale<=0 || depth<=min_scale || height<=0 || height>100 ||
     profile_z0<=0 || profile_z0>min_scale || source_z0<=0 || source_z0>=height ||
     any(!is.finite(roughness)) || any(roughness<=0 | roughness>=height))
    stop("Invalid spacing, scales or roughness/profile bounds")
  if(!is.logical(reflection)||length(reflection)!=1L||is.na(reflection)) stop("Logical reflection required")
  if(!is.numeric(max_cells)||length(max_cells)!=1L||!is.finite(max_cells)||max_cells<9)
    stop("Positive finite max_cells >= 9 required")
  if(length(dem)*(if(reflection) 9 else 1)>max_cells)
    stop("Terrain exceeds max_cells after reflection; explicitly raise the limit only with sufficient memory")
  row <- (y-y0)/dx; col <- (x-x0)/dx
  if(row<0 || col<0 || row>nrow(dem)-1 || col>ncol(dem)-1) stop("Target outside cell centres")
  h <- dem
  if(reflection) {
    ix <- function(n) {p<-((-n):(2*n-1))%%(2*(n-1));ifelse(p<n,p,2*(n-1)-p)+1L}
    h <- dem[ix(nrow(dem)),ix(ncol(dem)),drop=FALSE]
    row <- row+nrow(dem);col <- col+ncol(dem)
  }
  scales <- exp(seq(log(min_scale),log(depth),length.out=6L))
  cc <- modal_coefficients(h,dx,scales,profile_z0,16*min_scale)
  args <- c(list(h=as.vector(t(h)),ny=nrow(h),nx=ncol(h),dx=dx),cc,
    list(penalty=16*min_scale,height=height,z0=profile_z0,reference=100,
         row=row,col=col,rtol=1e-6,maxiter=600L,
         targets=matrix(c(row,col,height),nrow=1)))
  fields <- lapply(1:2,function(i) do.call(modal_native,
    c(args,list(u=as.numeric(i==1),v=as.numeric(i==2)))))
  residual <- vapply(fields,function(z) z$relative_residual,0.)
  if(any(!is.finite(residual)) || any(residual>1.01e-6))
    stop("Terrain solve did not meet the true-residual tolerance")
  profile <- log1p(height/profile_z0)/log1p(100/profile_z0)
  mat <- cbind(fields[[1]]$targets[1,1:2],fields[[2]]$targets[1,1:2])/profile
  rose <- memory_rose_native(roughness,x0,y0,dx,dx,x,y,source_z0,height=height)
  structure(list(matrix=mat,surface=rose,
    residual=residual,iterations=vapply(fields,function(z) z$iterations,0.),
    provenance=list(formulation="six-mode/exact-fetch-memory",height=height,
      depth=depth,source_z0=source_z0,profile_z0=profile_z0,min_scale=min_scale,
      penalty=16*min_scale,reflection=reflection,smoothing_m=0,threads=1L,
      observations_used=FALSE,package_version="0.1.1")),class="taruwind_modal_response")
}

#' Apply a current reduced-mode response to target-interpolated wind
#' @param object Response from [taruwind_modal_prepare()].
#' @param u,v Matching numeric component vectors in the grid coordinate frame,
#'   already vector-interpolated to the target at the prepared height.
#' @return A tibble of wind speeds in metres per second.
#' @export
taruwind_modal_predict <- function(object,u,v) {
  if(!inherits(object,"taruwind_modal_response")) stop("Prepared modal response required")
  if(!is.numeric(u)||!is.numeric(v)||length(u)!=length(v)) stop("Matching numeric u/v required")
  # Repeat a target vector at all four nodes: weights select it unchanged.
  speed <- wind_apply_native(matrix(rep(u,4L),ncol=4L),
    matrix(rep(v,4L),ncol=4L),rep(.25,4L),as.vector(t(object$matrix)),
    object$surface$factor,2L)
  tibble::tibble(speed=speed)
}
