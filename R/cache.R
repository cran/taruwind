#' Inspect downloaded cache files without changing them
#'
#' Reports files and their total size. Files are not removed, opened as data,
#' or uploaded. A caller-selected persistent cache survives R restarts;
#' the default cache is inside the current session's temporary directory.
#' @param cache_dir Cache root passed to [taruwind_downscale()].
#' @return A tibble with relative path, bytes and modification time. Attributes
#'   `cache_dir` and `total_bytes` identify the root and total file size.
#' @export
#' @examples
#' info = taruwind_cache_info()
#' attr(info, "total_bytes")
taruwind_cache_info = function(cache_dir=file.path(tempdir(), "taruwind-cache")) {
  if (!is.character(cache_dir) || length(cache_dir)!=1L || is.na(cache_dir) || !nzchar(cache_dir))
    stop("cache_dir must be a non-empty path")
  paths = if (dir.exists(cache_dir)) list.files(cache_dir,recursive=TRUE,
    all.files=TRUE,full.names=TRUE,no..=TRUE) else character()
  details = file.info(paths)
  keep = !is.na(details$isdir) & !details$isdir
  paths = paths[keep]; details = details[keep,,drop=FALSE]
  root = normalizePath(cache_dir,mustWork=FALSE)
  out = tibble::tibble(file=substring(paths,nchar(sub("/+$","",cache_dir))+2L),
    bytes=as.numeric(details$size),modified=details$mtime)
  attr(out,"cache_dir") = root
  attr(out,"total_bytes") = sum(out$bytes)
  out
}
