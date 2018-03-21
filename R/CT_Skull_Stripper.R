#' @title CT Skull Stripping Wrapper
#' @description Simple Wrapper for \code{\link{CT_Skull_Strip}} and
#' \code{\link{CT_Skull_Strip_robust}}, with simple flag for robustness
#' @param ... (character) File to be skull stripped or object of class
#' @param maskfile (character) Filename for mask (if \code{keepmask = TRUE}).
#' If \code{NULL}, then a temporary file will be generated.
#' nifti
#' @param robust (logical) output filename
#' @export
CT_Skull_Stripper <- function(
  ...,
  robust = FALSE
){

  L = list(...)
  nL = names(L)
  if (is.null(nL)) {
    nL = rep("", length = length(L))
  }
  if (robust) {
    n = names(formals(CT_Skull_Strip_robust))
    n = c("", n)
    ind = which(nL %in% n)
    L = L[ind]
    res = do.call("CT_Skull_Strip_robust", L)
  } else {
    n = names(formals(CT_Skull_Strip))
    n = c("", n)
    ind = which(nL %in% n)
    L = L[ind]
    res = do.call("CT_Skull_Strip", L)
  }
  return(res)
}




#' @rdname CT_Skull_Stripper
#' @export
CT_Skull_Stripper_mask <- function(
  ...,
  maskfile = NULL,
  robust = FALSE
){
  if (is.null(maskfile)) {
    maskfile = tempfile(fileext = ".nii.gz")
  }
  # making retimg and reorient for mask
  args = list(...)
  retimg = args$retimg
  if (is.null(retimg)) {
    retimg = TRUE
  }
  reorient = args$reorient
  if (is.null(reorient)) {
    reorient = FALSE
  }

  res = CT_Skull_Stripper(
    ..., robust = robust,
    maskfile = maskfile)
  if (retimg) {
    mask = readnii(maskfile,
                   reorient = reorient)
  }
  L = list(ss_image = res,
           mask = mask)
  return(L)
}
