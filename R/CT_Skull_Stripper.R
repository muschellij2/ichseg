#' @title CT Skull Stripping Wrapper
#' @description Simple Wrapper for \code{\link{CT_Skull_Strip}} and
#' \code{\link{CT_Skull_Strip_robust}}, with simple flag for robustness
#' @param ... (character) File to be skull stripped or object of class
#' nifti
#' @param robust (logical) output filename
#' @export
CT_Skull_Stripper <- function(
  ...,
  robust = FALSE
){
  
  if (robust){
    CT_Skull_Strip_robust(...)
  } else {
    CT_Skull_Strip(...)
  }
}


