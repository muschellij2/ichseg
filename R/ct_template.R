#' @title Read in CT Template
#' @description Wrapper for returning a \code{nifti} object of the CT Template
#' @param ... Arguments passed to \code{\link{ct_template_fname}}
#'
#' @return Object of class \code{\link{nifti}}
#' @export
ct_template = function(...){
  fname = ct_template_fname(...)
  res = readnii(fname)
  return(res)
}

#' @title Get CT Template Filename
#' @description Wrapper for returning the filename of the CT Template
#' @param type Type of image.  \code{image} is the original image,
#' \code{brain} is the skull stripped image,
#' \code{mask} is brain mask.
#' @return Character filename
#' @export
ct_template_fname = function(type = c("brain", "mask", "image")){
  stub = "scct_unsmooth"
  type = match.arg(type)
  brain = type %in% "brain"
  mask = type %in% "mask"
  if (brain | mask) {
    stub = paste0(stub, "_SS_0.01")
  }
  if (mask) {
    stub = paste0(stub, "_Mask")
  }
  stub = paste0(stub, ".nii.gz")
  res = system.file(stub,
    package = "ichseg")
  if (!file.exists(res)) {
    stop("Template File not found!")
  }
  return(res)
}