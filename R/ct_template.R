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
#' \code{mask} is brain mask. Values with \code{cormack} or \code{stripped}
#' are in 2 millimeter resolution, and smoothed.
#' @return Character filename
#' @export
#' @examples
#' types = c("brain", "mask", "image",
#' "cormack",
#' "stripped_cormack",
#' "stripped_hu",
#' "stripped_hu8")
#' sapply(types, ct_template_fname)
ct_template_fname = function(type = c("brain", "mask", "image",
                                      "cormack",
                                      "stripped_cormack",
                                      "stripped_hu",
                                      "stripped_hu8")){
  type = match.arg(type)
  if (type %in% c("cormack",
                  "stripped_cormack",
                  "stripped_hu",
                  "stripped_hu8")) {
    stub = paste0("scct_", type)
  } else {
    stub = "scct_unsmooth"
    brain = type %in% "brain"
    mask = type %in% "mask"
    if (brain | mask) {
      stub = paste0(stub, "_SS_0.01")
    }
    if (mask) {
      stub = paste0(stub, "_Mask")
    }
  }
  stub = paste0(stub, ".nii.gz")
  res = system.file(stub,
                    package = "ichseg")
  if (!file.exists(res)) {
    stop("Template File not found!")
  }
  return(res)
}