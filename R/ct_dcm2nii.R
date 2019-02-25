
#' CT DICOM to NIfTI conversion
#'
#' @param basedir (character) directory to get files
#' @param merge_files Should files be merged, passed do \code{\link{dcm2nii}}
#' options
#' @param verbose print diagnostic messages
#' @param ... Additional parameters passed to \code{\link{dcm2nii}}
#'
#' @return A list or singular \code{nifti} image
#' @export
#'
#' @importFrom dcm2niir dcm2nii check_dcm2nii
#' @importFrom neurobase rescale_img check_nifti
ct_dcm2nii = function(basedir = ".", merge_files = TRUE,
                      verbose = TRUE, ...) {
  out = dcm2nii(basedir, merge_files = TRUE, verbose = verbose,
                ...)
  res = check_dcm2nii(out)
  img = lapply(res, function(x){
    if (verbose) {
      message("# reading in image")
    }
    img = check_nifti(res)
    if (verbose) {
      message("# rescaling data")
    }
    img = rescale_img(img, min.val = -1024,
                      max.val = 3071)
  })
  if (length(res) == 1) {
    img = img[[1]]
  }
  return(img)
}