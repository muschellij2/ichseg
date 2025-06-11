
#' CT DICOM to NIfTI conversion
#'
#' @param basedir (character) directory to get files
#' @param merge_files Should files be merged, passed do \code{\link{dcm2nii}}
#' options
#' @param verbose print diagnostic messages
#' @param drop_dim passed to \code{\link{readnii}} for dropping empty
#' dimensions
#' @param ... Additional parameters passed to \code{\link{dcm2nii}}
#' @param ignore_roi_if_multiple additional argument to pass to
#' [dcm2niir::check_dcm2nii()] to remove ROI overlays if present
#' @param uncorrected passed to [dcm2niir::check_dcm2nii()] to grab the
#' "uncorrected" scan.  Do not use unless you understand the `dcm2niix`
#' correction process.
#'
#' @return A list or singular \code{nifti} image
#' @export
#'
#' @importFrom dcm2niir dcm2nii check_dcm2nii
#' @importFrom neurobase rescale_img check_nifti
ct_dcm2nii = function(basedir = ".", merge_files = TRUE,
                      verbose = TRUE,
                      drop_dim = TRUE, ...,
                      ignore_roi_if_multiple = FALSE,
                      fail_on_error = FALSE,
                      uncorrected = FALSE) {
  if (!merge_files) {
    warning(
      paste0(
        "ichseg < v0.19.0 had overridden merge_files = FALSE (bug),",
        " please be aware for reproducibility."
      )
    )
  }
  out = dcm2nii(basedir, merge_files = merge_files, verbose = verbose,
                ...)
  if (fail_on_error && out$result > 0) {
    stop("Error in result from dcm2nii and fail_on_error = TRUE")
  }
  res = check_dcm2nii(out, ignore_roi_if_multiple = ignore_roi_if_multiple,
                      uncorrected = uncorrected)
  img = lapply(res, function(x){
    if (verbose) {
      message("# reading in image")
    }
    img = check_nifti(res, drop_dim = drop_dim)
    if (verbose) {
      message("# rescaling data")
    }
    img = rescale_img(img, min.val = -1024,
                      max.val = 3071,
                      drop_dim = drop_dim)
  })
  if (length(res) == 1) {
    img = img[[1]]
  }
  attr(img, "dcm2nii_result") = out$result
  return(img)
}