
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
#' @param fail_on_error if the exit code from [dcm2niir::dcm2nii] is not
#' `0` and `fail_on_error = TRUE`, function will error out
#' @param cleanup delete any temporary directories created during conversion,
#' overrides `dcm2nii` cleanup option.
#'
#' @return A list or singular \code{nifti} image
#' @export
#'
#' @importFrom dcm2niir dcm2nii check_dcm2nii
#' @importFrom neurobase rescale_img check_nifti
ct_dcm2nii = function(basedir = ".", merge_files = TRUE,
                      verbose = TRUE,
                      drop_dim = TRUE, ...,
                      cleanup = FALSE,
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
  args = list(...)
  if (is.null(args$copy_files)) {
    # use the default of command
    args$copy_files = formals(dcm2niir::dcm2nii)$copy_files
  }
  if (cleanup && !args$copy_files) {
    stop(
      paste0(
        "cleanup = TRUE requires copy_files = TRUE in dcm2nii,",
        " otherwise original data may be deleted.."
      )
    )
  }
  out = dcm2nii(basedir,
                merge_files = merge_files, verbose = verbose,
                cleanup = FALSE,
                ...)
  img = ct_process_dcm2nii(
    out,
    verbose = verbose,
    drop_dim = drop_dim,
    cleanup = cleanup,
    ignore_roi_if_multiple = ignore_roi_if_multiple,
    fail_on_error = fail_on_error,
    uncorrected = uncorrected
  )
  return(img)
}

#' @export
#' @rdname ct_dcm2nii
#' @param out output from [ct_dcm2nii]
ct_process_dcm2nii = function(out,
                              verbose = TRUE,
                              drop_dim = TRUE,
                              cleanup = FALSE,
                              ignore_roi_if_multiple = FALSE,
                              fail_on_error = FALSE,
                              uncorrected = FALSE) {
  if (fail_on_error && out$result > 0) {
    stop("Error in result from dcm2nii and fail_on_error = TRUE")
  }
  dir_temp = attr(out$result, "copy_files_temporary_directory")
  res = check_dcm2nii(out, ignore_roi_if_multiple = ignore_roi_if_multiple,
                      uncorrected = uncorrected)
  names(res) = res
  img = lapply(res, function(x){
    if (verbose) {
      message("# reading in image")
    }
    img = check_nifti(res, drop_dim = drop_dim, fast = TRUE)
    if (verbose) {
      message("# rescaling data")
    }
    img = rescale_img(img, min.val = -1024,
                      max.val = 3071,
                      drop_dim = drop_dim)
  })
  if (length(res) == 1) {
    img = img[[1]]
    names(img) = res
  }
  if (cleanup) {
    if (verbose) {
      message("# cleaning up temporary files")
    }
    if (!is.null(dir_temp) && dir.exists(dir_temp)) {
      unlink(dir_temp, recursive = TRUE, force = TRUE)
    }
  }
  attr(img, "dcm2nii_result") = out$result
  attr(img, "files") = res
  return(img)
}