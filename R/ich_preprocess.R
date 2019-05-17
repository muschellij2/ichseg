#' @title Preprocess data for ICH Segmentation
#' @description Will do skull stripping and registration to the CT template
#'
#' @param img CT image, object of class \code{nifti} or
#' character filename
#' @param skull_strip Should the image be skull stripped? If not,
#' a mask must be specified.  If mask specified, will override this
#' option
#' @param mask binary brain mask, object of class \code{nifti} or
#' character filename
#' @param robust If skull stripping, should
#' \code{\link{CT_Skull_Strip_robust}}
#' be used vs. \code{\link{CT_Skull_Strip}}
#' @param template.file Template to register to (Skull stripped template)
#' @param outprefix Passed to \code{\link[extrantsr]{registration}} for
#' the output transformations
#' @param typeofTransform Transformation to the template
#' @param interpolator Interpolation to be done after transformation
#' @param outfile Output filename for transformed registered image
#' @param verbose Print diagnostic output
#' @param shiny Should shiny progress be called?
#' @param ... Additional options passsed to \code{\link{CT_Skull_Strip_robust}}
#' or \code{\link{CT_Skull_Strip}}
#' @param roi Filename of ROI, which will be transformed
#' @param n4_correct Should N4 bias-field correction be done after
#' skull-stripping
#'
#'
#' @return List of output images and transformations
#' @importFrom neurobase check_nifti mask_img window_img
#' @importFrom extrantsr registration ants_apply_transforms bias_correct
#' @export
ich_preprocess = function(
  img,
  skull_strip = TRUE,
  robust = TRUE,
  mask = NULL,
  n4_correct = FALSE,
  template.file = system.file(
    "scct_unsmooth_SS_0.01.nii.gz",
    package = "ichseg"),
  outprefix = NULL,
  typeofTransform = c("Rigid", "Affine"),
  interpolator = "Linear",
  outfile = NULL,
  verbose = TRUE,
  shiny = FALSE,
  roi = NULL,
  ...) {

  if (skull_strip & is.null(mask)) {
    if (shiny) {
      shiny::incProgress(message = "Skull Stripping Image")
    }
    if (robust) {
      ss = CT_Skull_Strip_robust(img, retimg = TRUE, ...)
    } else {
      ss = CT_Skull_Strip(img, retimg = TRUE, ...)
    }
    mask = ss > 0
  } else {
    stopifnot(!is.null(mask))
  }
  mask = check_nifti(mask)
  mask = mask > 0

  img = check_nifti(img)
  ss = mask_img(img, mask)
  ss = window_img(ss, window = c(0, 100), replace = "zero")


  if (n4_correct) {
    if (verbose) {
      msg = "# Running N4 Correction"
      message(msg)
    }
    n4_img = extrantsr::bias_correct(
      ss, correction = "N4",
      mask = mask,
      verbose = verbose > 1)
    ss = n4_img
  }


  if (is.null(outprefix)) {
    outprefix = tempfile()
  }
  if (shiny) {
    shiny::incProgress(message = "Rigid-Body Registration")
  }
  typeofTransform = match.arg(typeofTransform)
  res = registration(
    filename = ss,
    skull_strip = FALSE,
    correct = FALSE,
    outfile = outfile,
    retimg = TRUE,
    typeofTransform = typeofTransform,
    template.file = template.file,
    interpolator = interpolator,
    remove.warp = FALSE,
    outprefix = outprefix,
    verbose = verbose > 1)

  omask = ants_apply_transforms(fixed = template.file,
                                moving = mask,
                                # typeofTransform = typeofTransform,
                                interpolator = interpolator,
                                transformlist = res$fwdtransforms)
  if (!is.null(roi)) {
    roi_interpolator = "genericLabel"
    oroi = ants_apply_transforms(fixed = template.file,
                                 moving = roi,
                                 # typeofTransform = typeofTransform,
                                 interpolator = roi_interpolator,
                                 transformlist = res$fwdtransforms)
    res$roi_interpolator = roi_interpolator
  } else {
    oroi = NULL
  }

  res$mask = mask
  res$transformed_roi = oroi
  res$ss_image = ss
  res$transformed_image = res$outfile
  res$transformed_mask = omask
  res$outfile = NULL

  return(res)

}

#' @export
#' @rdname ich_preprocess
ich_cnn_preprocess = function(
  ...,
  template.file = system.file(
    "scct_unsmooth_SS_0.01_128x128x128.nii.gz",
    package = "ichseg")
) {
  args = list(...)
  nargs = names(args)
  if ("skull_strip" %in% nargs) {
    skull_strip = args$robust
    if (!skull_strip) {
      warning("ich_cnn_preprocess needs skull_strip = TRUE")
    }
    args$skull_strip = TRUE
  }
  if ("robust" %in% nargs) {
    robust = args$robust
    if (!robust) {
      warning("ich_cnn_preprocess needs robust = TRUE")
    }
    args$robust = TRUE
  }
  args$template.file = template.file
  args$smooth_before_threshold = TRUE
  args$smooth.factor = 1
  args$remove.neck = TRUE
  args$recog = FALSE
  args$nvoxels = 0
  res = do.call(ich_preprocess, args = args)
  return(res)
}
