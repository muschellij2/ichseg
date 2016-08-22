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
#' @param ... Additional options passsed to \code{\link{CT_Skull_Strip_robust}}
#' or \code{\link{CT_Skull_Strip}}
#'
#' @return List of output images and transformations
#' @importFrom fslr check_nifti mask_img window_img
#' @importFrom extrantsr registration ants_apply_transforms
#' @export
ich_preprocess = function(img,
                          skull_strip = TRUE,
                          robust = TRUE,
                          mask = NULL,

                          template.file = system.file(
                            "scct_unsmooth_SS_0.01.nii.gz",
                            package = "ichseg"),
                          outprefix = NULL,
                          typeofTransform = c("Rigid", "Affine"),
                          interpolator = "Linear",
                          outfile = NULL,
                          verbose = TRUE,
                          ...) {

  if (skull_strip & is.null(mask)) {
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

  if (is.null(outprefix)){
    outprefix = tempfile()
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
    verbose = verbose)

  omask = ants_apply_transforms(fixed = template.file,
                                moving = mask,
                                typeofTransform = typeofTransform,
                                interpolator = interpolator,
                                transformlist = res$fwdtransforms)
  res$mask = mask
  res$ss_image = ss
  res$transformed_image = res$outfile
  res$transformed_mask = omask
  res$outfile = NULL

  return(res)

}
