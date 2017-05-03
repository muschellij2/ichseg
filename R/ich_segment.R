#' @title Predict ICH Segmentation
#' @description Will preprocess and predict the ICH voxels
#'
#' @param img CT image, object of class \code{nifti} or
#' character filename
#' @param mask binary brain mask, object of class \code{nifti} or
#' character filename
#' @param model model to use for prediction,
#' either the random forest (rf) or logistic
#' @param save_imgs Logical to save all images that are created as
#' predictors
#' @param outdir Output directory of saved images, needs to be set
#' if \code{save_imgs = TRUE}
#' @param stub Basename to write image names if \code{save_imgs = TRUE}
#' @param verbose Print diagnostic output
#' @param shiny Should shiny progress be called?
#' @param ... Additional options passsed to \code{\link{ich_preprocess}}
#'
#' @return List of output prediction/probability images
#' @export
#' @importFrom fslr have.fsl
ich_segment = function(img,
                       mask = NULL,
                       model = c("rf", "logistic"),
                       save_imgs = FALSE,
                       outdir = NULL,
                       stub = NULL,
                       verbose = TRUE,
                       shiny = FALSE,
                       ...) {

  # if (!have_matlab()) {
  #   stop("MATLAB Path not defined!")
  # }

  if (!have.fsl()) {
    stop("FSL Path Not Found!")
  }

  model = match.arg(model)

  if (verbose) {
    msg = "# Processing The Data"
    message(msg)
  }
  if (shiny) {
    shiny::setProgress(message = msg, value = 0)
  }
  # orig.img = img
  preprocess = ich_preprocess(
    img = img,
    mask = mask,
    verbose = verbose,
    shiny = shiny,
    ...)

  timg = preprocess$transformed_image
  tmask = preprocess$transformed_mask > 0.5

  if (verbose) {
    msg = "# Making Predictors"
    message(msg)
  }
  if (shiny) {
    shiny::setProgress(message = msg, value = 2/3)
  }
  if (save_imgs){
    if (is.character(img)){
      if (is.null(stub)){
        stub = paste0(nii.stub(img, bn = TRUE), "_reg_")
      }
    }
  }
  img.pred = make_predictors(timg, mask = tmask,
                             roi = NULL, save_imgs = save_imgs,
                             stub = stub,
                             outdir = outdir,
                             verbose = verbose,
                             shiny = shiny)
  df = img.pred$df
  nim = img.pred$nim
  rm(list = "img.pred")
  gc()

  # data(MOD)
  ##############################################################
  # Making prediction images
  ##############################################################
  # grabbing the evnironment to extract exported stuff
  if (verbose) {
    msg = "# Running ich_predict"
    message(msg)
  }
  if (shiny) {
    shiny::setProgress(message = msg, value = 3/3 - 0.3)
  }
  L = ich_predict(df = df,
                  nim = nim,
                  model = model,
                  native_img = img,
                  native = TRUE,
                  verbose = verbose,
                  transformlist = preprocess$invtransforms,
                  interpolator = preprocess$interpolator,
                  shiny = shiny)
  L$preprocess = preprocess
  if (shiny) {
    shiny::setProgress(value = 3/3)
  }
  return(L)

}
