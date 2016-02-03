#' @title Predict ICH Segmentation
#' @description Will preprocess and predict the ICH voxels
#'
#' @param img CT image, object of class \code{nifti} or
#' character filename
#' @param mask binary brain mask, object of class \code{nifti} or
#' character filename
#' @param model model to use for prediction,
#' either the random forest (rf) or logistic
#' @param verbose Print diagnostic output
#' @param ... Additional options passsed to \code{\link{ich_preprocess}}
#'
#' @return List of output prediction/probability images
#' @import randomForest
#' @export
ich_segment = function(img,
                       mask = NULL,
                       model = c("rf", "logistic"),
                       verbose = TRUE,
                       ...) {

  if (!have_matlab()) {
    stop("MATLAB Path not defined!")
  }
  if (!have.fsl()) {
    stop("FSL Path Not Found!")
  }

  model = match.arg(model)

  if (verbose) {
    message("# Processing The Data")
  }
  # orig.img = img
  preprocess = ich_preprocess(
    img = img,
    mask = mask,
    verbose = verbose, ...)

  timg = preprocess$transformed_image
  tmask = preprocess$transformed_mask > 0.5

  if (verbose) {
    message("# Making Predictors")
  }
  img.pred = make_predictors(timg, mask = tmask,
                             roi = NULL, save_imgs = FALSE,
                             verbose = verbose)
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
    message("# Running ich_predict")
  }
  L = ich_predict(df = df,
                  nim = nim,
                  model = model,
                  native_img = img,
                  native = TRUE,
                  verbose = verbose,
                  transformlist = preprocess$invtransforms,
                  interpolator = preprocess$interpolator)
  L$preprocess = preprocess

  return(L)

}
