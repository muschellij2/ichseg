#' @title Predict ICH Images
#' @description This function will take the data.frame of predictors and
#' predict the ICH voxels from the model chosen.
#'
#' @param df \code{\link{data.frame}} of predictors.  If \code{multiplier}
#' column does not exist, then \code{\link{ich_candidate_voxels}} will
#' be called
#' @param nim object of class \code{\link{nifti}}, from
#' \code{\link{make_predictors}}
#' @param model model to use for prediction,
#' either the random forest (rf) or logistic
#' @param verbose Print diagnostic output
#' @param native Should native-space predictions be given?
#' @param native_img object of class \code{\link{nifti}}, which
#' is the dimensions of the native image
#' @param transformlist Transforms list for the transformations back to native space.
#' NOTE: these will be inverted.
#' @param interpolator Interpolator for the transformation back to native space
#' @param native_thresh Threshold for re-thresholding binary mask after
#' interpolation
#' @param ... Additional options passsed to \code{\link{ich_preprocess}}
#'
#' @return List of output registered and native space
#' prediction/probability images
#' @importFrom neurobase remake_img
#' @importFrom extrantsr ants_bwlabel
#' @import randomForest
#' @seealso \code{\link{ich_candidate_voxels}}
#' @export
ich_predict = function(df,
                       nim,
                       model = c("rf", "logistic"),
                       verbose = TRUE,
                       native = TRUE,
                       native_img = NULL,
                       transformlist = NULL,
                       interpolator = NULL,
                       native_thresh = 0.5,
                       ...) {

  # if (!have_matlab()) {
  #   stop("MATLAB Path not defined!")
  # }

  cn = colnames(df)
  if (!("multiplier" %in% cn)) {
    df$multiplier = ich_candidate_voxels(df)
  }
  df$Y = NULL
  cc = complete.cases(df)
  if (!all(cc)) {
    warning("NAs or missing in DF, removing")
    for (icn in seq(ncol(df))) {
      x = df[, icn]
      if (!(class(x) %in% c("factor", "character"))) {
        x[ !is.finite(x) ] = 0
      }
      df[, icn] = x
    }
  }
  if (verbose) {
    message("# Making Prediction Images")
  }
  env = as.environment("package:ichseg")

  # Getting modlist for model and cutoff
  modlist.name = paste0(model, "_modlist")
  modlist = env[[modlist.name]]
  mod = modlist$mod
  cutoff = modlist$mod.dice.coef[1, "cutoff"]

  rm(list = c("modlist"))

  # Getting smoothed cutoff
  smoothed_name = paste0("smoothed_", model, "_cutoffs")
  scutoffs = env[[smoothed_name]]
  smoothed_cutoff = scutoffs$mod.dice.coef[1, "cutoff"]
  rm(list = c("scutoffs", "smoothed_name"))

  p = switch(model,
             rf = predict(mod,
                          newdata = df[ df$multiplier, ],
                          type = "prob")[, "1"],
             logistic = predict(mod,
                                df[ df$multiplier, ],
                                type = "response"))

  mult_img = niftiarr(nim, df$multiplier)

  # p = predict(mod, df[ df$multiplier, ], type = "response")
  pimg = remake_img(p,
                    nim,
                    mult_img)

  mask = niftiarr(nim, df$mask)
  pimg = mask_img(pimg, mask)
  sm.pimg  = mean_image(pimg,
                        nvoxels = 1)
  sm.pimg[abs(sm.pimg) <
            .Machine$double.eps ^ 0.5 ] = 0
  sm.pimg = niftiarr(nim, sm.pimg)
  sm.pimg[is.na(sm.pimg)] = 0


  sm.pred = sm.pimg > smoothed_cutoff
  pred = pimg > cutoff

  # cc = spm_bwlabel(pred, k = 100)
  # scc = spm_bwlabel(sm.pred, k = 100)
  cc = ants_bwlabel(img = pred, k = 100, binary = TRUE)
  scc = ants_bwlabel(img = sm.pred, k = 100, binary = TRUE)

  ##############################################################
  # Back to Native Space!
  ##############################################################
  res = list(
    prediction_image = cc,
    smoothed_prediction_image = scc,
    probability_image = pimg,
    smoothed_probability_image = sm.pimg)

  ##############################################################
  # Inverted!
  ##############################################################
  native_res = NULL
  if (native) {
    stopifnot(!is.null(interpolator))
    stopifnot(!is.null(transformlist))
    native_res = lapply(res, function(x){
      ants_apply_transforms(fixed = native_img,
                            moving = x,
                            transformlist = transformlist,
                            interpolator = interpolator,
                            whichtoinvert = c(1)
      )
    })
    native_res$smoothed_prediction_image = neurobase::datatyper(
        native_res$smoothed_prediction_image > native_thresh
    )
    native_res$prediction_image = neurobase::datatyper(
      native_res$prediction_image > native_thresh
    )
  }
  res$cutoff = cutoff
  res$smoothed_cutoff = smoothed_cutoff

  L = list(registered_prediction = res,
           native_prediction = native_res)
  return(L)
}
