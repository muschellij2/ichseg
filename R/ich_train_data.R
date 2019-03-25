#' Create ICH Train Data
#'
#' @param imgs CT image, object of class \code{nifti} (list) or
#' character filename
#' @param verbose Print diagnostic output
#' @param cutoffs cutoffs to pass to \code{\link{ich_candidate_voxels}}
#' @param run_candidate Run \code{\link{ich_candidate_voxels}}
#' to get candidate voxels
#' @param rois ROI image, object of class \code{nifti} (list) or
#' character filename
#'
#' @param ... Additional options passsed to \code{\link{ich_preprocess}}
#'
#' @note This is a simple wrapper for \code{\link{ich_process_predictors}}
#' @return List of one \code{data.frame} with the predictors
#' @export
ich_train_data = function(
  imgs,
  run_candidate = TRUE,
  cutoffs = ichseg::est.cutoffs,
  rois = NULL,
  ...,
  verbose = TRUE) {

  imgs = checkimg(imgs)
  if (!is.null(rois)) {
    rois = checkimg(rois)
  } else {
    rois = lapply(seq_along(imgs), function(x) NULL)
  }
  stopifnot(length(imgs) == length(rois))

  results = mapply(function(img, roi) {
    L = ich_process_predictors(
      img = img,
      roi = roi,
      ...,
      verbose = verbose)
    df = L$img.pred$df
    df$filename = img
    if (run_candidate) {
      df$candidate = ich_candidate_voxels(df, cutoffs = cutoffs)
    }
    df
  }, imgs, rois, SIMPLIFY = FALSE)

  results = do.call("rbind", results)
  return(results)
}


#' ICH Train Model
#'
#' @param df A \code{data.frame} from \code{\link{ich_train_data}} or
#' the \code{$img.pred$df} output from  \code{\link{ich_process_predictors}}
#' output list.
#' @param do.trace passed to \code{\link{randomForest}} for training verbosity
#' @param ... additional arguments passed to \code{\link{randomForest}}
#'
#' @return A \code{randomForest} object
#' @importFrom stats as.formula
#' @export
ich_train_model = function(
  df,
  do.trace = TRUE, ...) {


  cn = c("moment1", "moment2", "skew", "kurtosis", "value", "thresh",
         "zscore1", "zscore2", "zscore3", "win_z", "pct_thresh", "prob_img",
         "pct_zero_neighbor", "any_zero_neighbor", "dist_centroid", "smooth5",
         "smooth10", "smooth20", "zscore_template", "flipped_value", "Y")
  keep_cn = intersect(cn, colnames(df))
  df = df[, keep_cn]
  df$Y = factor(df$Y)

  sds = apply(as.matrix(df), 2, sd)
  dropper = sds == 0
  novar = ""
  if (any(dropper)) {
    novar = keep_cn[dropper]
    novar = c("", novar)
    novar = paste0(novar, collapse = " - ")
  }
  formstr = paste0("Y ~ .", novar)
  formula = as.formula(formstr)


  rf.mod = randomForest(
    formula = formula,
    data = df,
    ...,
    do.trace = do.trace)

  return(rf.mod)
}

