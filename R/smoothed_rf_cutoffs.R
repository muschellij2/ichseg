#' @title Estimated Random Forest List of Cutoffs for Smoothed Data for ICH
#'
#' @description A list containing the random forest model cutoffs
#' for predicting ICH after smoothing and other factors for
#' prediction.
#'
#' @format A list with 7 elements, which are:
#' \describe{
#' \item{mod.pauc}{the partial AUC estimated with the training data }
#' \item{mod.pauc.cutoff}{the sens/spec for partial AUC cutoff }
#' \item{mod.sens.cutoff}{the sensitivity and cutoff estimated with the training datan }
#' \item{mod.dice.cutoff}{the Dice Similarity Index and cutoff estimated with the training data }
#' \item{mod.acc}{the accuracy and cutoff estimated with the training data}
#' }
"smoothed_rf_cutoffs"

#' @rdname smoothed_rf_cutoffs
"smoothed_big_rf_cutoffs"
