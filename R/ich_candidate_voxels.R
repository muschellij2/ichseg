#' @title Create ICH Segmentation Candidate Voxesl
#' @description Takes estimated cutoffs from paper from training data
#' and creates a logical candidate vector
#' @param df \code{data.frame} of image from \code{\link{make_predictors}}
#'
#' @return Logical vector
#' @export
ich_candidate_voxels = function(df){

#   fname = file.path(outdir,
#                     paste0("Reseg_Aggregate_data_cutoffs",
#                            adder, ".Rda"))
#
#   load(file = fname)

#   data(est.cutoffs)
#
#   med.ztemp = median(df$zscore_template[keep.ind])
#
#   df$gr_medztemp = (df$zscore_template >
#                       med.ztemp)

  keepnames = c("zscore3", "zscore2",
                "pct_thresh", "zscore_template")
  for (icut in keepnames) {
    qcuts = ichseg::est.cutoffs[, icut]
    colname = paste0(icut, ".cutoff")
    df[, colname] =
      df[, icut] >= qcuts[1] &
      df[, icut] <= qcuts[2]
  }

  df$include = df$value >= 30 &
    df$value <= 100
  df$zval = df[, "zscore3.cutoff"] &
    df$include &
    df$pct_thresh.cutoff
  df$zval2 = df[, "zscore2.cutoff"] &
    df$zval

  candidate = df$zval2
  return(candidate)
}

