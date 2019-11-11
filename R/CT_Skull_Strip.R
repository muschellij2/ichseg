#' @title CT Skull Stripping within R
#' @description Skull Stripping (using FSL's BET) a CT file using \code{fslr}
#' functions
#' @param img (character) File to be skull stripped or object of class
#' nifti
#' @param outfile (character) output filename
#' @param keepmask (logical) Should we keep the mask?
#' @param maskfile (character) Filename for mask (if \code{keepmask = TRUE}).
#' If \code{NULL}, then will do \code{paste0(outfile, "_Mask")}.
#' @param inskull_mesh (logical) Create inskull_mesh file from bet?
#' (Warning - will take longer)
#' This an exterior surface of the brain. (experimental)
#' Also, if \code{outfile} is \code{NULL}, then this will be created in
#' a temporary directory and not be retrieved.
#' @param retimg (logical) return image of class nifti
#' @param reorient (logical) If retimg, should file be reoriented when read in?
#' Passed to \code{\link{readNIfTI}}.
#' @param intern (logical) pass to \code{\link{system}}
#' @param betcmd (character) bet command to be used, see \code{\link{fslbet}}
#' @param opts (character) additional options to \code{\link{fslbet}}
#' @param presmooth (logical) indicator if pre-smoothing should be
#' done before BET
#' @param remask (logical) Mask the smoothed image with HU mask from initial
#' step?
#' @param refill (logical) indicator to post-smooth mask and then fill
#' @param refill.thresh (numeric) Value to threshold post-smoothed mask
#' @param sigma (integer) size of Gaussian kernel passed to
#' \code{\link{fslsmooth}} if \code{presmooth} is \code{TRUE}
#' @param lthresh (default: 0) Lower value to threshold CT
#' \code{\link{fslthresh}}
#' @param uthresh (default: 100) Upper value to threshold CT
#' \code{\link{fslthresh}}
#' @param verbose (logical) Should diagnostic output be printed?
#' @param ... additional arguments passed to \code{\link{fslbet}}.
#' @param smooth_before_threshold Should the image be smoothed before
#' thresholding?  This can be useful for bone-window scans.
#'
#' @return character or logical depending on intern
#' @importFrom fslr fslthresh fslmaths fslfill fslsmooth fslmask fslbet
#' @importFrom neurobase nii.stub
#' @importFrom fslr get.imgext
#' @export
CT_Skull_Strip <- function(
  img,
  outfile = NULL,
  keepmask = TRUE,
  maskfile = NULL,
  inskull_mesh = FALSE,
  retimg = TRUE,
  reorient = FALSE,
  intern=TRUE,
  betcmd = "bet2",
  opts = "-f 0.01 -v",
  presmooth = TRUE,
  remask = TRUE,
  refill = FALSE,
  refill.thresh = 0.75,
  sigma = 1,
  lthresh = 0,
  uthresh = 100,
  smooth_before_threshold = FALSE,
  verbose=TRUE,
  ...
){
  if (retimg) {
    if (is.null(outfile)) {
      outfile = tempfile()
    }
  } else {
    stopifnot(!is.null(outfile))
  }
  orig_img = img


  if (smooth_before_threshold) {
    if (verbose) {
      message(paste0("# Smoothing before Thresholding"))
    }
    pm_file = tempfile()
    result = fslsmooth(img,
                       outfile = pm_file,
                       sigma = sigma,
                       retimg = FALSE,
                       verbose = verbose)
    img = pm_file
  }

  if (verbose) {
    message(paste0("# Thresholding Image to ",
                   lthresh, "-", uthresh, "\n"))
  }

  tfile = tempfile()
  tcopy = tempfile(fileext = ".nii.gz")
  outfile = nii.stub(outfile)
  run = fslthresh(img, thresh = lthresh, uthresh = uthresh,
                  outfile = tfile,
                  retimg = FALSE,
                  intern = FALSE,
                  verbose = verbose)
  if (!any(c(readnii(tfile)) > 0)) {
    stop("No positive values in the thresholded output!")
  }
  if (verbose) {
    message(paste0("# Thresholding return: ", run, "\n"))
  }

  if (remask) {
    if (verbose) {
      message(paste0("# Absolute value so fslmaths -bin keeps all mask",
                     " even if lthresh < 0\n"))
    }
    absfile = tempfile()
    run = fslmaths(tfile, outfile = absfile,
                   retimg = FALSE,
                   intern = FALSE,
                   opts = "-abs",
                   verbose = verbose)


    if (verbose) {
      message(paste0("# Abs return: ", run, "\n"))
    }


    if (verbose) {
      message(paste0(
        "# Creating binary mask to remask after filling\n"))
    }
    bonefile = tempfile()
    #   fslbin(outfile, retimg = FALSE, outfile = bonefile, intern=FALSE)
    fslfill(file = absfile,
            bin = TRUE,
            outfile = bonefile,
            retimg = FALSE,
            intern = FALSE,
            verbose = verbose)
  }

  ### Must prefill for the presmooth - not REALLY necessary, but if you're
  ### smoothing, you likely have noisy data.

  #   if (isTRUE(presmooth)){
  #     if (!isTRUE(prefill)){
  #       warning("if presmooth, must have prefill = TRUE")
  #       prefill= TRUE
  #     }
  #   }
  #   #### prefill will mask the img
  #   if (isTRUE(prefill)){
  #     if (verbose){
  #       message(paste0("Remasking 0 - 100 after filling holes\n"))
  #     }
  #     run = fslmask(img,
  #                   mask = bonemask,
  #                   outfile = outfile,
  #                   retimg = FALSE,
  #                   intern = FALSE)
  #   }


  if (isTRUE(presmooth)) {

    if (verbose) {
      message(paste0("# Presmoothing image\n"))
    }
    run = fslsmooth(tfile,
                    outfile = tfile,
                    sigma = sigma,
                    retimg = FALSE,
                    intern = FALSE,
                    verbose = verbose)
    if (verbose) {
      message(paste0("# Pre Smoothing Diagnostic: ", run, "\n"))
    }

    if (remask) {
      if (verbose) {
        message(paste0("# Remasking Smoothed Image\n"))
      }
      run = fslmask(tfile,
                    mask = bonefile,
                    outfile = tfile,
                    retimg = FALSE,
                    intern = FALSE,
                    verbose = verbose)
      if (verbose) {
        message(paste0("# Pre Smoothing Diagnostic: ", run, "\n"))
      }

    }
  }
  if (verbose) {
    message(paste0("copied file: ", tcopy))
    file.copy(tfile, tcopy)
    file.copy(paste0(tfile, ".nii"), tcopy, overwrite = FALSE)
    file.copy(paste0(tfile, ".nii.gz"), tcopy, overwrite = FALSE)
  }

  #### Different options for bet
  if (inskull_mesh) {
    opts = paste0(opts, " -A")
    if (betcmd != "bet") {
      warning("Can't use bet2 with inskull mesh, using bet")
    }
    betcmd = "bet"
  }
  betintern = TRUE
  if (verbose) {
    message(paste0("# Running ", betcmd, "\n"))
    betintern = FALSE
  }
  runbet = fslbet(infile = tfile,
                  outfile = tfile,
                  retimg = FALSE, intern = betintern,
                  betcmd = betcmd,
                  opts = opts,
                  verbose = verbose, ...)

  if (verbose) {
    message(paste0("# ", betcmd, " output: ", runbet, "\n"))
  }
  ###############
  # Cleanup
  ##############
  ext = get.imgext()

  if (isTRUE(inskull_mesh)) {
    if (verbose) {
      message("# Deleting extraneous files\n")
    }
    end.names = c("inskull_mask", "outskin_mask", "outskin_mesh",
                  "outskull_mask", "outskull_mesh",
                  "skull_mask")
    end.names = paste0(end.names, ext)
    end.names = c(end.names,
                  "inskull_mesh.off",
                  "outskin_mesh.off",
                  "outskull_mesh.off",
                  "mesh.vtk")
    end.names = paste(tfile, end.names, sep = "_")
    file.remove(end.names)
  }

  #####################
  # Filling in holes of the mask (like choroid plexis and CSF)
  #####################
  if (verbose) {
    message("# Using fslfill to fill in any holes in mask \n")
  }
  if (is.null(maskfile)) {
    maskfile = nii.stub(outfile)
    maskfile = paste0(maskfile, "_Mask")
  }
  stopifnot(inherits(maskfile, "character"))
  # outmask = paste(outfile, "inskull_mask", sep="_")
  # outmask = paste0(outmask, ext)
  # file.rename(outmask, paste0(maskfile, ext))
  absfile = tempfile()
  run = fslmaths(tfile, outfile = absfile,
                 retimg = FALSE,
                 intern = FALSE,
                 verbose = verbose, opts = "-abs")

  runmask = fslfill(file = absfile, bin = TRUE,
                    outfile = maskfile,
                    retimg = FALSE,
                    intern = FALSE,
                    verbose = verbose)
  if (refill) {
    smfile = tempfile()
    fslmaths(maskfile, retimg = FALSE, outfile = smfile,
             opts = "-kernel boxv 7 -fmean",
             verbose = verbose)
    ### smooth and if you add to 0, then > .5 (if 1 before > 1 now, then bin)
    addopt = sprintf(" -add %s -thr %f -bin", smfile, refill.thresh)
    fslmaths(maskfile, retimg = FALSE, outfile = maskfile,
             opts = addopt,
             verbose = verbose)
  }

  if (verbose) {
    message(paste0("# fslfill output: ", runmask, "\n"))
  }

  if (verbose) {
    message("# Using the filled mask to mask original image\n")
  }
  res = fslmask(
    file = orig_img,
    mask = maskfile,
    outfile = outfile,
    retimg = retimg,
    intern = intern,
    reorient = reorient,
    verbose = verbose
  )

  if (!keepmask) {
    if (verbose) {
      message("# Removing Mask File\n")
    }
    maskfile = nii.stub(maskfile)
    maskfile = paste0(maskfile, ext)
    file.remove(maskfile)
  }
  return(res)
}
