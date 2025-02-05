#' @title Robust CT Skull Stripping
#' @description Skull Stripping (using FSL's BET) a CT file
#' using \code{fslr}
#' functions and robustified by registration
#' @param img (character) File to be skull stripped or object of class
#' nifti
#' @param outfile (character) output filename
#' @param keepmask (logical) Should we keep the mask?
#' @param maskfile (character) Filename for mask
#' (if \code{keepmask = TRUE}).
#' If \code{NULL}, then will do \code{paste0(outfile, "_Mask")}.
#' @param retimg (logical) return image of class nifti
#' @param reorient (logical) If retimg, should file be
#' reoriented when read in?
#' Passed to \code{\link{readNIfTI}}.
#' @param int Fractional Intensity passed to
#' \code{\link{CT_Skull_Strip}} and
#' subsequently \code{\link{fslbet}}.
#' @param lthresh (default: 0) Lower value to threshold CT
#' \code{\link{fslthresh}}
#' @param uthresh (default: 100) Upper value to threshold CT
#' \code{\link{fslthresh}}
#' @param nvoxels Number of voxels to dilate/erode.
#' See \code{\link{dil_ero}}.
#' If \code{nvoxels = 0}, then no smoothing is done.
#' @param remove.neck Run \code{\link{remove_neck}} to register
#' the template to a
#' thresholded image to remove neck slices.
#' @param remover if \code{remove.neck = TRUE}, then which function
#' would you like to use to remove the neck
#' @param smooth.factor Smoothing factor for \code{\link{fslbet}}.
#' See \code{-w}
#' option in \code{fslbet.help()}.
#' @param verbose (logical) Should diagnostic output be printed?
#' @param opts Not used
#' @param template.file Template to warp to original image
#' space, passed to
#' \code{\link{remove_neck}}
#' @param template.mask Mask of template to use as rough
#' brain mask, passed
#' to \code{\link{remove_neck}}
#' @param ... additional arguments passed to
#' \code{\link{CT_Skull_Strip}} or
#' \code{\link{remove_neck}}.
#' @param recog Re-estimate the center of gravity (COG) and
#' skull strip.
#' @param mask_to_background When masking, should the values outside the
#' mask be set to 0 (default) or -1024 (when TRUE)
#' @return Skull-stripped \code{nifti} object
#' @note This function first thresholds an image, runs a rigid
#' registration
#' (default in \code{\link{remove_neck}}) to drop any slices
#' below the transformed
#' skull stripped template to remove neck slices.  The neck-removed
#' image is
#' then skull stripped using defaults in \code{\link{CT_Skull_Strip}}.
#'  A new
#' center of gravity is estiamted using \code{\link{cog}}, then the
#' image is
#' skull stripped again using the new cog and the smoothness factor
#' (passed to \code{-w} argument in BET). After the skull stripped
#' mask is
#' created, the image is dilated and eroded using
#' \code{\link{dil_ero}} to
#' fill holes using a box kernel with the number of voxels
#' \code{nvoxels} in
#' all 3 directions.
#' @importFrom extrantsr remove_neck double_remove_neck
#' @importFrom neurobase check_outfile readnii writenii cog
#' @importFrom oro.nifti readNIfTI drop_img_dim
#' @export
CT_Skull_Strip_robust <- function(
  img,
  outfile = NULL,
  keepmask = TRUE,
  maskfile = NULL,
  retimg = TRUE,
  reorient = FALSE,
  int = "0.01",
  lthresh = 0,
  uthresh = 100,
  nvoxels = 5,
  remove.neck = TRUE,
  remover = c("remove_neck", "double_remove_neck"),
  smooth.factor = 2,
  recog = TRUE,
  verbose = TRUE,
  opts = NULL,
  mask_to_background = FALSE,
  template.file =
    system.file("scct_unsmooth_SS_0.01.nii.gz",
                package = "ichseg"),
  template.mask =
    system.file("scct_unsmooth_SS_0.01_Mask.nii.gz",
                package = "ichseg"),
  ...
){

  ###########################################
  # Checking outfile input/output
  ###########################################
  outfile = check_outfile(outfile = outfile,
                          retimg = retimg,
                          fileext = "")
  outfile = path.expand(outfile)
  tfile = tempfile()

  ### Working on maskfile
  if (is.null(maskfile)){
    maskfile = nii.stub(outfile)
    maskfile = paste0(maskfile, "_Mask")
  }
  maskfile = nii.stub(maskfile)
  stopifnot(inherits(maskfile, "character"))

  #############################
  # Threshold Image
  #############################
  img = check_nifti(img, reorient = reorient)
  thresh_img = niftiarr(img, img * (img > lthresh & img < uthresh))
  if (!any(c(thresh_img) > 0)) {
    stop("No positive values in the thresholded output!")
  }

  thresh_img = drop_img_dim(thresh_img)
  thresh = checkimg(thresh_img)
  rm(thresh_img)
  ## need to do rep.value = 0 because template is like that.

  #############################
  # Removing Neck
  #############################
  if (remove.neck) {
    if (verbose) {
      message(paste0("# Removing Neck\n"))
    }
    L = list(
      file = thresh,
      template.file = template.file,
      template.mask = template.mask,
      rep.value = 0,
      verbose = verbose,
      ret_mask = TRUE,
      swapdim = TRUE,
      ...)
    remover = match.arg(remover)
    neck_mask = do.call(remover, args = L)
    # neck_mask = remove_neck(thresh,
    #                         rep.value = 0,
    #                         template.file = template.file,
    #                         template.mask = template.mask,
    #                         ret_mask = TRUE,
    #                         swapdim = TRUE,
    #                         verbose = verbose,
    #                         ...)
  } else {
    neck_mask = niftiarr(img, array(1, dim = dim(img)))
  }
  add_value = 0
  if (mask_to_background) {
    add_value = 1024
  }
  noneck = mask_img(img + add_value, neck_mask) - add_value
  noneck = drop_img_dim(noneck)
  noneck = checkimg(noneck)

  rm(neck_mask)
  #############################
  # Skull Stripping no-neck image
  #############################
  if (verbose) {
    message(paste0("# Skull Stripping for COG\n"))
  }
  ss = CT_Skull_Strip(noneck, outfile = tfile, retimg = TRUE,
                      maskfile = maskfile,
                      lthresh = lthresh, uthresh = uthresh,
                      opts = paste("-f ", int,
                                   ifelse(verbose, " -v", "")),
                      verbose = verbose,
                      keepmask = TRUE, reorient = reorient, ...)

  ssmask = readNIfTI(maskfile,
                     reorient = reorient)
  if (recog) {
    #############################
    # Setting new center of gravity and rerunning with smoothness factor
    #############################
    cog = cog(ssmask,
              ceil = TRUE)
    if (verbose){
      message(paste0("# Skull Stripping with new cog\n"))
    }
    tmp_out = CT_Skull_Strip(noneck, outfile = tfile, retimg = TRUE,
                        opts = paste("-f ", int,
                                     ifelse(verbose, "-v", ""),
                                     "-w ", smooth.factor,
                                     paste(c("-c", cog),
                                           collapse=" ")),
                        maskfile = maskfile,
                        keepmask = TRUE,
                        reorient = reorient,
                        verbose = verbose,
                        ...)
    rm(tmp_out)
    ssmask = readNIfTI(maskfile,
                       reorient = reorient)
  }

  #############################
  # Filling the mask
  #############################
  if (nvoxels > 0){
    if (verbose) {
      message(paste0("# Filling Holes \n"))
    }
    ssmask = dil_ero(ssmask,
                     retimg = TRUE,
                     nvoxels = nvoxels,
                     verbose = verbose)
  }
  ss = mask_img(img + add_value, ssmask) - add_value

  ss = drop_img_dim(ss)

  writenii(ss,
           filename = outfile)

  #############################
  # Removing mask if keepmask = FALSE
  #############################
  if (!keepmask) {
    if (verbose) {
      message("# Removing Mask File\n")
    }
    maskfile = nii.stub(maskfile)
    ext = get.imgext()
    maskfile = paste0(maskfile, ext)
    file.remove(maskfile)
  }
  if (!retimg) {
    outfile = nii.stub(outfile)
    ext = get.imgext()
    outfile = paste0(outfile, ext)
    ss = outfile
  }
  return(ss)
}




#' @rdname CT_Skull_Strip_robust
#' @note \code{CT_Skull_Strip_register} removes the neck, registers the
#' image to the template, using a rigid-body transformation,
#' runs the skull stripper to get a mask, then transforms the mask
#' back to native space.
#' @export
CT_Skull_Strip_register <- function(
  img,
  outfile = NULL,
  keepmask = TRUE,
  maskfile = NULL,
  retimg = TRUE,
  reorient = FALSE,
  lthresh = 0,
  uthresh = 100,
  remove.neck = TRUE,
  remover = c("remove_neck", "double_remove_neck"),
  verbose = TRUE,
  mask_to_background = FALSE,
  ...
){

  ###########################################
  # Checking outfile input/output
  ###########################################
  outfile = check_outfile(outfile = outfile,
                          retimg = retimg,
                          fileext = "")
  outfile = path.expand(outfile)
  tfile = tempfile()

  ### Working on maskfile
  if (is.null(maskfile)) {
    maskfile = nii.stub(outfile)
    maskfile = paste0(maskfile, "_Mask")
  }
  maskfile = nii.stub(maskfile)
  stopifnot(inherits(maskfile, "character"))

  #############################
  # Threshold Image
  #############################
  if (verbose) {
    message(paste0("# Thresholding Image to ",
                   lthresh, "-", uthresh, "\n"))
  }
  img = check_nifti(img, reorient = reorient)
  thresh_img = niftiarr(img, img * (img > lthresh & img < uthresh))

  thresh_img = drop_img_dim(thresh_img)
  thresh = checkimg(thresh_img)
  rm(thresh_img)
  ## need to do rep.value = 0 because template is like that.

  #############################
  # Removing Neck
  #############################
  template.file =
    system.file("scct_unsmooth_SS_0.01.nii.gz",
                package = "ichseg")
  template.mask =
    system.file("scct_unsmooth_SS_0.01_Mask.nii.gz",
                package = "ichseg")
  if (remove.neck) {
    if (verbose) {
      message(paste0("# Removing Neck\n"))
    }
    L = list(
      file = thresh,
      template.file = template.file,
      template.mask = template.mask,
      rep.value = 0,
      verbose = verbose,
      ret_mask = TRUE,
      swapdim = TRUE,
      ...)
    remover = match.arg(remover)
    neck_mask = do.call(remover, args = L)
    # neck_mask = remove_neck(thresh,
    #                         rep.value = 0,
    #                         template.file = template.file,
    #                         template.mask = template.mask,
    #                         ret_mask = TRUE,
    #                         swapdim = TRUE,
    #                         verbose = verbose,
    #                         ...)
  } else {
    neck_mask = niftiarr(img, array(1, dim = dim(img)))
  }
  add_value = 0
  if (mask_to_background) {
    add_value = 1024
  }
  noneck = check_nifti(img)
  noneck = mask_img(noneck + add_value, neck_mask) - add_value
  noneck = drop_img_dim(noneck)
  noneck = checkimg(noneck)

  rm(neck_mask)

  if (verbose) {
    message("# Registering to Template")
  }
  template_with_skull =
    system.file("scct_unsmooth.nii.gz",
                package = "ichseg")
  reg = registration(
    noneck,
    template.file = template_with_skull,
    typeofTransform = "Rigid",
    interpolator = "NearestNeighbor",
    correct = FALSE,
    verbose = verbose)
  original_noneck = noneck
  noneck = reg$outfile

  #############################
  # Skull Stripping no-neck image
  #############################
  if (verbose) {
    message(paste0("# Skull Stripping in Template Space"))
  }
  template_ss = CT_Skull_Strip(
    noneck, outfile = tfile, retimg = TRUE,
    maskfile = maskfile,
    lthresh = lthresh, uthresh = uthresh,
    verbose = verbose,
    keepmask = TRUE, reorient = reorient, ...)
  template_ssmask = readnii(maskfile,
                            reorient = reorient)

  if (verbose) {
    message("# Inverting Transformation")
  }
  ssmask = ants_apply_transforms(
    fixed = original_noneck,
    moving = template_ssmask,
    transformlist = reg$fwdtransforms,
    whichtoinvert = 1,
    interpolator = "NearestNeighbor")
  ssmask = neurobase::copyNIfTIHeader(img = img, arr = ssmask)
  writenii(ssmask, maskfile)
  rm(template_ssmask)

  ss = mask_img(img + add_value, ssmask) - add_value
  ss = drop_img_dim(ss)

  writenii(ss, filename = outfile)

  #############################
  # Removing mask if keepmask = FALSE
  #############################
  if (!keepmask) {
    if (verbose) {
      message("# Removing Mask File\n")
    }
    maskfile = nii.stub(maskfile)
    ext = get.imgext()
    maskfile = paste0(maskfile, ext)
    file.remove(maskfile)
  }
  if (!retimg) {
    outfile = nii.stub(outfile)
    ext = get.imgext()
    outfile = paste0(outfile, ext)
    ss = outfile
  }
  return(ss)
}

#' @rdname CT_Skull_Strip_robust
#' @export
#' @param smooth_before_threshold Should the image be smoothed before
#' thresholding?  This can be useful for bone-window scans.
#' @param add_1024 Adding 1024 to the image *before* running the skull
#' stripping.  The values are subtracted after.  This has interplay with
#' \code{mask_to_background}
CT_Skull_Strip_smooth = function(
  img, ...,
  smooth_before_threshold = TRUE,
  smooth.factor = 1,
  remove.neck = TRUE,
  remover = c("remove_neck", "double_remove_neck"),
  recog = FALSE,
  nvoxels = 0,
  add_1024 = FALSE) {
  # grab all the args
  args = list(...)
  args$smooth_before_threshold = smooth_before_threshold
  args$smooth.factor = smooth.factor
  args$remove.neck = remove.neck
  args$remover = match.arg(remover)
  args$recog = recog
  args$nvoxels = nvoxels

  # read in the image, and add 1024, adjust the lthresh/uthresh
  if (add_1024) {
    val = 1024
    img = check_nifti(img, allow.array = TRUE)
    img = img + val
    args$lthresh = args$lthresh + val
    args$uthresh = args$uthresh + val
    mb = args$mask_to_background
    if (!is.null(mb)) {
      if (mb) {
        warning(paste0("add_1024 is TRUE, but mask_to_background ",
                       "TRUE, setting to FALSE"))
      }
    } else {
      mb = FALSE
    }
    # now 0 is truly 0
    args$mask_to_background = FALSE
  }
  args$img = img

  retimg = args$retimg
  # always return an image
  if (is.null(retimg)) {
    retimg = TRUE
  }
  # make sure you return an image and now we know for sure the outfile
  if (add_1024) {
    args$retimg = TRUE
    outfile = args$outfile
    if (is.null(outfile)) {
      outfile = tempfile(fileext = ".nii.gz")
    }
    args$outfile = outfile
  }
  # run the ss
  ss = do.call(CT_Skull_Strip_robust, args = args)
  # read in the mask
  if (add_1024) {
    mask = sub("[.]nii", "_Mask.nii", outfile)
    mask = readnii(mask)
    # remove the 1024
    if (mb) {
      # background values are -1024
      ss = mask_img(img, mask) - val
    } else {
      # background values are 0
      ss = mask_img(img - val, mask)
    }
    writenii(ss, outfile)
  }
  if (retimg) {
    return(ss)
  } else {
    return(outfile)
  }
}