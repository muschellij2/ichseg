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
#' @importFrom extrantsr remove_neck
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
  smooth.factor = 2,
  verbose=TRUE,
  opts = NULL,
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


  if (verbose){
    message(paste0("# Thresholding Image to ",
               lthresh, "-", uthresh, "\n"))
  }


  ### Working on maskfile
  if (is.null(maskfile)){
    maskfile = nii.stub(outfile)
    maskfile = paste0(outfile, "_Mask")
  }
  maskfile = nii.stub(maskfile)
  stopifnot(inherits(maskfile, "character"))

  #############################
  # Threshold Image
  #############################
  img = check_nifti(img, reorient = reorient)
  thresh = niftiarr(img, img * (img > lthresh & img < uthresh))

  thresh = drop_img_dim(thresh)
  ## need to do rep.value = 0 because template is like that.

  #############################
  # Removing Neck
  #############################
  if (remove.neck) {
    if (verbose) {
      message(paste0("# Removing Neck\n"))
    }
    neck_mask = remove_neck(thresh,
                            rep.value = 0,
                            template.file = template.file,
                            template.mask = template.mask,
                            ret_mask = TRUE,
                            swapdim = TRUE,
                            ...)
  } else {
    neck_mask = niftiarr(img, array(1, dim = dim(img)))
  }
  noneck = mask_img(img, neck_mask)
  noneck = drop_img_dim(noneck)

  #############################
  # Skull Stripping no-neck image
  #############################
  if (verbose) {
    message(paste0("# Skull Stripping for COG\n"))
  }
  ss = CT_Skull_Strip(noneck, outfile = outfile, retimg = TRUE,
                      maskfile = maskfile,
                      lthresh = lthresh, uthresh = uthresh,
                      opts = paste("-f ", int,
                                   ifelse(verbose, " -v", "")),
                      keepmask = TRUE, reorient = reorient, ...)

  ssmask = readNIfTI(maskfile,
                     reorient = reorient)

  #############################
  # Setting new center of gravity and rerunning with smoothness factor
  #############################
  cog = cog(ssmask,
            ceil=TRUE)
  if (verbose){
    message(paste0("# Skull Stripping with new cog\n"))
  }
  ss = CT_Skull_Strip(noneck, outfile = outfile, retimg = TRUE,
                      opts = paste("-f ", int,
                                   ifelse(verbose, "-v", ""),
                                   "-w ", smooth.factor,
                                   paste(c("-c", cog),
                                         collapse=" ")),
                      maskfile = maskfile,
                      keepmask = TRUE,
                      reorient = reorient,
                      ...)

  ssmask = readNIfTI(maskfile,
                     reorient = reorient)

  #############################
  # Filling the mask
  #############################
  if (verbose) {
    message(paste0("# Filling Holes \n"))
  }
  if (nvoxels > 0){
    ssmask = dil_ero(ssmask,
                    retimg = TRUE,
                    nvoxels = nvoxels)
  }

  ss = mask_img(img, ssmask)
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

  return(ss)
}