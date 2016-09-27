#' @title Create Mask of the Face
#'
#' @description Creates a rough mask of the face of a brain image
#'
#' @param file File for face masking - either filename or class nifti
#' @param mask file or \code{nifti} to mask the \code{file}
#' @param robust If \code{mask = NULL}, then \code{robust} is
#' passed to \code{\link{CT_Skull_Stripper}}
#' @param template.file Template to warp to original image space
#' @param template.mask Mask of template to use as rough face mask.  If
#' \code{template.file} is not specified, \code{template.mask_inds}
#' must be
#' @param template.mask_inds List of length 3 for indices of
#' \code{template.file} to indicate the mask.
#' @param typeofTransform Transformation for template to image, passed to
#' \code{\link{ants_regwrite}}.
#' @param swapdim Should the dimensions be swapped before registration,
#' and then reset after
#' @param verbose Print out diagnostic messages
#' @param ... not used
#' @export
#' @return Object of class nifti
#' @importFrom neurobase check_mask_fail
#' @importFrom extrantsr rpi_orient reverse_rpi_orient
face_mask <- function(
  file,
  mask = NULL,
  robust = TRUE,
  template.file =
    system.file(
      "scct_unsmooth_SS_0.01.nii.gz",
      package = "ichseg"),
  template.mask = NULL,
  template.mask_inds = list(50:130, 170:217, 1:15),
  typeofTransform = "Affine",
  swapdim = TRUE,
  verbose = TRUE,
  ...){


  if (is.null(mask)) {
    if (verbose) {
      message(paste0("# Skull Stripping \n"))
    }
    mask = tempfile(fileext = ".nii.gz")
    ss = CT_Skull_Stripper(file,
                           maskfile = mask,
                           verbose = verbose,
                           robust = robust)
    rm(ss)
  }
  mask = check_nifti(mask)
  check_mask_fail(mask)

  if (is.null(template.mask) &&
      is.null(template.mask_inds)) {
    stop("Need template.mask or template.mask_inds")
  }

  template.file = checkimg(template.file)
  if (is.null(template.mask)) {
    nim = check_nifti(template.file)

    ###############################
    # using inds
    ###############################
    template.mask = niftiarr(nim, 0)
    inds = expand.grid(template.mask_inds)
    inds = as.matrix(inds)
    template.mask[inds] = 1
    # template.mask[50:130, 170:217, 1:15] = 1
    template.mask = cal_img(template.mask)
  }
  check_mask_fail(template.mask)

  file = check_nifti(file)
  # ofile = tempfile(fileext = '.nii.gz')
  img = mask_img(file, mask)

  if (swapdim) {
    if (verbose) {
      message(paste0("# Swapping Dimensions \n"))
    }
    L = rpi_orient(file)
    file = L$img
    sorient = L$orientation
    ori = L$convention
  }

  if (verbose) {
    message("# Template Registration to image\n")
  }

  # Registration
  rev_reg = registration(
    filename = template.file,
    skull_strip = FALSE,
    correct = FALSE,
    template.file = img,
    typeofTransform = "Affine")

  if (verbose) {
    message(paste0("# Applying Transforms to template mask \n"))
  }
  # Apply Transform to Mask image
  mask_trans = ants_apply_transforms(
    fixed = img,
    moving = template.mask,
    transformlist = rev_reg$fwdtransforms,
    interpolator = "nearestNeighbor")

  if (verbose) {
    message(paste0("# Applying Mask to Original Image \n"))
  }
  ######################################
  # Applying the mask to the image
  ######################################
  ind = which(mask_trans > 0.5, arr.ind = TRUE)
  minz = ceiling(mean(ind[,"dim3"]))
  zs = seq(minz)
  miny = min(ind[,"dim2"])
  ys = seq(miny, dim(mask_trans)[2])
  xs = unique(ind[,"dim1"])
  inds = expand.grid(xs, ys, zs)
  inds = as.matrix(inds)
  newimg = niftiarr(mask, 0)
  newimg[inds] = 1
  newimg = cal_img(newimg)


  if (swapdim) {
    if (verbose) {
      message(paste0("# Swapping Dimensions Back\n"))
    }
    newimg = reverse_rpi_orient(
      file = newimg,
      convention = ori,
      orientation = sorient, verbose = verbose)
  }
  return(newimg)
}
