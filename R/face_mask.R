#' @rdname face_mask
#' @title Create Mask of the Face
#' @aliases ct_face_mask,mri_face_mask
#' @description Creates a rough mask of the face from a head scan
#'
#' @param file File for face masking - either filename or class nifti
#' @param mask file or \code{nifti} to mask the \code{file}
#' @param robust If \code{mask = NULL}, then \code{robust} is
#' passed to \code{\link{CT_Skull_Stripper}}
#' @param template.file Template to warp to original image space
#' @param template.face_mask Mask of template to use as rough face mask.  If
#' \code{template.file} is not specified, \code{template.face_mask_inds}
#' must be
#' @param template.face_mask_inds List of length 3 for indices of
#' \code{template.file} to indicate the mask.
#' @param typeofTransform Transformation for template to image, passed to
#' \code{\link{ants_regwrite}}.
#' @param swapdim Should the dimensions be swapped before registration,
#' and then reset after
#' @param verbose Print out diagnostic messages
#' @param ... arguments passed to \code{\link{CT_Skull_Stripper}}
#' @export
#' @return Object of class nifti
#' @importFrom neurobase check_mask_fail
#' @importFrom fslr rpi_orient reverse_rpi_orient
#' @examples \dontrun{
#' file = "~/Desktop/Desktop/scratch/100-318_20070723_0957_CT_3_CT_Head-.nii.gz"
#' mask = NULL
#' robust = FALSE
#' face = ct_face_mask(
#'    file = file,
#'    robust = FALSE,
#'     template.mask = system.file("scct_unsmooth_SS_0.01_Mask.nii.gz",
#'                   package = "ichseg")
#'    )
#'  img = readnii(file)
#'  rimg = randomize_mask(img, mask = face)
#' }
#'#'
ct_face_mask <- function(
  file,
  mask = NULL,
  robust = TRUE,
  template.file =
    system.file(
      "scct_unsmooth_SS_0.01.nii.gz",
      package = "ichseg"),
  template.face_mask = NULL,
  template.face_mask_inds = list(50:130, 170:217, 1:15),
  typeofTransform = "Affine",
  swapdim = TRUE,
  verbose = TRUE,
  ...){


  if (is.null(mask)) {
    if (verbose) {
      message(paste0("# Skull Stripping \n"))
    }
    mask = tempfile(fileext = ".nii.gz")
    ss = CT_Skull_Stripper(
      img = file,
      maskfile = mask,
      verbose = verbose,
      robust = robust,
      template.file = template.file,
      ...)
    rm(ss)
  }
  mask = check_nifti(mask)
  check_mask_fail(mask)

  if (is.null(template.face_mask) &&
      is.null(template.face_mask_inds)) {
    stop("Need template.face_mask or template.face_mask_inds")
  }

  template.file = checkimg(template.file)
  if (is.null(template.face_mask)) {
    nim = check_nifti(template.file)

    ###############################
    # using inds
    ###############################
    template.face_mask = niftiarr(nim, 0)
    inds = expand.grid(template.face_mask_inds)
    inds = as.matrix(inds)
    template.face_mask[inds] = 1
    # template.face_mask[50:130, 170:217, 1:15] = 1
    template.face_mask = cal_img(template.face_mask)
  }
  check_mask_fail(template.face_mask)

  file = check_nifti(file)
  # ofile = tempfile(fileext = '.nii.gz')
  img = mask_img(file, mask)
  rm(list = c("file", "mask"))
  if (swapdim) {
    if (verbose) {
      message(paste0("# Swapping Dimensions \n"))
    }
    L = rpi_orient(img)
    img = L$img
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
    typeofTransform = typeofTransform,
    verbose = verbose)

  if (verbose) {
    message(paste0("# Applying Transforms to template mask \n"))
  }

  # Apply Transform to Mask image
  mask_trans = ants_apply_transforms(
    fixed = img,
    moving = template.face_mask,
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

  newimg = niftiarr(img, 0)
  newimg[inds] = 1
  newimg = cal_img(newimg)


  if (swapdim) {
    if (verbose) {
      message(paste0("# Swapping Dimensions Back\n"))
    }
    newimg = reverse_rpi_orient(
      file = newimg,
      convention = ori,
      orientation = sorient,
      verbose = verbose)
  }
  return(newimg)
}


#' @rdname face_mask
#' @export
#' @importFrom fslr mni_fname
#' @examples \dontrun{
#' library(fslr)
#' library(extrantsr)
#' mri = "~/Desktop/Desktop/scratch/SUBJ0001-01-MPRAGE.nii.gz"
#'
#' template.file = mni_fname(brain = TRUE)
#' tmask = mni_fname(brain = TRUE, mask = TRUE)
#'
#' template.face_mask_inds = list(50:130, 170:217, 1:15)
#' brain = fslbet_robust(mri,
#' remove.neck = TRUE,
#' remover = "double_remove_neck",
#' template.file = template.file,
#' template.mask = tmask)
#' mask = brain > 0
#' img = brain
#' template.face_mask = NULL
#' verbose = TRUE
#' face = mri_face_mask(
#'    file = img,
#'    mask = mask,
#'    template.file = template.file
#'    )
#' }
#'
mri_face_mask <- function(
  ...,
  mask = NULL,
  robust = FALSE,
  template.file = mni_fname(brain = TRUE)
){


  L = list(...)
  L$robust = robust
  L$mask = mask
  L$template.file = template.file

  if (is.null(mask)) {
    func = function(L, arg, opt) {
      nL = names(L)
      if (!arg %in% nL) {
        L[arg] = opt
      }
      return(L)
    }
    L = func(L, "presmooth", FALSE)
    L = func(L, "remask", FALSE)
    L = func(L, "inskull_mesh", FALSE)
    L = func(L, "opts", "-v")
  }
  res = do.call("ct_face_mask", args = L)
  return(res)
}
