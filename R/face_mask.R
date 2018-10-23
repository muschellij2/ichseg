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
#' @param extend_mask after transformation, should the mask be extended to the
#' front of the image to ensure all face has been removed?
#' @param typeofTransform Transformation for template to image, passed to
#' \code{\link{ants_regwrite}}.
#' @param swapdim Should the dimensions be swapped before registration,
#' and then reset after
#' @param skull_strip Should the data require skull stripping?
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
#'    robust = FALSE
#'    )
#'  img = readnii(file)
#'  rimg = randomize_mask(img, mask = face)
#' }
#'
ct_face_mask <- function(
  file,
  skull_strip = TRUE,
  mask = NULL,
  robust = TRUE,
  template.file =
    system.file(
      ifelse(skull_strip,
             "scct_unsmooth_SS_0.01.nii.gz",
             "scct_unsmooth.nii.gz"),
      package = "ichseg"),
  template.face_mask = NULL,
  template.face_mask_inds = list(50:130, 170:217, 1:15),
  extend_mask = TRUE,
  typeofTransform = "Affine",
  swapdim = TRUE,
  verbose = TRUE,
  ...){

  if (skull_strip) {
    mask = .make_ss_mask(file = file,
                         mask = mask,
                         verbose = verbose,
                         robust = robust,
                         template.file = template.file, ...)
  }

  template.face_mask = .make_template_mask(
    template.file = template.file,
    template.mask = template.face_mask,
    template.inds = template.face_mask_inds)

  L = .mask_reg(file = file,
                mask = mask,
                verbose = verbose,
                swapdim = swapdim,
                template.file = template.file,
                typeofTransform = typeofTransform,
                template.mask = template.face_mask)

  mask_trans = L$mask_trans
  img = L$img


  ######################################
  # Applying the mask to the image
  ######################################
  mask_trans = mask_trans > 0.5
  any_in_mask = any(mask_trans)
  ind = which(mask_trans, arr.ind = TRUE)
  if (extend_mask) {
    if (any_in_mask) {
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
    } else {
      warning("No registered object in mask found - cannot extend!")
      newimg = mask_trans
    }
  } else {
    newimg = mask_trans
  }


  if (swapdim) {
    if (verbose) {
      message(paste0("# Swapping Dimensions Back\n"))
    }
    sorient = L$sorient
    ori = L$ori
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
  skull_strip = TRUE,
  mask = NULL,
  robust = FALSE,
  template.file = mni_fname(brain = skull_strip)
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
    uthresh = L$uthresh
    if (is.null(uthresh)) {
      file = L$file
      if (is.null(file)) {
        file = L[[1]]
      }
      uthresh = fslr::fslmax(file)
      L$uthresh = uthresh
    }

    L = func(L, "presmooth", FALSE)
    L = func(L, "remask", FALSE)
    L = func(L, "inskull_mesh", FALSE)
    L = func(L, "opts", "-v")
  }
  res = do.call("ct_face_mask", args = L)
  return(res)
}
