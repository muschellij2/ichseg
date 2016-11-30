#' @rdname ear_mask
#' @title Create Mask of the Ears
#' @aliases ct_ear_mask,mri_ear_mask
#' @description Creates a rough mask of the ear from a head scan
#'
#' @param file File for ear masking - either filename or class nifti
#' @param mask file or \code{nifti} to mask the \code{file}
#' @param robust If \code{mask = NULL}, then \code{robust} is
#' passed to \code{\link{CT_Skull_Stripper}}
#' @param template.file Template to warp to original image space
#' @param template.ear_mask Mask of template to use as rough ear mask.  If
#' \code{template.file} is not specified, \code{template.left_ear_inds} and
#' \code{template.right_ear_inds}
#' must be
#' @param template.left_ear_inds List of length 3 for indices of
#' \code{template.file} to indicate the left-ear mask.
#' @param template.right_ear_inds List of length 3 for indices of
#' \code{template.file} to indicate the right-ear mask.
#' @param ... arguments passed to \code{\link{ct_face_mask}}
#' @export
#' @return Object of class nifti
#' @examples \dontrun{
#' file = "~/Desktop/Desktop/scratch/100-318_20070723_0957_CT_3_CT_Head-.nii.gz"
#' mask = NULL
#' robust = FALSE
#' face = ct_ear_mask(
#'    file = file,
#'    robust = FALSE,
#'     template.mask = system.file("scct_unsmooth_SS_0.01_Mask.nii.gz",
#'                   package = "ichseg")
#'    )
#'  img = readnii(file)
#'  rimg = randomize_mask(img, mask = face)
#' }
#'#'
ct_ear_mask = function(
  file,
  template.file =
    system.file(
      "scct_unsmooth_SS_0.01.nii.gz",
      package = "ichseg"),
  template.ear_mask = NULL,
  template.left_ear_inds = list(170:180, 60:110, 1:60),
  template.right_ear_inds = list(1:10, 60:110, 1:60),
  ...
){

  if (is.null(template.ear_mask) &&
      (is.null(template.left_ear_inds) ||
       is.null(template.right_ear_inds))
      ){
    stop("Need template.ear_mask or template.left/right_ear_inds")
  }

  template.file = checkimg(template.file)
  if (is.null(template.ear_mask)) {
    nim = check_nifti(template.file)

    ###############################
    # using inds
    ###############################
    template.ear_mask = niftiarr(nim, 0)
    # right ear
    inds = expand.grid(template.right_ear_inds)
    inds = as.matrix(inds)
    template.ear_mask[inds] = 1

    # left ear
    inds = expand.grid(template.left_ear_inds)
    inds = as.matrix(inds)
    template.ear_mask[inds] = 1
    template.ear_mask = cal_img(template.ear_mask)
  }
  check_mask_fail(template.ear_mask)

  res = ct_face_mask(
      file = file,
      template.file = template.file,
      template.face_mask = template.ear_mask,
      template.face_mask_inds = NULL,
      ...)

  return(res)

}
