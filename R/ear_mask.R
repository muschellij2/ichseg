#' @rdname ear_mask
#' @title Create Mask of the Ears
#' @aliases ct_ear_mask,mri_ear_mask
#' @description Creates a rough mask of the ear from a head scan
#'
#' @param file File for ear masking - either filename or class nifti
#' @param type is it a CT or MRI?  Will determine if \code{\link{ct_face_mask}} or
#' \code{\link{mri_face_mask}} is called
#' @param template.file Template to warp to original image space
#' @param template.ear_mask Mask of template to use as rough ear mask.  If
#' \code{template.file} is not specified, \code{template.left_ear_inds} and
#' \code{template.right_ear_inds}
#' must be
#' @param template.left_ear_inds List of length 3 for indices of
#' \code{template.file} to indicate the left-ear mask.
#' @param template.right_ear_inds List of length 3 for indices of
#' \code{template.file} to indicate the right-ear mask.
#' @param extend_mask passed to \code{\link{ct_face_mask}} or
#' \code{\link{mri_face_mask}}, but should be \code{FALSE} otherwise the
#' mask will be extended forward
#' @param ... arguments passed to \code{\link{ct_face_mask}} or
#' \code{\link{mri_face_mask}}
#' @export
#' @return Object of class nifti
#' @examples \dontrun{
#' file = "~/Desktop/Desktop/scratch/100-318_20070723_0957_CT_3_CT_Head-.nii.gz"
#' mask = NULL
#' robust = FALSE
#' ears = ct_ear_mask(
#'    file = file,
#'    robust = FALSE
#'  )
#'  img = readnii(file)
#'  rimg = randomize_mask(img, mask = ears)
#' }
ear_mask = function(
  file,
  type = c("ct", "mri"),
  template.file =
    system.file(
      "scct_unsmooth_SS_0.01.nii.gz",
      package = "ichseg"),
  template.ear_mask = NULL,
  template.left_ear_inds = list(170:180, 60:110, 1:60),
  template.right_ear_inds = list(1:10, 60:110, 1:60),
  extend_mask = FALSE,
  ...
){

  if (is.null(template.ear_mask) &&
      (is.null(template.left_ear_inds) ||
       is.null(template.right_ear_inds))
  ) {
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

  type = match.arg(type)
  func = paste0(type, "_face_mask")
  L = list(
    file = file,
    template.file = template.file,
    template.face_mask = template.ear_mask,
    template.face_mask_inds = NULL,
    extend_mask = extend_mask,
    ... = ...)
  res = do.call(func, args = L)

  return(res)

}

#' @rdname ear_mask
#' @export
ct_ear_mask = function(...) {
  ear_mask(..., type = "ct")
}

#' @rdname ear_mask
#' @export
mri_ear_mask = function(
  ...,
  template.file = mni_fname(brain = TRUE)
) {
  ear_mask(...,
           template.file = template.file,
           type = "mri")
}
