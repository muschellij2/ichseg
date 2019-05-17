#' Create CT Head Mask
#'
#' @param img CT image of class \code{\link{nifti}}
#' @param fill_size Size of fill (in voxels),
#' passed to \code{\link[extrantsr]{filler}}
#'
#' @return An object of class \code{\link{nifti}}
#' @export
create_ct_head_mask = function(img, fill_size = 11) {
  img = neurobase::check_nifti(img)
  head_mask = img >= -100 & img <= 1000
  cc = extrantsr::largest_component(head_mask)
  filled = extrantsr::filler(cc, fill_size = fill_size)
  return(filled)
}

#' @export
#' @rdname create_ct_head_mask
segment_ct_head = function(img, fill_size = 11) {
  img = neurobase::check_nifti(img)
  mask = create_ct_head_mask(img, fill_size = fill_size)
  head_img = neurobase::mask_img(img + 1024, mask)
  head_img = head_img - 1024
  return(img)
}