#' @rdname ear_mask
#' @title Create Mask of the Ears
#' @aliases ct_ear_mask,mri_ear_mask
#' @description Creates a rough mask of the ear from a head scan
#'
#' @param file File for face masking - either filename or class nifti
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
#' @param extend_left after transformation, should the mask be extended to the
#' left of the image to ensure all left ear has been removed?
#' @param extend_right after transformation, should the mask be extended to the
#' right of the image to ensure all right ear has been removed?
#' @param typeofTransform Transformation for template to image, passed to
#' \code{\link{ants_regwrite}}.
#' @param swapdim Should the dimensions be swapped before registration,
#' and then reset after
#' @param skull_strip Should the data require skull stripping if
#' no mask is provided?
#' @param verbose Print out diagnostic messages
#' @param ... arguments passed to \code{\link{CT_Skull_Stripper}}
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
ct_ear_mask = function(
  file,
  skull_strip = TRUE,
  mask = NULL,
  robust = TRUE,
  template.file =
    system.file(
      "scct_unsmooth_SS_0.01.nii.gz",
      package = "ichseg"),
  template.ear_mask = NULL,
  template.left_ear_inds = list(170:180, 60:110, 1:60),
  template.right_ear_inds = list(1:10, 60:110, 1:60),
  extend_left = TRUE,
  extend_right = TRUE,
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

  template.left_ear_mask = .make_template_mask(
    template.file = template.file,
    template.mask = template.ear_mask,
    template.inds = template.left_ear_inds)
  template.right_ear_mask = .make_template_mask(
    template.file = template.file,
    template.mask = template.ear_mask,
    template.inds = template.right_ear_inds)

  template.ear_mask = template.left_ear_mask | template.right_ear_mask

  L = .mask_reg(file = file,
                mask = mask,
                verbose = verbose,
                swapdim = swapdim,
                template.file = template.file,
                typeofTransform = typeofTransform,
                template.mask = template.ear_mask)

  mask_trans = L$mask_trans
  img = L$img


  ######################################
  # Applying the mask to the image
  ######################################
  mask_trans = mask_trans > 0.5
  any_in_mask = any(mask_trans)
  ind = which(mask_trans, arr.ind = TRUE)
  dimg = dim(img)
  xdim = dimg[1]
  midx = xdim/2
  if (extend_left || extend_right) {
    if (any_in_mask) {
      df = data.frame(ind)
      df$left = c("left", "right")[ (df$dim1 > midx) + 1]
      ss = split(df, df$left)
      ss = lapply(ss, function(x) {
        x$left = NULL
        x
      })
      inds = NULL
      ind = ss$left
      if (!is.null(ind) && nrow(ind) > 0) {
        xs = seq(1, max(ind$dim1)) # for left
        ys = unique(ind$dim2)
        zs = unique(ind$dim3)
        run_inds = expand.grid(xs, ys, zs)
        run_inds = as.matrix(run_inds)
      }
      inds = rbind(inds, run_inds)

      ind = ss$right
      if (!is.null(ind) && nrow(ind) > 0) {
        xs = seq(min(ind$dim1), xdim) # for right
        ys = unique(ind$dim2)
        zs = unique(ind$dim3)
        run_inds = expand.grid(xs, ys, zs)
        run_inds = as.matrix(run_inds)
      }
      inds = rbind(inds, run_inds)

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



#' @rdname ear_mask
#' @export
mri_ear_mask = function(
  ...,
  skull_strip = TRUE,
  mask = NULL,
  robust = FALSE,
  template.file = mni_fname(brain = skull_strip)
){

  L = list(...)
  L$robust = robust
  L$skull_strip = skull_strip
  L$mask = mask
  L$template.file = template.file

  if (is.null(mask) & skull_strip) {
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
  res = do.call("ct_ear_mask", args = L)
  return(res)
}
