
.make_ss_mask = function(file, mask, verbose, robust, template.file, ...) {
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
  return(mask)
}

.make_template_mask = function(template.file,
                               template.mask,
                               template.inds) {
  if (is.null(template.mask) &&
      is.null(template.inds)) {
    stop("Need template.face_mask or template.face_mask_inds")
  }

  template.file = checkimg(template.file)
  if (is.null(template.mask)) {
    nim = check_nifti(template.file)

    ###############################
    # using inds
    ###############################
    template.mask = niftiarr(nim, 0)
    inds = expand.grid(template.inds)
    inds = as.matrix(inds)
    template.mask[inds] = 1
    # template.face_mask[50:130, 170:217, 1:15] = 1
    template.mask = cal_img(template.mask)
  }
  template.mask = check_nifti(template.mask, allow.array = TRUE)
  check_mask_fail(template.mask)
  return(template.mask)
}



.mask_reg = function(file, mask, verbose, swapdim,
                     template.file, typeofTransform,
                     template.mask) {
  file = check_nifti(file)
  # ofile = tempfile(fileext = '.nii.gz')
  if (!is.null(mask)) {
    img = mask_img(file, mask)
  } else {
    img = file
  }
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
    moving = template.mask,
    transformlist = rev_reg$fwdtransforms,
    interpolator = "nearestNeighbor")

  if (verbose) {
    message(paste0("# Applying Mask to Original Image \n"))
  }
  L = list(mask_trans = mask_trans,
           img = img)
  L$fwdtransforms = rev_reg$fwdtransforms
  if (swapdim) {
    L$sorient = sorient
    L$ori = ori
  }
  return(L)
}