#' @title Fill image holes with dilation then erosion
#' @description This function calls \code{mean_image} to dilate an image, then calls
#' it again to erode it.
#' @param file (character) filename of image to be filled
#' @param outfile (character) name of resultant filled file
#' @param nvoxels (integer) Number of voxels to smooth over, creates vxvxv box.
#' @param zeropad (logical) Perform \code{zero_pad} before running.
#' @param remove.ends (logical) Remove top and bottom dilation.
#' @param tol (double) Tolerance for thresholding after \code{\link{fft}}.
#' @param refill (logical) Run \code{\link{fslfill}} after dilation/erosion.
#' @param retimg (logical) return image of class nifti
#' @param reorient (logical) If retimg, should file be reoriented when read in?
#' Passed to \code{\link{readNIfTI}}.
#' @param intern (logical) pass to \code{\link{system}}
#' @param verbose (logical) print out command before running
#' @param ... additional arguments passed to \code{\link{readNIfTI}}.
#' @return character or logical depending on intern
#' @importFrom neurobase zero_pad
#' @note This function binarizes the image before running.
#' @export
dil_ero = function(file,
                    outfile = NULL,
                    nvoxels = 3,
                    zeropad = TRUE,
                    remove.ends = FALSE,
                    tol = .Machine$double.eps^0.5,
                    refill = TRUE,
                    retimg = FALSE,
                    reorient = FALSE,
                    intern=TRUE, verbose = TRUE,
                    ...){
  have.outfile = TRUE

  if (retimg){
    if (is.null(outfile)) {
      outfile = tempfile()
      have.outfile = FALSE
    }
  } else {
    stopifnot(!is.null(outfile))
  }

  file = check_nifti(file, reorient = reorient)
  bin = niftiarr(file, file > 0)
  dimg = dim(bin)
  ##### should make for all max
  ind = which(bin >0, arr.ind=TRUE)
  ind = ind[ (ind[, "dim3"] %in% c(1, dimg[3])) |
               (ind[, "dim1"] %in% c(1, dimg[1])) |
               (ind[, "dim2"] %in% c(1, dimg[2])) ,]
  nind = nrow(ind)

  img = bin
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }
  stopifnot(length(nvoxels) %in% c(1, 3))
  if (length(nvoxels) == 1){
    stopifnot(is.wholenumber(nvoxels))
    n <- length((-nvoxels):nvoxels)
    kdim = c(n, n, n)
  }
  if (length(nvoxels) == 3){
    kdim = (nvoxels*2) + 1
    stopifnot(is.wholenumber(kdim))
  }

  if (zeropad) img = zero_pad(img, kdim = kdim)
  dilated = mean_image(img, nvoxels = nvoxels, verbose = verbose) > tol
  dil = mean_image(1-dilated, nvoxels = nvoxels, verbose = verbose)
  dil = dil > tol
  dil = 1 - dil
  if (zeropad) dil = zero_pad(dil, kdim = kdim, invert=TRUE)
  dil = niftiarr(bin, dil)

  if (remove.ends) {
    #### making the ends correct - boundary problem
    dil@.Data[,,1] = array(0, dim=dimg[c(1,2)])
    dil@.Data[,,dimg[3]] = array(0, dim=dimg[c(1,2)])

    dil@.Data[,1,] = array(0, dim=dimg[c(1,3)])
    dil@.Data[,dimg[2],] = array(0, dim=dimg[c(1,3)])

    dil@.Data[1,,] = array(0, dim=dimg[c(2,3)])
    dil@.Data[dimg[1],,] = array(0, dim=dimg[c(2,3)])

    if (nind >0 ){
      dil@.Data[ ind ] = 1
    }
  }
  if (refill) {
    dil = fslfill(file = dil,
                  retimg=TRUE,
                  intern=intern, verbose = verbose)
  }
  dil = cal_img(dil)
  if (have.outfile){
    gzipped = grepl("gz$", get.imgext())
    writenii(dil, filename = outfile, gzipped = gzipped)
  }
  if (retimg){
    return(dil)
  }
  return(outfile)
}
