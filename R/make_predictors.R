#' @title Make CT Predictors
#'
#' @description Create a set of predictors for ICH segmentation
#' for CT
#' @param img Filename of image intensities
#' @param mask Filename of brain mask
#' @param roi Filename of ROI for Y
#' @param nvoxels Voxel neighborhood
#' @param moments Moments of neighborhood to take
#' @param center Center the moments
#' @param lthresh Lower threshold for neighborhood setting
#' @param uthresh Upper threshold for neighborhood setting
#' @param sigmas Sigma values for Gaussian smoothed images (in mm)
#' @param save_imgs Logical to save all images that are created as
#' predictors
#' @param outdir Output directory of saved images, needs to be set
#' if \code{save_imgs = TRUE}
#' @param stub Basename to write image names if \code{save_imgs = TRUE}
#' @param overwrite If \code{save_imgs} is \code{TRUE},
#' then should
#' the files be overwritten?  If not, then files will be read
#' instead instead of code being re-run.
#' @param template.file Template to register to (CT Template)
#' @param mean.img Mean image in template space for z-scoring
#' @param sd.img SD image in template space for z-scoring
#' @param zscore.typeofTransform type of transform for z-scoring
#' @param zscore.interpolator type of interpolator for z-scoring
#' @param flip.typeofTransform type of transform for flipped difference
#' @param flip.interpolator type of interpolator for flipped difference
#' @param low_thresh Threshold for forcing values to zero
#' @param verbose Logical indicator if output messages should be
#' printed
#' @param ... options passed to \code{\link{get_neighbors}}
#' @importFrom fslr readnii checkimg fslerode zscore_img
#' @importFrom oro.nifti zero_trans cal_img voxdim pixdim convert.datatype convert.bitpix
#' @importFrom extrantsr zscore_template otropos reg_flip
#' @export
#' @return List of a data.frame of Predictors and set of
#' indices to
#' keep in mask and an empty nifti object for plotting.
#' Also the number of voxels of the roi that were not in the
#' mask
make_predictors <- function(img, mask, roi = NULL,
                            nvoxels = 1, moments = 1:4,
                            center = c(FALSE, TRUE, TRUE, TRUE),
                            lthresh = 40, uthresh = 80,
                            sigmas = c(5, 10, 20),
                            save_imgs = TRUE,
                            outdir = NULL,
                            stub = NULL,
                            overwrite = FALSE,
                            template.file = system.file(
                              "scct_unsmooth_SS_0.01.nii.gz",
                              package = "ichseg"),
                            mean.img = system.file(
                              "Mean_Image.nii.gz",
                              package = "ichseg"),
                            sd.img = system.file(
                              "SD_Image.nii.gz",
                              package = "ichseg"),
                            zscore.typeofTransform = "SyN",
                            zscore.interpolator = "Linear",
                            flip.typeofTransform = "Affine",
                            flip.interpolator = "LanczosWindowedSinc",
                            low_thresh = 1e-13,
                            verbose= TRUE,
                            ...) {
  make_fname = function(addstub){
    fname = addstub
    fname = file.path(outdir, paste0(stub, "_", fname, ".nii.gz"))
    fname
  }

  write_img = function(arr, addstub){
    fname = addstub
    fname = file.path(outdir, paste0(stub, "_", fname))
    if ( !inherits( arr, "nifti")){
      mom = nim
      mom@.Data = array(arr, dim = dim(mom))
    } else {
      mom = arr
    }
    mom = cal_img(mom)
    mom = zero_trans(mom)
    mom = datatyper(mom,
                    datatype = convert.datatype()$FLOAT32,
                    bitpix = convert.bitpix()$FLOAT32)
    writenii(mom, filename = fname)
  }


  make_moment_list = function(){
    mom.exist = sapply(moments, function(moment){
      addstub = paste0("moment", moment)
      fname = make_fname(addstub)
      file.exists(fname)
    })

    ### if all data does not exist or rewrite - remake data
    if (!all(mom.exist) | overwrite) {
      if (verbose){
        message("# Running Local_Moment")
      }
      moms = local_moment(
        img,
        mask = mask,
        nvoxels = nvoxels,
        moment = moments, center = center,
        invert = FALSE,
        ...)
      if (save_imgs) {
        for (imom in seq_along(moments)) {
          moment = moments[imom]
          addstub = paste0("moment", moment)
          fname = make_fname(addstub)
          mom.img = moms[[imom]]
          write_img(mom.img, addstub)
        }
      }
    }
    ### if all data exists and not to rewrite - just read in
    if (all(mom.exist) & !overwrite) {
      moms = lapply(moments, function(moment){
        addstub = paste0("moment", moment)
        fname = make_fname(addstub)
        mom.img = readnii(fname, reorient = FALSE)
        return(mom.img)
      })
    }
    if (verbose) {
      message("# Creating Moment Matrix")
    }
    mat = matrix(NA,
                 ncol = length(moms),
                 nrow = prod(dim(mask))
                 )
    for (imom in seq_along(moms)) {
      mat[, imom] = c(moms[[imom]])
    }
    # moms = sapply(moms, c)
    rm(list = "moms"); gc(); gc();
    colnames(mat) = paste0("moment", moments)
    mat = as.data.frame(mat)
    return(mat)
  }

  img.fname = checkimg(img)
  img = check_nifti(img)

  mask.fname = checkimg(mask)
  orig.mask = mask = check_nifti(mask)

  if (!is.null(roi)) {
    #     roi.fname = roi
    roi = check_nifti(roi)
  }

  # stopifnot(class(roi) == "character")
  #   stopifnot(class(img) == "character")
  #   stopifnot(class(mask) == "character")


  if (save_imgs){
    stopifnot(!is.null(outdir))
    stopifnot(!is.null(stub))
  }
  if (is.null(outdir)){
    outdir = tempdir()
  }
  if (is.null(stub)){
    stub = nii.stub(img.fname, bn = TRUE)
  }

  if (verbose) {
    message("\n# Reading Images\n")
  }
  # img = readnii(img.fname, reorient= FALSE)
  orig.img = img
  # stub = nii.stub(img.fname, bn=TRUE)

  nim = niftiarr(img, array(NA, dim = dim(orig.img)))
  nim = datatyper(nim,
                  datatype = convert.datatype()$FLOAT32,
                  bitpix = convert.bitpix()$FLOAT32)
  dimg = dim(nim)
  vdim = voxdim(nim)

  #   if (is.null(roi)){
  #     roi = niftiarr(img, array(0, dim = dim(nim)))
  #   }

  if (verbose) {
    message("# Eroding Mask\n")
  }
  addstub = "usemask"
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite){
    mask = readnii(fname, reorient=FALSE)
  } else {
    # erode the mask
    mask = fslerode(file = mask.fname,
                    kopts = "-kernel box 3x3x1",
                    reorient = FALSE, retimg = TRUE)
    #### may add this - think about it
    #     mask = fslfill(mask, bin=TRUE, retimg = TRUE)
    mask = mask > 0
    mask = datatyper(mask)
    if (save_imgs){
      write_img(mask, addstub)
    }
  }

  orig.masked.img = mask_img(orig.img, orig.mask)
  masked.img = mask_img(orig.img, mask)
  if (!is.null(roi)){
    miss.roi = sum(mask == 0 & roi > 0 )
  } else {
    miss.roi = NULL
  }


  keep.ind = which(mask > 0)

  if (verbose) {
    message("# Getting Moments\n")
  }
  ################################################
  # Making Moment Images
  ################################################

  df = make_moment_list()
  ###########
  # Creating Skew/Kurtosis
  ###########
  df$skew = df$moment3/df$moment2 ^ {3/2}
  df$kurtosis = df$moment4/df$moment2 ^ {2}
  df$moment3 = df$moment4 = NULL

  ###########
  # Writing Skew
  ###########
  addstub = "skew"
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite) {

  } else {
    skew = niftiarr(img, df$skew)
    if (save_imgs) {
      write_img(skew, addstub)
    }
  }

  ###########
  # Writing Kurtosis
  ###########
  addstub = "kurtosis"
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite) {

  } else {
    kurtosis = niftiarr(img, df$kurtosis)
    if (save_imgs) {
      write_img(kurtosis, addstub)
    }
  }
  rm(list = c("kurtosis", "skew")); gc(); gc();

  ################################################
  # Making Percent threshold image
  ################################################
  if (verbose) {
    message(paste0("# Getting thresholded from ",
                   lthresh, " to ", uthresh, "\n"))
  }
  addstub = paste0("thresh_", lthresh, "_", uthresh)
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite){
    thresh_img = readnii(fname, reorient=FALSE)
  } else {
    thresh_img = niftiarr(img, img >= lthresh & img <= uthresh)
    if (save_imgs){
      write_img(thresh_img, addstub)
    }
  }

  df$value = c(masked.img)
  df$thresh = c(thresh_img)

  cn = colnames(df)

  ################################################
  # Making Z-score Images
  ################################################
  if (verbose) {
    message(paste0("# Getting Z-scored images\n"))
  }
  zimgs = lapply(1:3, function(i){
    addstub = paste0("zscore", i)
    fname = make_fname(addstub)
    if (file.exists(fname) & !overwrite){
      img = readnii(fname, reorient=FALSE)
    } else {
      img = zscore_img(masked.img, mask = mask, margin = i)
      if (save_imgs){
        write_img(img, addstub)
      }
    }
    return(img)
  })
  zimgs = sapply(zimgs, c)
  colnames(zimgs) = paste0("zscore", 1:3)
  zimgs = as.data.frame(zimgs)
  for (i in 1:3) {
    df[,  paste0("zscore", i)] = zimgs[, paste0("zscore", i)]
  }
  rm(list = c("zimgs")); gc(); gc();

#   wmean2 = function(img, mask, trim = 0.2){
#     x = img[ mask == 1]
#     mn = psych::winsor.mean(x, trim = trim)
#     s = psych::winsor.sd(x, trim = trim)
#     z = (x-mn)/s
#     img[mask == 1] = z
#     img[mask == 0] = 0
#     img = cal_img(img)
#     return(img)
#   }

  wmean = function(img, mask, trim = 0.2){
    x = img[ mask == 1 ]
    stopifnot(length(trim) == 1)
    stopifnot(trim > 0)
    stopifnot(trim <= 0.5)
    qtrim <- quantile(x,
                      c(trim, 0.5, 1 - trim),
                      na.rm = TRUE)
    xbot <- qtrim[1]
    xtop <- qtrim[3]

    if (trim < 0.5) {
      x[x < xbot] <- xbot
      x[x > xtop] <- xtop
    } else {
      x[!is.na(x)] <- qtrim[2]
    }

    mn = mean(x, na.rm=TRUE)
    s = sd(x, na.rm=TRUE)
    img = (img-mn)/s
    img = mask_img(img, mask)
    img = finite_img(img)
    return(img)
  }

  if (verbose) {
    message("# Getting Winsorized Image\n")
  }
  addstub = paste0("win_z")
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite){
    wmean_img = readnii(fname, reorient=FALSE)
  } else {
    wmean_img = wmean(orig.img, mask = mask, trim = 0.2)
    if (save_imgs){
      write_img(wmean_img, addstub)
    }
  }
  df$win_z = c(wmean_img)
  rm(list = c("wmean_img")); gc(); gc();

  ################################################
  # Making Percent threshold image
  ################################################
  if (verbose) {
    message(paste0("# Getting Percent thresholded from ",
                   lthresh, " to ", uthresh, "\n"))
  }
  addstub = paste0("pct_thresh_", lthresh, "_", uthresh)
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite){
    pct_thresh = readnii(fname, reorient=FALSE)
  } else {
    pct_thresh = local_moment(thresh_img, mask = mask,
                              nvoxels = nvoxels,
                              moment = 1, center = FALSE,
                              ...)[[1]]
    if (save_imgs){
      write_img(pct_thresh, addstub)
    }
  }

  df$pct_thresh = c(pct_thresh)
  rm(list = c("thresh_img", "pct_thresh")); gc(); gc();


  ################################################
  # Making Probability
  ################################################
  if (verbose) {
    message(paste0("# Getting Top Probability Segmentation from Atropos\n"))
  }
  addstub = paste0("prob_img")
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite){
    prob_img = readnii(fname, reorient = FALSE)
  } else {
    window.masked.img = window_img(masked.img, window = c(0, 100))
    seg = otropos( window.masked.img, i = "KMeans[4]", v = 1)
    prob_img = seg$probabilityimages[[3]] + seg$probabilityimages[[4]]
    rm(list = c("seg")); gc(); gc();
    if (save_imgs){
      write_img(prob_img, addstub)
    }
  }

  df$prob_img = c(prob_img)
  rm(list = c("prob_img")); gc(); gc();


  ################################################
  # Making Percent that are 0
  ################################################
  if (verbose) {
    message(paste0("# Getting Percent 0\n"))
  }
  ###### changed to masked.img
  thresh_0 = niftiarr(masked.img, masked.img == 0)
  addstub = "pct_zero_neighbor"
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite){
    pct_zero_neighbor = readnii(fname, reorient=FALSE)
  } else {
    pct_zero_neighbor = local_moment(thresh_0, mask = NULL,
                                     nvoxels = nvoxels,
                                     moment = 1, center = FALSE,
                                     ...)[[1]]
    pct_zero_neighbor = mask_img(pct_zero_neighbor, mask)
    if (save_imgs){
      write_img(pct_zero_neighbor, addstub)
    }
  }
  df$pct_zero_neighbor = c(pct_zero_neighbor)
  rm(list = c("pct_zero_neighbor", "thresh_0")); gc(); gc();

  if (verbose) {
    message(paste0("# Getting Any 0 Neighbors\n"))
  }
  addstub = "any_zero_neighbor"
  df$any_zero_neighbor = (df$pct_zero_neighbor > 0) *1

  ################################################
  # Making Distance to centroid
  ################################################
  if (verbose) {
    message(paste0("# Getting Distance to centroid\n"))
  }
  addstub = "dist_centroid"
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite){
    dist.img = readnii(fname, reorient=FALSE)
  } else {
    centroid = t(which(mask > 0, arr.ind=TRUE))
    all.ind = expand.grid(lapply(dimg, seq))
    colnames(all.ind) = paste0("dim", seq(length(dimg)))
    all.ind = t(as.matrix(all.ind))
    all.ind = all.ind * vdim
    centroid = centroid * vdim
    dist.img = t(all.ind - rowMeans(centroid))
    rm(list = c("all.ind")); gc(); gc();

    dist.img = sqrt(rowSums(dist.img^2))
    dist.img = niftiarr(img, array(dist.img, dim =dimg))
    dist.img = mask_img(dist.img, mask)
    dist.img = datatyper(dist.img,
                         datatype= convert.datatype()$FLOAT32,
                         bitpix= convert.bitpix()$FLOAT32)
    if (save_imgs){
      write_img(dist.img, addstub)
    }
  }
  df$dist_centroid = c(dist.img)
  rm(list = c("dist.img")); gc(); gc();

  ################################################
  # Making 10mm and 20mm smoother
  ################################################
  make_smooth_img = function(sigma){
    if (verbose) {
      message(paste0("# Getting Smooth ", sigma, "\n"))
    }
    addstub = paste0("smooth", sigma)
    if (save_imgs){
      fname = make_fname(addstub)
    } else {
      fname = tempfile()
      if (file.exists(fname)){
        file.remove(fname)
      }
    }
    if (file.exists(fname) & !overwrite) {
      smooth.img = readnii(fname, reorient = FALSE)
    } else {
      smooth.img = fslsmooth(img.fname, sigma=sigma,
                             mask = mask, retimg = TRUE,
                             outfile = fname)
    }
    return(c(smooth.img))
  }
  #   df$smooth2 = make_smooth_img(sigma=2)
  #   df$smooth5 = make_smooth_img(sigma=5)
  smooths = sapply(sigmas, make_smooth_img)
  colnames(smooths) = paste0("smooth", sigmas)
  df = cbind(df, smooths)
  rm(list = c("smooths")); gc(); gc();




  if (verbose) {
    message(paste0("# Z-score to template \n"))
  }
  addstub = "zscore_template"
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite) {
    zscore = readnii(fname, reorient = FALSE)
  } else {
    zscore = extrantsr::zscore_template(img = orig.masked.img,
                                        template.file = template.file,
                                        mean.img = mean.img,
                                        sd.img = sd.img,
                                        typeofTransform = zscore.typeofTransform,
                                        interpolator = zscore.interpolator,
                                        verbose = verbose)
    if (save_imgs){
      write_img(zscore, addstub)
    }
  }
  df$zscore_template = c(zscore)
  rm(list = c("zscore")); gc(); gc();

  if (verbose) {
    message(paste0("# Flipped Difference \n"))
  }
  addstub = "flipped_value"
  fname = make_fname(addstub)
  if (file.exists(fname) & !overwrite) {
    flipped_value = readnii(fname, reorient=FALSE)
  } else {
    flipper = extrantsr::reg_flip(
      t1 = orig.masked.img,
      mask = orig.mask,
      template.file = template.file,
      typeofTransform = flip.typeofTransform,
      interpolator = flip.interpolator,
      verbose = verbose)
    flipper = flipper$t1
    ##########################
    # Take difference
    ##########################
    flipped_value = orig.masked.img - flipper
    rm(list = c("flipper")); gc(); gc();
    if (save_imgs) {
      write_img(flipped_value, addstub)
    }
  }
  df$flipped_value = c(flipped_value)
  rm(list = c("flipped_value")); gc(); gc();

  if (verbose) {
    message(paste0("# Thresholding small values \n"))
  }


  for (icn in seq(ncol(df))) {
    x = df[, icn]
    if (!(class(x) %in% c("factor", "character"))) {
      x[ !is.finite(x) ] = 0
    }
    df[, icn] = x
  }

  df = as.matrix(df)
  low = abs(df) < low_thresh
  df[ low ] = 0

  df = data.frame(df, stringsAsFactors = FALSE)
  if (!is.null(roi)) {
    df$Y = c(roi)
  } else {
    df$Y = NA
  }
  df$mask = c(mask)


  return(list(df = df, keep.ind = keep.ind, nim = nim,
              miss.roi = miss.roi))
}


