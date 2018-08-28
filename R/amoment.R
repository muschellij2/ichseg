
#' @title Calculate Moments of Neighborhood
#'
#' @description This function calculates Local Moments (mean, standard deviation, skew) for an array.
#' @param image input image
#' @param window window (in width) for the neighborhood
#' @param nvoxels window (in voxels) for the neighborhood 1 results in a 3x3 cube
#' @param moment vector moments taken (1- mean, 2-sd, 3-skew)
#' @param mask array or object of class nifti of same size as image
#' @param only.mask Should objects outside the mask (i.e. zeros) be
#' counted the moment?  Default is FALSE so edges are weighted to 0
#' @param center vector of indicator of central moment.
#' if TRUE mean image is subtracted.  Same length as moment
#' @param invert Standardize the values by the power: 1/moment
#' @param mean_image mean image to be subtracted. If not supplied, and central = TRUE, local_moment_edge is run with mom = 1
#' @param na.rm remove NAs from the moment image calculation
#' @param remask set areas outside of mask to 0
#' @param ... Arguments passed to \code{\link{get_neighbors}}
#' @importFrom magic ashift
#' @importFrom neurobase niftiarr datatyper
#' @export
#' @return List of arrays same lenght as moment
#' @examples
#' x = array(rnorm(1000), dim=c(10, 10, 10))
#' mask = abs(x) < 1
#' mean.x = local_moment(x, nvoxels=1, moment = 1, mask=mask,
#' center = FALSE,
#' remask = FALSE)[[1]]
#' var.x = local_moment(x, nvoxels=1, moment = 2, mask=mask, center = TRUE,
#' mean_image = mean.x, remask=FALSE)[[1]]
#'
#' ### check that x[2,2,2] mean is correct
#' check = x[1:3,1:3,1:3]
#' ## masking
#' vals = check[abs(check) < 1]
#' m = mean(vals)
#' all.equal(m, mean.x[2,2,2])
#' n = length(vals)
#' v = var(vals) * (n-1)/n
#' var.x[2,2,2]
#' all.equal(v, var.x[2,2,2])
#'
local_moment <- function(
  image,
  window = NULL,
  nvoxels=NULL,
  moment,
  mask = NULL,
  only.mask = FALSE,
  center=is.null(mean_image),
  invert = FALSE,
  # the mean
  mean_image=NULL, # mean image to be subtracted.  If not supplied
  # local_moment_edge is run with mom = 1
  na.rm=TRUE, # remove NAs from the image (mask set to 0),
  remask = TRUE,  # set areas outside of mask to 0
  ...
  ) {

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }

  ## if mask is null - whole image
  dimg = dim(image)[1:3]
  if (is.null(mask)) {
    mask = array(1, dim=dimg)
  }

  lmom = length(moment)
  lcent = length(center)
  stopifnot(lmom == lcent)

  neigh = get_neighbors(image, window = window, nvoxels=nvoxels,
                mask = mask, ...)

  neighbors = neigh$neighbors
  indices = neigh$indices

  mn = rowMeans(neighbors, na.rm=na.rm)
  img_list = vector("list", length = lmom)
  for (imom in seq(lmom)){
    cent = center[imom]
    mom = moment[imom]
    moment_image = neighbors
    if (isTRUE(cent)){
      if (mom != 1) moment_image = moment_image - mn
    }

    moment_image = moment_image^mom
    moment_image = rowMeans(moment_image, na.rm=na.rm)

    realpow <- function(x, pow) {
      sgn = sign(x)
      x = abs(x)
      x = x^{pow}
      x = sgn * x
    }
    ### put on same scale
    if (invert) {
      moment_image = realpow(moment_image, pow = 1/mom)
    }

    moment_image = array(moment_image, dim=dimg)
    if (remask) {
      moment_image = moment_image * mask
    }
    if ( inherits(image, "nifti") ){
      moment_image = niftiarr(image, moment_image)
      moment_image = datatyper(moment_image,
                              datatype= convert.datatype()$FLOAT32,
                              bitpix= convert.bitpix()$FLOAT32)
    }
    img_list[[imom]] = moment_image
  }

  # moment <- array(moment, dim = dim(image))
  return(img_list)
}


#' @title Calculate Moments of Neighborhood
#'
#' @description This function calculates Local Moments (mean, standard deviation, skew) for an array.
#' @param image input image
#' @param window window (in width) for the neighborhood
#' @param nvoxels window (in voxels) for the neighborhood 1 results in a 3x3 cube
#' @param mask array or object of class nifti of same size as image
#' @param rm.value remove the voxel itself before taking the moment
#' @param check.wrap Logical - check wrapround and put \code{rep.value} in
#' for the wrapped values
#' @param rep.value Replace wrapped values (edge of image) to this value
#' @param ... Not used
#' @importFrom magic ashift
#' @export
#' @return Array with same dimensions as image
#' @examples
#' x = array(rnorm(1000), dim=c(10, 10, 10))
#' neigh = get_neighbors(x, nvoxels = 1)
get_neighbors <- function(
  image,
  window = NULL,
  nvoxels=NULL,
  mask = NULL,
  rm.value = FALSE, # remove the voxel itself before returning
  check.wrap = TRUE,
  rep.value = 0,
  ...
) {

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }

  ## if mask is null - whole image
  dimg = dim(image)[1:3]
  if (is.null(mask)) {
    mask = array(1, dim=dimg)
  }

  mypermutations = function(winds){
    windlist = lapply(1:3, function(x) winds)
    eg = expand.grid(windlist)
    eg = eg[ order(eg$Var1, eg$Var2, eg$Var3), ]
    eg = as.matrix(eg)
    rownames(eg) = NULL
    colnames(eg) = NULL
    eg
  }

  ### initalize array
  ### different ways of parameterizing the "window"
  if (is.null(nvoxels)) {
    stopifnot(!is.null(window))
    stopifnot(length(window) == 1)
    if ((window %% 2) != 1) {
      stop("window must be odd number")
    }
    winds = (-window/2 + .5):(window/2 - .5)
#     indices <- gtools::permutations(window, 3,
#                             v= (-window/2 + .5):(window/2 - .5),
#                             repeats.allowed=TRUE)
  }  else {
    stopifnot(is.wholenumber(nvoxels))
    stopifnot(is.null(window))

    winds <- (-nvoxels):nvoxels
#     indices <- gtools::permutations(length(winds), 3, v = winds,
#                             repeats.allowed=TRUE)
  }
  indices = mypermutations(winds)

  if (rm.value){
    allzero = apply(indices == 0, 1, all)
    indices = indices[!allzero,]
  }

  image = image * mask

  nruns <- nrow(indices)
  mat = matrix(NA, nrow=prod(dimg), ncol=nruns)
  for ( i in 1:nruns){
    shifter <- ashift(image, indices[i,])
    mat[, i] = shifter
  }

  if (check.wrap){
    tall.ind = t(as.matrix(expand.grid(dim1=seq(dimg[1]),
                                      dim2=seq(dimg[2]),
                                      dim3=seq(dimg[3]))))
    dimg.mat = matrix(dimg, ncol=prod(dimg), nrow = 3)
    for ( i in 1:nruns){
      all.ind2 = tall.ind + indices[i,]
      outside = {all.ind2 < 1} + all.ind2 > dimg.mat
      any.outside = colSums(outside) > 0
      mat[any.outside,i] = rep.value
    }
  }

  return(list(neighbors=mat, indices = indices))
}


#' @title Calculate Mean Image by FFT
#'
#' @description This function calculates the mean image using fft
#' @param x 3D array
#' @param nvoxels window (in voxels) for the neighborhood
#' (1 results in a 3x3 cube)
#' @param shift Should results be shifted back?
#' @param verbose Diagnostic outputing
#' @importFrom magic ashift
#' @importFrom stats fft
#' @export
#' @return Array of size of x
mean_image = function(x, nvoxels, shift = TRUE, verbose = TRUE){

  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }
  stopifnot(length(nvoxels) %in% c(1, 3))

  if (length(nvoxels) == 1){
    stopifnot(is.wholenumber(nvoxels))
    n <- length((-nvoxels):nvoxels)

    h = array(1, dim=c(n, n, n))
    h = h / length(h)
  }
  if (length(nvoxels) == 3){
    dims = (nvoxels*2) + 1
    stopifnot(is.wholenumber(dims))
    h = array(1, dim=dims)
    h = h / length(h)
  }
  if (verbose){
    cat(paste0("Dimension of kernel is ", paste(dim(h), collapse = "x"), "\n"))
  }

  # % x - 3dim matrix
  # % h - smoothing kernel - 3dim matrix size(h) =< size(x)

  x[is.na(x)]=0;

  dx = dim(x)
  RX = dx[1];
  CX = dx[2];
  SX = dx[3];

  dh = dim(h)
  RH = dh[1];
  CH = dh[2];
  SH = dh[3];

  z = array(0, dim = dx);
  z[1:RH,1:CH,1:SH] = h;

  if (verbose){
    cat("Creating Kernel fft\n")
  }
  H = fft(z);
  rm(z)
  if (verbose){
    cat("Creating Image fft\n")
  }
  X = fft(x);
  rm(x)
  if (verbose){
    cat("Convolution\n")
  }
  X = H * X
  rm(H)
  y = fft(X, inverse = TRUE)
  y = y / length(y)
  y = Re(y)

#   print(ceiling(RH/2))
#   print(ceiling(CH/2))
#   print(ceiling(SH/2))
  #   [p+1:m 1:p]
  if (verbose){
    cat("Shifting\n")
  }
  if (shift) {
    y = ashift(y, v = -(ceiling(dim(h)/2)-1))
  }

  #   y=circshift(y,c(-ceiling(RH/2), -ceiling(CH/2), -ceiling(SH/2)))
  # compensation for group delay
  return(y)
}

