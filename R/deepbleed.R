#' DeepBleed Model
#'
#' @param outdir Output directory for `DeepBleed` model
#'
#' @note \url{https://github.com/muschellij2/deepbleed}
#'
#' @return A list of the output images and predictions.
#' @export
#' @rdname deepbleed
#'
#' @examples
#' \donttest{
#' destfile = file.path(tempdir(), "01.tar.xz")
#' dl = download.file(
#'   "https://archive.data.jhu.edu/api/access/datafile/1311?gbrecs=true",
#'   destfile = destfile)
#' res = untar(tarfile = destfile, exdir = tempdir())
#' fname = file.path(tempdir(), "01", "BRAIN_1_Anonymized.nii.gz")
#' tdir = tempfile()
#' dir.create(tdir)
#' download_deepbleed_model(outdir = tdir)
#' mod = load_deepbleed_model(outdir = tdir)
#' if (fslr::have_fsl()) {
#' predict_deepbleed(fname, outdir = tdir)
#' }
#' }
download_deepbleed_model = function(outdir = NULL) {
  if (is.null(outdir)) {
    outdir = system.file(package = "ichseg")
  }
  fnames = c("_index", "_data-00000-of-00002",
             "checkpoint", "_data-00001-of-00002")
  real_fnames = sub("_", ".", fnames)

  outfiles = file.path(outdir, real_fnames)
  if (!all(file.exists(outfiles))) {
    url = paste0("https://www.dropbox.com/s/v2ptd9mfpo13gcb/",
                 "mistie_2-20200122T175000Z-001.zip?dl=1")
    tfile = tempfile(fileext = ".zip")
    dl = utils::download.file(url, destfile = tfile)

    ofiles_list = utils::unzip(
      tfile,
      exdir = outdir,
      list = TRUE,
      junkpaths = TRUE)
    ofiles = utils::unzip(tfile, exdir = outdir, junkpaths = TRUE)
    stopifnot(all(basename(ofiles) == fnames))
    file.rename(ofiles, outfiles)
  }
  stopifnot(all(file.exists(outfiles)))

  outdir = path.expand(outdir)
  outdir = normalizePath(outdir)
  if (!grepl("/$", outdir)) {
    outdir = paste0(outdir, "/")
  }
  return(outdir)
}

#' @rdname deepbleed
#' @export
load_deepbleed_model = function(outdir = NULL) {
  outdir = download_deepbleed_model(outdir)
  path = system.file("deepbleed", package = "ichseg")
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("You need the reticulate package for deepbleed")
  }
  model = reticulate::import_from_path("models.vnet", path)
  vnet = model$VNet()
  vnet$load_weights(outdir)
  vnet
}

#' @rdname deepbleed
#' @param image image to segment using `DeepBleed` model
#' @param mask brain mask image
#' @param ... additional arguments to send to
#' \code{\link{CT_Skull_Stripper_mask}}
#' @export
predict_deepbleed = function(image,
                             mask = NULL,
                             ...,
                             outdir = NULL) {

  vnet = load_deepbleed_model(outdir = outdir)
  image = check_nifti(image)
  if (is.null(mask)) {
    mask = CT_Skull_Stripper_mask(image, ...)
    mask = mask$mask
  }
  mask = check_nifti(mask)
  ss = mask_img(image, mask)
  template.file = system.file(
    'scct_unsmooth_SS_0.01_128x128x128.nii.gz',
    package = 'ichseg')
  reg = extrantsr::registration(
    ss,
    template.file = template.file,
    typeofTransform = "Rigid",
    affSampling = 64)
  image = reg$outfile
  image = array(image, dim = c(1L, dim(image), 1L))
  prediction = vnet$predict(image)

  arr = drop(prediction)
  arr = neurobase::copyNIfTIHeader(arr =  arr, reg$outfile)
  native = extrantsr::ants_apply_transforms(
    fixed = ss,
    moving = arr,
    interpolator = "nearestNeighbor",
    transformlist = reg$invtransforms,
    whichtoinvert = 1)
  L = list(
    brain_mask = mask,
    skull_stripped = ss,
    template_space = image,
    template_prediction = arr,
    native_prediction = native)
  return(L)

}