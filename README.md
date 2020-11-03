
<!-- README.md is generated from README.Rmd. Please edit that file -->

ichseg <img src="man/figures/logo.png" align="right" height="139" />
====================================================================

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/muschellij2/ichseg.svg?branch=master)](https://travis-ci.com/muschellij2/ichseg)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/muschellij2/ichseg?branch=master&svg=true)](https://ci.appveyor.com/project/muschellij2/ichseg)
<!-- badges: end -->

The goal of `ichseg` is to perform preprocessing on computed tomography
(CT) scans, including skull stripping. Computes predictors of
intracerebral hemorrhage (ICH) and uses these to predict a binary
hemorrhage mask from the data.

Citing
------

To cite `ichseg`, you can run:

    citation("ichseg")

    Muschelli J, Sweeney EM, Ullman NL, Vespa P, Hanley DF, Crainiceanu CM
    (2017). "PItcHPERFeCT: Primary Intracranial Hemorrhage Probability
    Estimation using Random Forests on CT." _NeuroImage: Clinical_, *14*,
    379-390.

    A BibTeX entry for LaTeX users is

      @Article{muschelli2017pitchperfect,
        title = {{PItcHPERFeCT}: Primary Intracranial Hemorrhage Probability Estimation using Random Forests on {CT}},
        author = {John Muschelli and Elizabeth M Sweeney and Natalie L Ullman and Paul Vespa and Daniel F Hanley and Ciprian M Crainiceanu},
        journal = {NeuroImage: Clinical},
        volume = {14},
        pages = {379--390},
        year = {2017},
        publisher = {Elsevier},
      }

Installation
------------

You can install `ichseg` from github with:

    # install.packages("devtools")
    devtools::install_github("muschellij2/ichseg")

Requirements
------------

These functions require a working installation of FSL
(<a href="https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation" class="uri">https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation</a>),
which can be installed via Neurodebian as well:
<a href="http://neuro.debian.net/pkgs/fsl-complete.html" class="uri">http://neuro.debian.net/pkgs/fsl-complete.html</a>.

Prediction
----------

In order to segment ICH from an image, use the `ich_segment` function:

    ichseg::ich_segment(img = "/path/to/ct/scan")
