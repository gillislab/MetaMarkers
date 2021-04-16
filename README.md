# MetaMarkers: robust marker signatures from single-cell data

To install the package, run the following lines in R:

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("gillislab/MetaMarkers")
```

The installation process is usually fast, but can take up to 10 minutes starting from an empty R installation.

You can access vignettes of the package in [HTML format](./doc) or in an interactive [Rmd format](./vignettes). The demo takes approximately 2 minutes to run (not counting installation time). We propose 3 vignettes:

 - [Computation of meta-analytic markers](http://htmlpreview.github.io/?https://github.com/gillislab/MetaMarkers/doc/MetaMarkers.html)
 - [Cell type annotation](http://htmlpreview.github.io/?https://github.com/gillislab/MetaMarkers/doc/Annotation.html)
 - [Hierarchical cell type annotation](http://htmlpreview.github.io/?https://github.com/gillislab/MetaMarkers/doc/HierarchicalAnnotation.html)

A MetaMarkers webserver is under construction. On the MetaMarkers webserver, you can upload your dataset
and obtain annotations for 85 neuronal cell types defined by the BICCN based on meta-analytic markers.
To access the demo, visit:

 - http://milton.cshl.edu/MetaMarkers/

The package has been tested on recent R versions (3.6 and 4.0) on Linux and Windows and is expected to work on any decently recent R installation.
