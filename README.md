# Synchrony Causes Major Cycles in Deer Populations and Deer-Vehicle Collisions Across Wisconsin: Introduction to the Repository of All Analyses Supporting the Paper

Thomas L. Anderson, Southeast Missouri State University  
Lawrence W. Sheppard, University of Kansas  
Jonathan A. Walter, University of Virginia  
Robert E. Rolley, Wisconsin Department of Natural Resources  
Daniel C. Reuman, University of Kansas  

## Introduction

This repository can be used to reproduce the complete analyses behind the paper "Synchrony Causes Major Cycles in Deer Populations and Deer-Vehicle Collisions Across Wisconsin" and to recompile the paper itself. Data are also included in the repository. 

## How to compile

Knit the makefile.Rmd using R markdown. If all dependencies are in place (see next section) this should re-compute all analyses from data to paper, resulting in three pdfs: MainText.pdf (the main text of the paper), SuppMat.pdf (the supporting information file for the paper), and makefile.pdf (notes on the compilation process - can be useful for error mitigation in the event of failure). 

The knit may take a several hours or a few days, depending on your computer speed, the value of `numsurrog` in SuppMat.Rmd, and other factors. Subsequent knits, if any, can be faster because packages will be installed (see below) and because intermediate results are cached.

If you try to knit MainText.Rmd or SuppMat.Rmd directly, you may have some success, but cross-document references and other features may fail so this is not recommended.

To compile the documents from the command line, use the following: Rscript -e "library(knitr); knit('makefile.Rmd')".

## Dependencies

### Core dependencies

R, R markdown, R studio, latex and bibtex. 

### Dependencies on the R `checkpoint` package

This codes uses the R `checkpoint` package. This is set up in the master file `makefile.Rmd` in the code chunk called `checkpoint_package`, which contains a line of code specifying a date.

    checkpoint("2018-08-31",checkpointLocation = "./")

The `checkpoint` package then automatically scans through other files looking for other required R packages. It then downloads and installs the newest versions of those packages available on the given date. This helps ensure that re-compiling the document uses _exactly_ the same code that was originally used. This can take some time on first run (you are warned) but it is faster on subsequent runs because the packages are already installed. This also means that R package dependencies should only be the `checkpoint` package, since that package should scan for other packages and install them locally. Quite a few MB disk space are used (200-300).

### Dependencies on `pandoc`

The open source program `pandoc` converts documents from one format to another. Here, the `knitr` package uses it to convert the markdown files into `latex` format so that they can then be turned into PDF files. Installers for multiple operating systems are available here: https://pandoc.org/installing.html.

### Dependencies on `pdflatex`

The makefile makes a system call to `pdflatex`, so software supporting that needs to be installed:

  * On Windows, you can use Miktex (https://miktex.org/howto/install-miktex), 
  * On Linux, install latex (e.g., sudo apt-get install texlive), and
  * On Mac, use the MacTeX installer (http://www.tug.org/mactex/)

### Additional dependencies?

If you find additional dependencies were needed on your system, please let us know: reuman@ku.edu. The compilation process was tested by Reuman on Ubuntu 16.04 using R version 3.4.4 and R studio version 1.1.423; and by Anderson on Windows 7 using R ersion 3.3.3 and R studio version 1.0.136. It has not been tested on Mac. We have endeavored to list all dependencies we can think of above, but we have only compiled on our own machines, so we cannot guarantee that additional dependencies will not also be needed on other machines. This repository is intended to record a workflow, and is not designed or tested for distribution and wide use on multiple machine. It is not guaranteed to work on the first try without and hand-holding on arbitrary computing setups.

## Intermediate files:

Knitting the makefile automatically produces a lot of 'intermediate' files. Files ending in `.tex` are the converted documents from `.Rmd` including all the R code output and the rest (files ending `.log`, `.aux`, `.lof`, `.lot`, `.toc`  and `.out` ) are intermediate files that `pdflatex` uses to keep track of various parts of the document. Some of these can be useful for diagnosing problems, if any. 

## Acknowlegements

D. Lyden, J. Steinglein and the WIDNR for data access; and D. Storm, T. Van Deelen, and J. 
Millspaugh for initial consultation on deer populations in Wisconsin. This material is based upon 
work supported by the National Science Foundation under Grant Numbers 17114195 and 1442595. Any 
opinions, findings, and conclusions or recommendations expressed in this material are those of 
the author(s) and do not necessarily reflect the views of the National Science Foundation. The work was also supported by the James S McDonnell Foundation, the United States Department of Agriculture award 2016-67012-24694, and The Nature Conservancy. 

