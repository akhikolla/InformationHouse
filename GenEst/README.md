# GenEst
<img src = 'inst/app/www/GenEst.png' height = '80' align="right" />

## GenEst: Generalized Fatality Estimator    

**GenEst** is a tool for estimating mortalities from efficiency, persistence,
and carcass data.

## DISCLAIMER

This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.

## Installation
Setup and installation require several steps. Do not skip any steps.

### Updated version of R (>= 3.5.0, released on 23 April 2018):
R is free and open source software for statistical computing. If R is not installed on your computer or if your version of R is <3.5.0, download and install the latest version from https://cran.r-project.org/, following the instructions provided at the site. In particular, "Download" and then "install R for the first time" (if working in Windows), or "Download" and then follow the further instructions on the subsequent web page (if working on Mac OS or Linux-like OS). If you already have an older copy of R installed on your computer, the new version will be installed alongside the old. Unless you know a reason why you want to keep both versions, it is usually a good idea to uninstall the old version to avoid confusion and clutter. 

NOTE TO EXPERIENCED R USERS: When you install a new version of R, packages that you previously installed under an older version may not be immediately available to the new R. If not, the easiest way to make them available is to copy the package folders in your old "library" folder into the "library" folder in your new R installation. Then, enter `update.packages()` in R. If asked about a CRAN mirror, choose the nearest location. If you are working in Windows OS and are asked whether you want to install packages "from source", choose "No".


### Third-party packages: 
Several third-party pacakges are required; all are free and open source and available from CRAN. The easiest way to install them is to run the following commands in R (with guidance concerning potential dialog boxes given below the commands):

```

package_req <- c("corpus", "DT", "gsl", "gtools", "htmltools", "htmlwidgets", "lubridate", "MASS", "matrixStats", "mvtnorm", "Rcpp", "shiny", "shinyjs",  "survival")
package_new <- package_req[!(package_req %in% installed.packages()[,"Package"])] 
if(length(package_new) > 0) install.packages(package_new)
if (packageVersion("htmlwidgets") < "1.5") install.packages("htmltools")
if (packageVersion("shiny") < "1.4.0") install.packages("shiny")
```
-- If asked about a "CRAN mirror", choose the nearest location.

-- If asked whether you want to use a "personal library", choose "Yes"

-- If you are on Windows and are asked whether you want to install packages and their dependencies "from source", choose "No" (unless you are ready to go to lunch, in which case, you can select "Yes" and the installation may well be done by the time you get back).

### GenEst: 
Click on "Tags" under the "Repository" tab on the left sidebar at https://code.usgs.gov/ecosystems/GenEst and then click the link for the specific release you want. 

-- For Windows, download the compressed folder GenEst_1.x.x.zip (do not unzip) and note where it is stored. You will install from the local .zip folder. 

-- For Mac OS or Unix-like OS, download the compressed file GenEst_1.x.x.tar.gz and note where it is stored. You will install from the local .tar.gz file. 

If you are working directly in R (not R Studio), run the following command:
```
install.packages(file.choose()) # and navigate to the package archive file you just downloaded: GenEst_1.x.x.xxx
```

If you are working in R Studio:

Click "Install" in the Packages pane.

Select "Package Archive File (.zip; .tar.gz)" as "Install from:" in the dialog box.

Browse to where you saved the zip file, and open it so it appears in the "Package archive" space.

Click the Install button on the dialog box.

## Getting Started
### Graphical user interface (GUI): easy-to-use buttons and menus

To start the GUI, open R and enter the command:
```
library(GenEst)
runGenEst()
```

Download the User Guide from a link near the bottom of the "Help" page in the app or from https://pubs.usgs.gov/tm/7c19/tm7c19.pdf

### R command line: more functionality and flexibility
```
library(GenEst)
browseVignettes("GenEst")
?GenEst
```
Also, help files for GenEst functions are accessible in the standard R way, for example:
```
?pkm
```
## Further Reading
GenEst User Guide: https://doi.org/10.3133/tm7C19

GenEst Statistical Models:  https://doi.org/10.3133/tm7A2
