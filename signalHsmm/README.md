[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/signalHsmm)](https://cran.r-project.org/package=signalHsmm)
[![Downloads](http://cranlogs.r-pkg.org/badges/signalHsmm)](https://cran.r-project.org/package=signalHsmm)
[![Build Status](https://api.travis-ci.org/michbur/signalHsmm.png)](https://travis-ci.org/michbur/signalHsmm)

<img src="https://github.com/michbur/signalHsmm/blob/master/inst/signal_gui/logo.png" alt="signalHsmm" style="height: 200px;"/>

Predict Presence of Signal Peptides
-------------------------

signalHsmm predicts presence of signal peptides in eukaryotic proteins using hidden semi-Markov models. The implemented algorithm can be accessed as a web-based service http://smorfland.uni.wroc.pl/shiny/signalHsmm/.

Local instance of signalHsmm
------------------------
signalHsmm can be also used locally as the R package. It can be installed from CRAN using:

```R
install.packages("signalHsmm")
```

You can install the latest development version of the package directly from github:

```R
source("https://install-github.me/michbur/signalHsmm")
```

After the installation, the GUI can be accessed locally:

```R
library(signalHsmm)
gui_signalHsmm()
```
All signalHsmm functionalities can be also invoked in the batch mode, for example:

```R
run_signalHsmm(benchmark_dat[1:10])
```
