## ----setup, include=FALSE, echo=FALSE-----------------------------------------
# knitr::knit_hooks$set(optipng = knitr::hook_optipng)
# knitr::opts_chunk$set(optipng = '-o7')

knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(fig.width = 12)
knitr::opts_chunk$set(fig.height = 6)

library(immunarch)
# source("../R/testing.R")
# immdata = load_test_data()
data(immdata)

## ---- eval=F------------------------------------------------------------------
#  # To load the data from a single file (note that you don't need to specify the data format):
#  file_path = paste0(system.file(package="immunarch"), "/extdata/io/Sample1.tsv.gz")
#  immdata <- repLoad(file_path)

## ---- eval=F------------------------------------------------------------------
#  # For instance you have a following structure in your folder:
#  # >_ ls
#  # immunoseq1.txt
#  # immunoseq2.txt
#  # immunoseq3.txt
#  # metadata.txt

## ---- eval=F------------------------------------------------------------------
#  # To load the whole folder with every file in it type:
#  file_path = paste0(system.file(package="immunarch"), "/extdata/io/")
#  immdata <- repLoad(file_path)
#  print(names(immdata))
#  
#  # In order to do that your folder must contain metadata file named
#  # exactly "metadata.txt".
#  
#  # In R, when you load your data:
#  # > immdata <- repLoad("path/to/your/folder/")
#  # > names(immdata)
#  # [1] "data" "meta"
#  
#  # Suppose you do not have "metadata.txt":
#  # > immdata <- repLoad("path/to/your/folder/")
#  # > names(immdata)
#  # [1] "data" "meta"
#  

## -----------------------------------------------------------------------------
as_tibble(data.frame(Sample = c("immunoseq1", "immunoseq2", "immunoseq3"), stringsAsFactors = F))

## ---- eval=F------------------------------------------------------------------
#  # Your list of repertoires in immunarch's format
#  DATA
#  # Metadata data frame
#  META
#  
#  # Create a temporary directory
#  dbdir = tempdir()
#  
#  # Create a DBI connection to MonetDB in the temporary directory.
#  con = DBI::dbConnect(MonetDBLite::MonetDBLite(), embedded = dbdir)
#  
#  # Write each repertoire to MonetDB. Each table has corresponding name from the DATA
#  for (i in 1:length(DATA)) {
#    DBI::dbWriteTable(con, names(DATA)[i], DATA[[i]], overwrite=TRUE)
#  }
#  
#  # Create a source in the temporary directory with MonetDB
#  ms = MonetDBLite::src_monetdblite(dbdir = dbdir)
#  res_db = list()
#  
#  # Load the data from MonetDB to dplyr tables
#  for (i in 1:length(DATA)) {
#    res_db[[names(DATA)[i]]] = dplyr::tbl(ms, names(DATA)[i])
#  }
#  
#  # Your data is ready to use
#  list(data = res_db, meta = META)

## ----basic-data---------------------------------------------------------------
top(immdata$data[[1]])

## ---- eval=FALSE--------------------------------------------------------------
#  coding(immdata$data[[1]])

## ---- eval=FALSE--------------------------------------------------------------
#  noncoding(immdata$data[[1]])

## ---- eval=FALSE--------------------------------------------------------------
#  nrow(inframes(immdata$data[[1]]))

## ---- eval=FALSE--------------------------------------------------------------
#  nrow(outofframes(immdata$data[[1]]))

## -----------------------------------------------------------------------------
filter(immdata$data[[1]], V.name == 'TRBV10-1')

## -----------------------------------------------------------------------------
ds = repSample(immdata$data, "downsample", 100)
sapply(ds, nrow)

## -----------------------------------------------------------------------------
ds = repSample(immdata$data, "sample", .n = 10)
sapply(ds, nrow)

