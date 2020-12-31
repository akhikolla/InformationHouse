#--------------------------------------------------
# BoostMLR news
BoostMLR.news <- function(...) {
  newsfile <- file.path(system.file(package="BoostMLR"), "NEWS")
  file.show(newsfile)
}