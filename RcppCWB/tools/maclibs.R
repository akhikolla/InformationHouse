instdir <- commandArgs()[6]
target_tarball <- file.path(instdir, "libglib.tar.gz")
target_dir <- file.path(instdir, "macos")
  
if (getRversion() < "3.3.0") setInternet2()
  
download.file(
  "https://github.com/PolMine/libglib/archive/master.tar.gz",
  target_tarball, method = "libcurl", quiet = FALSE
)
  
dir.create(target_dir, showWarnings = FALSE)
  
untar(target_tarball, exdir = target_dir)
unlink(target_tarball)
  
pkgconfig_file <- file.path(target_dir, "libglib-master/pkgconfig/glib-2.0.pc")
pc <- readLines(pkgconfig_file)
pc[1] <- sprintf("prefix=%s/macos/libglib-master", getwd())
writeLines(text = pc, con = pkgconfig_file)

