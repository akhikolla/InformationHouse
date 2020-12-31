if( !file.exists("../windows/libcl-master/include/cl.h") ){
  
  if (getRversion() < "3.3.0") setInternet2()
  download.file(
    "https://github.com/PolMine/libcl/archive/master.zip",
    "libcl.zip", method = "wininet", quiet = TRUE
    )
  dir.create("../windows", showWarnings = FALSE)
  unzip("libcl.zip", exdir = "../windows")
  unlink("libcl.zip")
  
}
