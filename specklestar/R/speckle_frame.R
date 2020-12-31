#' Get selected speckle frame
#' Get specified speckle frame as matrix from file
#'
#' @param data_file a character string with the path name to a file.
#' @param frame an integer.
#' @return 512 x 512 matrix with given frame.
#' @examples
#' \dontrun{
#' # On Unix-like operating systems only
#' # Read frame number 2 from file to matrix
#' obj_filename <- system.file("extdata", "ads15182_550_2_frames.dat", package = "specklestar")
#' frame2 <- speckle_frame(obj_filename, 2)
#' }
#' @export
speckle_frame <- function(data_file = file.choose(), frame = 1) {
  tmp_file <- paste(tempdir(), '/tmp.dat', sep = '')
  system(sprintf("dd if=%s of=%s bs=512*512*2 skip=%d count=1", data_file, tmp_file, frame - 1),
         ignore.stderr = TRUE)

  file_connector = file(tmp_file, "rb")
  frame <- readBin(file_connector, integer(), endian = "little", n = 512 * 512, size = 2)
  frame <- matrix(frame, 512, 512)

  close(file_connector)
  file.remove(tmp_file)

  return(frame)
}
