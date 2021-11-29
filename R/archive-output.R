upload_outputs_archive <- function(...) {
  zip("outputs.zip", "output")
  piggyback::pb_upload("outputs.zip")
  fs::file_delete("outputs.zip")
}

get_outputs_archive <- function(dir = ".", ...) {
  piggyback::pb_download("outputs.zip")
  unzip("outputs.zip", exdir = dir)
  fs::file_delete("outputs.zip")
}
