.onLoad <- function(...){
  quietly <- getOption('quietly')
  options(quietly = T)
  pkg_info <- "matRicom v1.0.1"
  packageStartupMessage(pkg_info)
  options(quietly = quietly)
}
