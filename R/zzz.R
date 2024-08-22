.onLoad <- function(...){
  quietly <- getOption('quietly')
  options(quietly = T)
  pkg_info <- "MatriCom 1.0"
  packageStartupMessage(pkg_info)
  options(quietly = quietly)
}
