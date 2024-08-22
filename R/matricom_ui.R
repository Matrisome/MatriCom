#' matricom.app: launch matricom as an interactive, local app
#'
#' @return a local instance of the matricom shinyapp, without connection to the (online-only) OA data
#' @export
matricom.app <- function(){
  library(shiny)
  library(Matrix)
  library(DT)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(shinyWidgets)
  library(shinybusy)
  library(shinyjs)
  library(shinyhelper)
  library(shinycssloaders)
  library(shinyBS)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratObject)
  library(scCustomize)
  library(sp)
  library(uwot)
  library(shinydashboard)
  library(scater)
  library(googledrive)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(gridExtra)
  library(writexl)
  library(BiocManager)
  library(packcircles)
  library(ggiraph)
  library(plotly)
  library(viridis)
  library(collapsibleTree)
  library(igraph)
  library(scales)
  library(qs)
  library(igraph)
  library(stringr)
  library(reshape2)
  library(htmltools)

  # where is webApp? Find and run from there
  app_dir <- system.file("webApp", package = "MatriCom")

  if (app_dir == "") {
    stop(
      "Could not find the app directory. Try re-installing `MatriCom`.",
      call. = FALSE
    )
  }

  runApp(appDir = app_dir, launch.browser = T)
}

