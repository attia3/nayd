install.packages(c("rgee","reticulate","geojsonio","V8","sf","terra"))

install.packages(c("usethis", "devtools", "roxygen2", "testthat"))

#- the correct order here 
Sys.setenv(
  RETICULATE_PYTHON = "C:/Users/ahmed.attia/AppData/Local/r-miniconda/envs/rgee311/python.exe"
)
library(reticulate)
py_config()
library(rgee)
ee_Initialize()
ee_Authenticate()

