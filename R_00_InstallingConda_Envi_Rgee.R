install.packages(c("rgee","reticulate","geojsonio","V8","sf","terra"))

install.packages(c("usethis", "devtools", "roxygen2", "testthat"))

#- the correct order here 
# Optional: manually point reticulate to the Python env used by rgee
Sys.setenv(
  RETICULATE_PYTHON = "~/AppData/Local/r-miniconda/envs/rgee_py311/python.exe"
)
library(reticulate)
py_config()

library(rgee)
ee_Initialize()
ee_Authenticate()

