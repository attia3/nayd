
library(usethis)

usethis::use_description(
  fields = list(
    Package = "NAYD",
    Title   = "NASS-anchored Yield Disaggregation Workflow (NAYD)",
    Version = "0.0.0.9000",
    `Authors@R` = 'person("Ahmed", "Attia",
      email = "ahmed.attia3@outlook.com",
      role  = c("aut", "cre")
    )',
    Description = paste(
      "Tools to disaggregate USDA NASS county-level yields to field- and",
      "cluster-scale production units using CDL, NDVI/ET weights, and",
      "NASS-anchored constraints."
    ),
    License   = "MIT + file LICENSE",
    Encoding  = "UTF-8",
    Roxygen   = "list(markdown = TRUE)",
    RoxygenNote = "7.3.2",
    LazyData  = "true"
  )
)
usethis::use_package("rlang", type = "Imports")
# Make a blank NAMESPACE that roxygen will manage
usethis::use_namespace()
usethis::use_mit_license("Ahmed Attia")

remove.packages("NAYD")
# 2. Make sure the directory is really gone (belt and suspenders)
lib <- .libPaths()[1]
unlink(file.path(lib, "NAYD"), recursive = TRUE, force = TRUE)
usethis::proj_set("C:/NAYD")  # or click the project in RStudio

devtools::document()
devtools::check()
devtools::install()
