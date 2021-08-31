library("devtools")
library(roxygen2)

setwd("..")
db.15 <- readRDS("Cdiff_DB_lite.15.RDS")
setwd("Desktop/Cdiff_fragmentR")
getwd()


create("CdiffFragR")

file.copy("Functions.8.30.R",
          "CdiffFragR/R/Functions.8.30.R", overwrite = T)

setwd("CdiffFragR")
document()

use_data(db.15, "Cdiff_DB_lite.15.rds")






