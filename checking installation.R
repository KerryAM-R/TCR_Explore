# usethis::create_package("../TCR.Explore")

devtools::document()
devtools::check()
devtools::build()


usethis::use_roxygen_md()
require('usethis')
use_cc0_license()

require(TCR.Explore)

remotes::install_github("KerryAM-R/TCR_Explore")

TCR.Explore::runApp_TCR_EXPLORE_V1()
