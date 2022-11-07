


#' SSD data
#'
#' A dataset containing the L.Posthuma SSDs
#'
#' @format A data frame with 12386 rows and 38 variables:
#' \describe{
#'   \item{test_id}{test id}
#'   \item{chemical_name}{chemial name}
#'   ...
#' }
#' @source "Posthuma et al. (2019): data-raw/data/Copy of etc4373-sup-0002-supmat.csv"
"Data.SSD"







#' Base Waterbase Dataset for HI, MCR analysis.
#'
#' A dataset containing the waterbase data
#'
#' @format A data frame with 616795 rows and 12 variables:
#' \describe{
#'   \item{test_id}{test id}
#'   \item{chemical_name}{chemial name}
#'   ...
#' }
#' @source \url{https://cfpub.epa.gov/ecotox/}
"Data.AggBySiteID.2"




#' Dataset for Chemical classification.
#'
#' A dataset containing the chemical classification data
#'
#' @format A data frame with 334 rows and 18 variables:
#' \describe{
#'   \item{CAS}{test id}
#'   \item{Chem.Group}{Industrial, Metal, PAH, Pest, Pharma}
#'   ...
#' }
#' @source Check with repos owner
"chemclass"





#' Dataset for Driver chemicals.
#'
#' A dataset containing the driver data
#'
#' @format A vector of CAS numbers that is frequency in drivers greater than 5%.
#' \describe{
#'   \item{CAS}{test id}
#' }
#' @name driver
#' @source \url{https://github.com/Bayer-Group/WaterBase}
"driver5"

#' @rdname driver
"driver10"


#' Dataset for Station information
#'
#' A dataset containing the station information
#'
#' @format A vector of CAS numbers that is frequency in drivers greater than 5%.
#' \describe{
#'   \item{LOQ.type}{LOQ type}
#'   \item{lat}{latitude}
#' }
#' @name stations1
#' @source \url{https://github.com/Bayer-Group/WaterBase}
"stations1_Mean"

#' @rdname stations1
"stations1_Max"



#' cas list for N sample>=1 per year and Nchemicals>=1 per year .
#'
#' A dataset containing the chemical classification data
#'
#' @format A data frame with 334 rows and 18 variables:
#' \describe{
#'   \item{CAS}{test id}
#'   \item{Chem.Group}{Industrial, Metal, PAH, Pest, Pharma}
#'   ...
#' }
#' @source \url{https://cfpub.epa.gov/ecotox/}
"allcas"
