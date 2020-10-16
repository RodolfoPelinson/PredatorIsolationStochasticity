#' Abundances for the second and third surveys
#'
#'
#' @format Numerical vector w/ with 44 observations
#' @source {Experimental Data}
"com_SS2_SS3_abundance"


#' Pond IDs for the second and third surveys
#'
#' @format Factor w/ 24 levels and 44 observations
#' @source {Experimental Data}
"ID_SS2_SS3"

#' Sampling Survey identities for the second and third surveys
#'
#' @format Factor w/ 2 levels: "1","2", and 44 observations
#' @source {Experimental Data}
"SS_SS2_SS3"

#' Isolation identities for the second and third surveys
#'
#' @format Factor w/ 3 levels: "30","120", "480" and 44 observations
#' @source {Experimental Data}
"isolation_SS2_SS3"

#' Fish treatment identities for the second and third surveys
#'
#' @format Factor w/ 2 levels: "absent","present", and 44 observations
#' @source {Experimental Data}
"fish_SS2_SS3"

#' Fish, isolation and surveyidentities combined for the second and third surveys
#'
#' @format Factor w/ 12 levels and 44 observations
#' @source {Experimental Data}
"All"


#' Traits and other relevant information considering the second and third surveys only
#'
#' Traits and other relevant information about each taxa sampled
#'
#' @format A data frame with 36 rows or 'observations' and 16 variables:
#' \describe{
#'   \item{species}{Factor variable: Ids of the sampled taxa}
#'   \item{trophic}{Factor variable: Trophic level of the sampled taxa}
#'   \item{microhabitat}{Factor variable: Microhabitat level of the sampled taxa}
#'   \item{volume}{Numeric variable: Volume of the largest individual of each sampled taxa}
#'   \item{mass}{Numeric variable: Wet mass of the largest individual of each sampled taxa}
#'   \item{volume_log}{Numeric variable: Log transformed volume of the largest individual of each sampled taxa}
#'   \item{family}{Character variable: Family of each sampled taxa}
#'   \item{order}{Character variable:  Order of each sampled taxa}
#'   \item{total_ab}{Numeric variable: Total abundance of each sampled taxa considering all surveys}
#'   \item{total_ab_log}{Numeric variable: Log transformed total abundance of each sampled taxa considering all surveys}
#'   \item{total_ab_SS1}{Numeric variable: Total abundance of each sampled taxa considering the first survey}
#'   \item{total_ab_SS2}{Numeric variable: Total abundance of each sampled taxa considering the second survey}
#'   \item{total_ab_SS3}{Numeric variable: Total abundance of each sampled taxa considering the third survey}
#'   \item{total_ab_SS1_log}{Numeric variable: Log transformed total abundance of each sampled taxa considering the first survey}
#'   \item{total_ab_SS2_log}{Numeric variable: Log transformed total abundance of each sampled taxa considering the second survey}
#'   \item{total_ab_SS3_log}{Numeric variable: Log transformed total abundance of each sampled taxa considering the third survey}
#' }
#' @source {Experimental Data}
"Trait_SS2_SS3"
