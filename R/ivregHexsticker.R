#' View the Official Hex Sticker for the ivreg Package
#' 
#' Open the official hex sticker for the \pkg{ivreg} package in your browser.
#' 
#' @aliases ivregHexsticker
#' @examples 
#' \dontrun{
#' ivregHexsticker()
#' }
#' @export
ivregHexsticker <- function(){
    utils::browseURL(paste0("file://", system.file("etc", "ivreg-hex.pdf", package="ivreg")))
}
