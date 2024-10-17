#' create long data.frame from array of posterior samples meeting a certain
#' filter criteria 
#' 
#' @param posterior array of posterior draws
#' @param filter_criteria string used to filter to only desired posterior draws
#'                        (should be of the format `paste(var_name, "\\[")`)
filter_posterior <- function(posterior, filter_criteria){
  if(length(dim(posterior)) == 2) {posterior_filtered <- posterior[,grepl(filter_criteria, dimnames(posterior)[[2]])]}
  else if(length(dim(posterior)) == 3) {posterior_filtered <- posterior[,,grepl(filter_criteria, dimnames(posterior)$parameters)]}
  posterior_long <- as.data.frame(posterior_filtered) %>%
    mutate(draw_id = 1:nrow(as.data.frame(posterior))) %>%
    melt(c("draw_id")) %>%
    mutate(t = as.integer(substr(as.character(variable), nchar(filter_criteria), nchar(as.character(variable))-1)))
  return(posterior_long)
}

