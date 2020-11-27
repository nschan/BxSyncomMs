#' Overlay corrnet
#'
#' @details Overlays two (correlation) networks. Typically, these are the result of computing two networks, one for each condition
#' @param network1 first network
#' @param network2 second network
#' @param treat1 treatment for first network, used for formatting
#' @param treat2 treatment for second network, used for formatting
#' @param cutoff correlation coefficient cutoff
#' @param return_corrs return correlation table, default: FALSE
#' @param detailed return detailed comparisons (default: TRUE), or just edge presence absence information, ignoring the correlation sign (TRUE)?
#' @author Niklas Schandry


overlay_corrnet <- function(network1 = NULL, 
                            network2 = NULL,
                            treat1 = "treat1",
                            treat2 = "treat2",
                            cutoff = 0.5 ,
                            return_corrs = FALSE,
                            detailed = TRUE) {
  require(igraph)
  
  full_net <- full_join(network1 %>% filter(abs(r) > cutoff),
                        network2 %>% filter(abs(r) > cutoff),
                        by = c("x", "y"), 
                        suffix =c("_1", "_2")) # This is much easier than tryin to use the treat1 and treat2 vars here, also the vars are not helpful at this stage.
  
  full_net %<>%
    mutate(corr_1 = case_when(r_1 < 0 ~ "Negative",
                              r_1 > 0 ~ "Positive",
                              TRUE ~ "NA"),
           corr_2 = case_when(r_2 < 0 ~ "Negative",
                              r_2 > 0 ~ "Positive",
                              TRUE ~ "NA") )
 if(!detailed){
  full_net %<>%  mutate(
    edge = case_when(
      corr_1 == corr_2 ~ "No change",
      ((corr_1 == "Negative" & corr_2 == "Positive") | (corr_1 == "Positive" & corr_2 == "Negative")) ~ "Sign change",
      (corr_2 == "NA") & (corr_1 != "NA") ~ paste(treat1, "Only"),
      (corr_2 != "NA") & (corr_1 == "NA") ~ paste(treat2, "Only"),
      TRUE ~ "check me")) %>% 
    distinct()
 } else {
   full_net %<>%  mutate(
     edge = case_when(
       (corr_1 == "Positive") & (corr_2 == "Positive") ~ "Positive (both)",
       (corr_1 == "Negative") & (corr_2 == "Negative") ~ "Negative (both)",
       (corr_1 == "Negative") & (corr_2 == "Positive")  ~ paste("Negative in", treat1, 
                                                                "Positive in", treat2),
       (corr_1 == "Positive") & (corr_2 == "Negative") ~ paste("Positive in",treat1,                                                                        "Negative in", treat2),
       (corr_2 == "NA") & (corr_1 == "Positive") ~ paste(treat1, "Positive"),
       (corr_2 == "NA") & (corr_1 == "Negative") ~ paste(treat1, "Negative"),
       (corr_1 == "NA") & (corr_2 == "Positive") ~ paste(treat2, "Positive"),
       (corr_1 == "NA") & (corr_2 == "Negative") ~ paste(treat2, "Negative"),
       TRUE ~ "check me")) %>% 
     distinct()
 }
  
  
  
  if(return_corrs) {
    return(full_net)
  }
  
  full_net <- full_net %>%
    dplyr::select(x,y, edge) %>% 
    igraph::graph_from_data_frame() %>% 
    ggnetwork::ggnetwork()
  return(full_net)
}