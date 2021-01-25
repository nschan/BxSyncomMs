#' Fix node placement
#'
#' @details This function takes a network that should be facetted in a plot and takes care of node placement.
#' @details This is a very specific function for the dataset here.
#' @details The dataframe needs to contain Treatment, Timepoint and Syncom columns or this breaks.
#' @details The dataframe _has_ to be constructed with ggnetwork::ggnetwork(arrow.gap = FALSE)
#' 
#' @param network multiple networks that were placed via igraph and then turned into a dataframe using ggnetwork::ggnetwork
#' @param Syncoms Syncoms in this network
#' @param Treatments Treatments in this network
#' @param Timepoints Timepoints

#' @author Niklas Schandry

fix_node_placement <- function(network = NULL,
                               Syncoms = list("Tolerant", "Random"),
                               Treatments = list("APO", "BOA"),
                               Timepoints =list("24h","96h")) {

  # Extract edges and nodes  
placed_edges <- network %>% filter(!is.na(Treatment))
placed_nodes <- network %>% filter(is.na(Treatment))

# The below works by comparing every node (x and y mapping) with the ends of the edges.
# This requires
placed_nodes <- bind_rows(
  bind_rows(
  lapply(Syncoms, function(syncom){
    lapply(Treatments, function(treat){
      lapply(Timepoints, function(tp){
  placed_nodes %>% 
  mutate(Syncom = case_when(x %in% c(placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>% 
                                         filter(Syncom == syncom) %$%
                                         x,
                                         placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         xend) ~ syncom,
                              TRUE ~ Syncom),
         Treatment = case_when(x %in% c(placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         x,
                                         placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         xend) ~ treat,
                              TRUE ~ Treatment),
         Timepoint = case_when(x %in% c(placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         x,
                                         placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         xend) ~ tp,
                              TRUE ~ Timepoint)) %>% 
  filter(!is.na(Syncom)) %>% 
  filter(!is.na(Timepoint)) %>% 
  filter(!is.na(Treatment))})
    })
    }) %>% unlist(recursive = FALSE)),
bind_rows(
  lapply(Syncoms, function(syncom){
    lapply(Treatments, function(treat){
      lapply(Timepoints, function(tp){
    placed_nodes %>% 
  mutate(Syncom = case_when(y %in% c(placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>% 
                                         filter(Syncom == syncom) %$%
                                         y,
                                         placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         yend) ~ syncom,
                              TRUE ~ Syncom),
         Treatment = case_when(y %in% c(placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         y,
                                         placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         yend) ~ treat,
                              TRUE ~ Treatment),
         Timepoint = case_when(y %in% c(placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         y,
                                         placed_edges %>% 
                                         filter(Treatment == treat) %>% 
                                         filter(Timepoint == tp) %>%
                                         filter(Syncom == syncom) %$%
                                         yend) ~ tp,
                              TRUE ~ Timepoint)) %>% 
  filter(!is.na(Syncom)) %>% 
  filter(!is.na(Timepoint)) %>% 
  filter(!is.na(Treatment))})
    })
    }) %>%
    unlist(recursive = FALSE))) %>% 
   distinct()


bind_rows(placed_nodes,
  placed_edges) %>% 
  distinct()
}

## No timepoint version

fix_node_placement_no_time <- function(network = NULL,
                               Syncoms = list("Tolerant", "Random"),
                               Treatments = list("APO", "BOA")) {
  
  # Extract edges and nodes  
  placed_edges <- network %>% filter(!is.na(Treatment))
  placed_nodes <- network %>% filter(is.na(Treatment))
  
  # The below works by comparing every node (x and y mapping) with the ends of the edges.
  # This requires
  placed_nodes <- bind_rows(
    bind_rows(
      lapply(Syncoms, function(syncom){
        lapply(Treatments, function(treat){
            placed_nodes %>% 
              mutate(Syncom = case_when(x %in% c(placed_edges %>% 
                                                   filter(Treatment == treat) %>% 
                                                   filter(Syncom == syncom) %$%
                                                   x,
                                                 placed_edges %>% 
                                                   filter(Treatment == treat) %>% 
                                                   filter(Syncom == syncom) %$%
                                                   xend) ~ syncom,
                                        TRUE ~ Syncom),
                     Treatment = case_when(x %in% c(placed_edges %>% 
                                                      filter(Treatment == treat) %>% 
                                                      filter(Syncom == syncom) %$%
                                                      x,
                                                    placed_edges %>% 
                                                      filter(Treatment == treat) %>% 
                                                      filter(Syncom == syncom) %$%
                                                      xend) ~ treat,
                                           TRUE ~ Treatment)) %>% 
              filter(!is.na(Syncom)) %>% 
              filter(!is.na(Treatment))})
        }) %>% 
        unlist(recursive = FALSE)),
    bind_rows(
      lapply(Syncoms, function(syncom){
        lapply(Treatments, function(treat){
            placed_nodes %>% 
              mutate(Syncom = case_when(y %in% c(placed_edges %>% 
                                                   filter(Treatment == treat) %>% 
                                                   filter(Syncom == syncom) %$%
                                                   y,
                                                 placed_edges %>% 
                                                   filter(Treatment == treat) %>% 
                                                   filter(Syncom == syncom) %$%
                                                   yend) ~ syncom,
                                        TRUE ~ Syncom),
                     Treatment = case_when(y %in% c(placed_edges %>% 
                                                      filter(Treatment == treat) %>% 
                                                      filter(Syncom == syncom) %$%
                                                      y,
                                                    placed_edges %>% 
                                                      filter(Treatment == treat) %>% 
                                                      filter(Syncom == syncom) %$%
                                                      yend) ~ treat,
                                           TRUE ~ Treatment)) %>% 
              filter(!is.na(Syncom)) %>% 
              filter(!is.na(Treatment))})
      }) %>%
        unlist(recursive = FALSE))) %>% 
    distinct()
  
  
  bind_rows(placed_nodes,
            placed_edges) %>% 
    distinct()
}



