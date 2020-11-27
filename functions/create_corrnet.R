#' Create correlation netowrk
#'
#' @details This function creates a OTU level correlation networks.
#' The treatment string needs to be part of the sample_name of the phyloseq sample_data slot
#' 
#' @param phyl phyloseq object
#' @param treatvar treatment variable, has to be part of sample data
#' @param treatment the treatment of interest
#' @param OTUcol which taxonomic rank
#' @param mincounts minimum counts required to retain taxon, default is 99
#' @param method the correlation coefficient calculated, default is spearman, see corrr::correlate
#' @param use see corrr::correlate, how to deal with missing data?
#' @author Niklas Schandry

create_corrnet <- function(phyl,
                           treatvar, 
                           treatment,
                           mincounts = 99,
                           OTUcol = "Strain",
                           method = "spearman",
                           use = "pairwise.complete.obs"){
  counts <- phyl %>% 
    phyloseq::phyloseq_to_deseq2(reformulate(treatvar)) %>% 
    DESeq2::DESeq() %>% 
    DESeq2::counts() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("Amplicon")
  
  ## Add genus information
  counts %<>%  dplyr::left_join(phyl %>%    
                                  phyloseq::tax_table() %>%
                                  as.data.frame() %>%
                                  tibble::rownames_to_column("Amplicon") %>%
                                  dplyr::select(Amplicon, eval(OTUcol)), by = "Amplicon")  %>%
    # Drop amplicon column
    dplyr::select(-Amplicon) %>%
    # Use genus as rownames
    column_to_rownames(eval(OTUcol))
  
  ## Move sample names into sample column
  counts %<>% 
    t() %>% # t() transposes. Transposition is switching rows are columns (sort of like rotating the table 90°, but not really)
    as.data.frame() %>% 
    tibble::rownames_to_column("Sample") %>% 
    dplyr::mutate(Sample = str_replace(Sample, "\\.","_")) # This should have no effect, is a remnant from course development. Also does not break anything
  
  # Format network 1
  ## Pull treatment from sample name
  counts %<>% dplyr::mutate(Treatment =
                              case_when(stringr::str_detect(Sample, treatment) ~ paste(treatment), 
                                        # Changed this so the treatment is detected from the sample name.
                                        # Maybe sub-optimal in general, but fine for this course..
                                        TRUE ~ "Other"))
  
  # Format network 1
  
  network <- counts %>%
    dplyr::filter(Treatment == treatment) %>% 
    dplyr::select(-Sample, -Treatment)
  
  
  network <- network[, colSums(network) > mincounts] 
  
  network <- network %>% 
    corrr::correlate(method = method, use = use) %>% 
    corrr::stretch(remove.dups = T) %>% 
    na.omit()
  
  network
}