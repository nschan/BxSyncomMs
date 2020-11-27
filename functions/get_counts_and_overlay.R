#' Get counts and overlay
#'
#' @details This function creates two OTU level correlation networks, one for each treatment. 
#' Unless return_corrs is true, these tables are passed to overlay_corrnet. 
#' The treat1 and treat2 strings need to be part of the sample_name of the phyloseq sample_data slot
#' 
#' @param phyl phyloseq object
#' @param treatvar treatment variable, has to be part of sample data
#' @param treat1 first treatment
#' @param treat2 second treatemnt
#' @param OTUcol which taxonomic rank
#' @param mincounts minimum counts required to retain taxon
#' @param mincorr minimum correlation coefficient to retain correlation
#' @param method the correlation coefficient calculated, default is spearman, see corrr:correlate
#' @param use see corrr::correlate, how to deal with missing data?
#' @param return_corrs return correlation table, default: FALSE (passed to overlay_corrent)
#' @author Niklas Schandry

get_counts_and_overlay <- function(phyl, # Phyloseq object
                                   treatvar = "sample_treatment", # Column that contains the treatments
                                   treat1, # First treatment
                                   treat2, # Second treatment
                                   OTUcol = "Genus", # Column that defines OTU
                                   mincounts = 99, # Minimum total counts per treatment to retain OTU for correlations
                                   mincorr = 0.5, # Correlation cutoff
                                   method = "spearman",
                                   use = "pairwise.complete.obs",
                                   return_corrs = FALSE
) {
  # Extract counts
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
  
  ## Pull treatment from sample name
  counts %<>% dplyr::mutate(Treatment =
                              case_when(stringr::str_detect(Sample, treat1) ~ paste(treat1), 
                                        # Changed this so the treatment is detected from the sample name.
                                        # Maybe sub-optimal in general, but fine for this course..
                                        stringr::str_detect(Sample, treat2) ~ paste(treat2),
                                        TRUE ~ "Other"))
  
  # Format network 1
  
  network1 <- counts %>%
    dplyr::filter(Treatment == treat1) %>% 
    dplyr::select(-Sample, -Treatment)
  
  network1 <- network1[, colSums(network1) > mincounts] 
  
  network1 <- network1 %>% 
    corrr::correlate(method = method, use = use) %>% 
    corrr::stretch(remove.dups = T) %>% 
    na.omit()
  
  # Format network2
  
  network2 <- counts %>%
    dplyr::filter(Treatment == treat2) %>% 
    dplyr::select(-Sample, -Treatment)
  
  network2 <- network2[, colSums(network2) > mincounts] 
  
  network2 <- network2 %>% 
    corrr::correlate(method = method, use = use) %>% 
    corrr::stretch(remove.dups = T) %>% 
    na.omit()
  
  # Pass to overlay function
  overlay_corrnet(network1,
                  network2,
                  treat1 = treat1,
                  treat2 = treat2,
                  cutoff = mincorr,
                  return_corrs = return_corrs)
  
  
}
