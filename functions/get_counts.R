get_counts <- function(phyl,
                       treatvar, 
                       OTUcol = "Strain"){
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
    # Drop amplicon colum
    dplyr::select(-Amplicon) %>%
    # Use genus as rownames
    column_to_rownames(eval(OTUcol))
 counts
 }