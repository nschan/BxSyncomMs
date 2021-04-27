#' Get counts
#'
#' @details This function gets OTU level estimated counts from DESeq2
#' The treatment string needs to be part of the sample_name of the phyloseq sample_data slot
#' 
#' @param phyl phyloseq object
#' @param treatvar treatment variable, has to be part of sample data
#' @param OTUcol which taxonomic rank
#' @author Niklas Schandry

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
    # Drop amplicon column
    dplyr::select(-Amplicon) %>%
    # Use genus as rownames
    column_to_rownames(eval(OTUcol))
  counts
}

#' Create correlations
#'
#' @details This function creates a OTU level correlations on estimated counts from DESeq2
#' The treatment string needs to be part of the sample_name of the phyloseq sample_data slot
#' 
#' @param phyl phyloseq object
#' @param treatvar treatment variable, has to be part of sample data
#' @param treatment the treatment of interest
#' @param OTUcol which taxonomic rank
#' @param mincounts minimum counts required to retain taxon, default is 99
#' @param method the correlation coefficient calculated, default is spearman, see correlation::correlation
#' @param p.adj see correlation::correlation, how to adjust p?
#' @param partial see correlation::correlation, compute partial correlations?
#' @param return_counts return count table from DESeq without computing correlations?
#' @author Niklas Schandry

phyloseq_corrs <- function(phyl,
                           treatvar, 
                           treatment,
                           mincounts = 99,
                           OTUcol = "Strain",
                           method = "spearman",
                           p.adj = "holm",
                           partial = FALSE){
  counts <- get_counts(phyl,
                       treatvar,
                       OTUcol)

  ## Move sample names into sample column
  counts %<>%
    t() %>% # t() transposes. Transposition is switching rows are columns (sort of like rotating the table 90°, but not really)
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(Sample = str_replace(Sample, "\\.","_")) # This should have no effect, is a remnant from course development. Also does not break anything
  counts %<>% dplyr::mutate(Treatment = "Other")
  # Format table
  ## Pull treatment from sample name
  for(i in 1:length(treatment)){
  counts %<>% dplyr::mutate(Treatment =
                              case_when(stringr::str_detect(Sample, treatment[i]) ~ paste(treatment[i]), 
                                        TRUE ~ Treatment))
  }
  out_df <- counts %>%
    dplyr::filter(Treatment %in% treatment) %>% 
    dplyr::select(-Sample, -Treatment)
  
  
  out_df <- out_df[, colSums(out_df) >= mincounts] 
  
  out_df <- out_df %>%
    correlation::correlation(method = method, p_adjust = p.adj, partial = partial)

  out_df
}
