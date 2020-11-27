#' Compare conditions
#' @details This function takes a phyloseq object and performs a OTU-level comparison between two conditions returns log2FC, does filtering on adjusted pval, controlled via alpha
#' @param phyl a phyloseq object
#' @param treatvar character, the treatment variable. Has to be part of the phyloseq object sample_data
#' @param treat1 first treatment; will be numerator
#' @param treat2 second treatment; will be denominator
#' @param OTUcol which OTU assignment level; default "Genus"
#' @param return_unfiltered logical, return without p-value filtering
#' @param alpha significance level, ignored if return_unfiltered is TRUE
#' @param plot logical, return plot? Default: TRUE; when FALSE returns table.
#' @param fit_type fitType used by the DESeq function.
#' @author Niklas Schandry 

compare_conditions <- function(phyl = NULL, 
                               # Phyloseq object
                               treatvar = NULL, 
                               # Column that contains the treatments
                               treat1 = NULL,
                               treat2 = NULL,
                               OTUcol = "Genus", 
                               # Column that defines OTU
                               return_unfiltered = FALSE,
                               alpha = 0.05,
                               plot = TRUE,
                               fit_type = c("parametric",
                                            "local", 
                                            "mean")) {
  if(is.null(phyl)){
    stop("Phyloseq object missing")
  }
  if(is.null(treat1) | is.null(treat2)){
    stop("Treatment1 and / or treatment2 are missing") 
    # Is this a good way to handle this?
    # What could be improved?
  }
  ## Extract counts, by treatvar (reformulate will turn this into ~Syncom for the default value)
  result_tab = phyl %>% 
    phyloseq::phyloseq_to_deseq2(reformulate(treatvar)) %>% 
    DESeq2::DESeq(fitType = fit_type) %>% 
    DESeq2::results(contrast = c(paste(treatvar), paste(treat1), paste(treat2)),
                    alpha = alpha,
                    pAdjustMethod = "holm",
                    tidy = TRUE) %>% 
    dplyr::left_join(as.data.frame(tax_table(phyl)) %>% 
                       rownames_to_column("row"),
                     by = "row") %>% 
    dplyr::select(-row)
  
  if(return_unfiltered) {
    return(result_tab)
  }
  
  result_tab = result_tab %>%
    dplyr::select(log2FoldChange, lfcSE, padj, eval(OTUcol)) %>% 
    filter(padj < alpha)
  if(plot) {
    result_tab %>% 
      ggplot(aes( x = log2FoldChange, y = Genus)) +
      geom_point() +
      geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.2) +
      geom_vline(aes(xintercept = 0), lty = 3) +
      ylab(element_blank()) +
      theme(text = element_text(family = "NimbusMon"))
  } else {
    result_tab
  }
}