#' Compare conditions. 
#' This function compares the relative abundance between two treatments, via DESeq2.
#' @details This function takes a phyloseq object and performs a OTU-level comparison between two conditions returns log2FC, does filtering on adjusted pval, controlled via alpha
#' @param phyl a phyloseq object
#' @param treatvar character, the treatment variable. Has to be part of the phyloseq object sample_data
#' @param treat1 first treatment; will be numerator
#' @param treat2 second treatment; will be denominator
#' @param OTUcol which OTU assignment level; default "Genus"
#' @param split_var a variable used to split the dataframe. If set this will loop through the levels of the split var
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
                               split_var = NULL,
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
  if(!is.null(split_var)){
    cat(paste("split_var is", eval(split_var), "\n"))
    samp_dat <- sample_data(phyl) %>% 
      as("data.frame")
    split_levels <-  samp_dat[,(eval(split_var))] %>% 
      unique()
    result_tab <- data.frame()
    for(curr_level in split_levels) 
    {
      cat("Looping ...\n")
      cat(paste(c(curr_level,"\n")))
      tmp_dat <- phyl %>% 
      sample_data() %>% 
      as("data.frame") 
      tmp_dat$split_by = tmp_dat[,eval(split_var)]
      
      tmp_dat <- tmp_dat %>% dplyr::filter(split_by == paste(curr_level))
      
      tmp_phyl <- phyl
      sample_data(tmp_phyl) <- sample_data(tmp_dat)
      
      tmp_results <- tmp_phyl %>% 
      phyloseq::phyloseq_to_deseq2(reformulate(treatvar)) %>% 
      DESeq2::DESeq(fitType = fit_type) %>% 
      DESeq2::results(contrast = c(paste(treatvar[1]), paste(treat1), paste(treat2)),
                      alpha = alpha,
                      pAdjustMethod = "holm",
                      tidy = TRUE) %>% 
      dplyr::left_join(as.data.frame(tax_table(phyl)) %>% 
                         rownames_to_column("row"),
                       by = "row") %>% 
      dplyr::select(-row) 
      tmp_results[eval(split_var)] <- rep(curr_level, nrow(tmp_results))
      result_tab <- rbind(result_tab, tmp_results)
    }
  } else{

  ## Extract counts, by treatvar
  result_tab = phyl %>% 
    phyloseq::phyloseq_to_deseq2(reformulate(treatvar)) %>% 
    DESeq2::DESeq(fitType = fit_type) %>% 
    DESeq2::results(contrast = c(paste(treatvar[1]), paste(treat1), paste(treat2)),
                    alpha = alpha,
                    pAdjustMethod = "holm",
                    tidy = TRUE) %>% 
    dplyr::left_join(as.data.frame(tax_table(phyl)) %>% 
                       rownames_to_column("row"),
                     by = "row") %>% 
    dplyr::select(-row)
  }
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