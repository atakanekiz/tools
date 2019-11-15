# Function to facilitate DE analysis in nanostring experiments (and other experiments such as RNAseq)

# Improvements
## add plot friendly gene column in results
## make summarizedexperiment friendly (pre processing before the function maybe?)
## integrate with hciR package's plotting functions? (plot_pca maybe?)
## heatmap of DE genes
## think about and ensure accuracy of multiple contrast designs
## Think about blocking in design matrix (edgeR manual) to find subtle patterns while controlling for batch/tissue effects (relevant?)
## Blocking/paired analyses for comparing the same samples types (matched healthy and tumor) from individual patients


limmar <- function(dat, # log2 transformed expression data. Genes in rows, samples in columns
                   sample_df, # Data frame with metadata associated with samples. Samples in rows metadata variables are in columns. 
                   # gene_df, # Data frame with metadata associated with genes
                   grouping_var, # Which sample metadata column do you want to group on for differential expression analyses? 
                                 # Create a pasted column with multiple grouping variables if needed (e.g. uninv.ln.brisk, skin.prim.brisk)
                   contrasts, # A character vector of contrasts to be extracted (e.g. c("brisk-nonbrisk", "brisk-absent"))
                   low_expr_filtering = F, # discard lowly expressed genes
                   lfc_testing_cutoff = 1.2,
                   show_summary = T, # Text output of numbers of DE genes
                   de_table = T, # show a table of DE genes
                   de_table_padj_cutoff = 0.05, # Show genes with Padj < 0.05
                   # plot_de_genes = T, # Plot differentially expressed genes between two grps
                   # report_genes =T, # Print upregulated and downregulated genes in console
                   show_venn = F, # Show common DE genes between contrasts
                   venn_subsets, # which DE analyses to compare against each other list of number indices for contrasts (eg list(c(1,2), c(1,3)))
                   extract_res = NULL, # Tabulate the results of DE analysis. Default (null) is the first contrast 
                   interactive_plot = F, # Glimma powered html outputs
                   meta_vars = NULL # For interactive Glimma MDS plot which parameters do you want to color code. A character vector matching colnames
){
  
  
  require(edgeR)
  require(limma)
  # require(ggpubr)
  require(dplyr)
  require(DT)
  
  
  grp <- sample_df[,grouping_var]
  
  design <- model.matrix(~ 0 + grp)
  
  colnames(design) <- gsub("grp", "", colnames(design))
  
  
  
  cont <- makeContrasts(
    contrasts = contrasts,
    levels = colnames(design)
  )
  
  
  # Setup DGE object
  
  dge <- DGEList(counts= dat,
                 samples = sample_df, group = grp)
  
  
  if(low_expr_filtering){
    
    keep <- filterByExpr(dge)   # can add group specification here later. Think about nonuniform gene expression in data subsets.
    dge <- dge[keep, , keep.lib.sizes=F]
    
  }
  
  
  # Linear fit and DE
  
  vfit <- lmFit(dge$counts, design)
  vfit <- contrasts.fit(vfit, contrasts=cont)
  
  
  
  # logFC thresholded results
  tfit <- treat(vfit, lfc=lfc_testing_cutoff)
  dt <- decideTests(tfit)
  
  if(show_summary) print(summary(dt))
  
  
  if(is.null(extract_res)) extract_res <- contrasts[1]
  
  # Save results as a global variable
  de_res <<- topTreat(tfit, coef=extract_res, n=Inf )
  
    
    de_res_sign <- de_res %>%
      add_column("Gene" = rownames(de_res), .after = 0) %>%
      filter(adj.P.Val < de_table_padj_cutoff) %>%
      mutate_if(is.numeric, signif, digits=3) %>%
      arrange(desc(logFC))
    
    sign_genes <<- de_res_sign$Gene
    
    
    # if(de_table) print(DT::datatable(de_res_sign))
  
    # if(plot_de_genes) {
    #   
    #   plot_genes <- de_res_sign$Gene
    #  
    #   dat_subset <- merged %>%
    #     filter(hist.tissue_desc %in% c("uninv_ln", "skin_prim"))
    #   
    #   
    #   ggdotplot(dat_subset, "hist.tissue_desc", plot_genes)
    #    
    # }
    
  
  if(show_venn){
    
    for(i in 1:length(venn_subsets)){
      
      print(vennDiagram(dt[, venn_subsets[[i]]], circle.col = c("red", "blue")))
      
      
    }
    
  }
  

  
  
  if(interactive_plot){
    
    require(Glimma)
    
    glMDSPlot(dge, groups=dge$samples[, meta_vars], launch = F, html= paste(extract_res, "MDSplot"))
    
    glMDPlot(tfit, status = dt[,extract_res], counts=dge$counts, launch = F, groups = dge$samples$group, html = paste(extract_res, "MDplot"))
    
    
  }
  
  
  
    if(de_table) print(htmltools::tagList(DT::datatable(de_res_sign)))
  
  
}




# # Simulate data
# 
# set.seed(123)
# dat <- matrix(rnorm(100000, mean = 20, sd = 2), ncol = 10)
# gene_names <- stringi::stri_rand_strings(dim(dat)[1], 5)
# rownames(dat) <- gene_names
# dat[1:500,1:3] <- dat[1:500,1:3] + 5
# dat[600:1000,4:6] <- dat[600:1000,4:6] - 5
# dat[1:500,7:10] <- dat[1:500,7:10]  + 2
# 
# 
# 
# colnames(dat) <- paste0("runid_", 1:10)
# 
# 
# 
# # Setup gene and sample metadata
# 
# 
# sample_df <- data.frame(var1 = rnorm(10), var2=stringi::stri_rand_strings(10, 3), 
#                         immune = c(rep(c("brisk", "nonbrisk"), each=3), rep("absent", 4))
#                         # group=grp,
# )
# 
# 
# genes <- data.frame(Symbol = gene_names,
#                     Location = rnorm(dim(dat)[2], mean = 1000, sd=500))
# 
# contrasts <- c("nonbrisk-absent", "brisk-nonbrisk", "brisk-absent")
# grouping_var <- "immune"
# 
# venn_subsets <- list(c(1,2), c(1,3))
# 
# show_venn =T




# 
# # SUMMARIZEDEXPERIMENT IMPLEMENTATION?
# 
# sumex <- SummarizedExperiment::SummarizedExperiment(assays = as.matrix(dat),
#                                                     colData = sample_df)
# 
# 
# hciR::plot_pca(sumex, intgroup = "hist.til_score", tooltip = c("clin.gender", 
#                                                                "hist.til_score",
#                                                                "clin.resp_class",
#                                                                "hist.tissue_desc"))
