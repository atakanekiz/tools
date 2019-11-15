# Create the KMplotter function to speed up plotting
KMplotter <- function(dat,   # data frame containing samples in rows, genes/grouping variables in columns
                      feature, # gene or variable to plot survival for
                      top_cutoff = 0.5, # choose top percent of the gene expression (quantile)
                      bottom_cutoff = 0.5, # choose bottom percent of gene expression (quantile)
                      categorization_val = NULL, # user can provide a numerical value to binarize gene expression value
                      grouping_var = NULL, # categorical variable that can divide the samples into subgroups (e.g. immune vs keratin subtype)
                      remove_groups = NULL, # discard groups within the grouping_var column that you want to omit (eg MITF low)
                      ... # arguments to pass to ggsurvplot
){
  
  
  require(dplyr)
  require(survival)
  require(survminer)
  require(tibble)
  
  
  
  if(!is.numeric(dat[, feature])){
    
    message("Selected feature is not numeric. Plotting categories instead.")
    
  } else {
    
    if(!is.null(categorization_val)){
      
      
      dat[, feature] <- cut(dat[, feature], breaks=c(-Inf, categorization_val, Inf), 
                            labels = c("Low", "High"))
      
      
    } else {
      
      if(top_cutoff == 0.5 & bottom_cutoff==0.5){
        
        top_segment <- quantile(dat[,feature], top_cutoff)
        
        dat[, feature] <- cut(dat[, feature], breaks=c(-Inf, top_segment, Inf), 
                              labels = c("Low", "High"))
        
        
      } else {
      
      top_segment <- quantile(dat[,feature], top_cutoff)
      bottom_segment <- quantile(dat[,feature], bottom_cutoff)
      
      
      dat[, feature] <- cut(dat[, feature], breaks=c(-Inf, bottom_segment, top_segment, Inf), 
                            labels = c("Low", "Mid", "High"))
      
      dat <- dat[dat[,feature]!= "Mid",] 
      
      }
      
      
    }
    
  }
  
  if(!is.null(grouping_var)){
    
    formula <- as.formula(paste0("Surv(days_to_event, vital_status)~",feature,"+", grouping_var))
    
    n_caption <- data.frame(table(dat[, grouping_var], dat[, feature]))
    
    n_caption <- add_column(n_caption, n="n =", .after = 2)
    
    n_caption <- paste(apply(n_caption, 1, paste, collapse = " "), collapse = "\n")
    
    
  } else {
    
    formula <- as.formula(paste0("Surv(days_to_event, vital_status)~",feature))
    
    n_caption <- summary(dat[, feature])
    
    n_caption <- paste0(names(n_caption), ":", n_caption, collapse = " ")
        
  }
  
  
  sfit <- surv_fit(formula, dat)
  logrank_p <- surv_pvalue(sfit)$pval
  
  ggsurvplot(sfit, pval = T, censor.shape = "|", legend.title = "Categories",
             caption = n_caption,
             font.caption=c(10, "bold.italic", "gray40"), ...)+
    labs(title = feature)
  
  
}
