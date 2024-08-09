##########################################################################################
# Function for Wilcoxon Test with multiple sampling correction and vulcano plot 
#############################################################################################
# Marisa Keller 
# Mai 2021 

# requries 
# library(pROC)
# library(tidyverse)
# library(EnhancedVolcano)

# colnames for the featTable must correspond to the rownames of the metaTable

# Wilcoxon Test ----
calWilcox <- function(featTable, metaTable, fileName,resultfolder, case, ctr){
  # prepare matrix
  #featTable <- featTable[which(rowSums(featTable >= 1e-05) >= 3), ]
  featTable <- featTable[!rownames(featTable)=="-1", , drop = FALSE]
  
  df = featTable[,colnames(featTable) %in% rownames(metaTable)]
  metaTable$ID = rownames(metaTable)
  p.val <- matrix(NA, nrow=nrow(featTable), ncol=1, 
                  dimnames=list(row.names(featTable)))
  fc <- p.val
  aucs.mat <- p.val
  aucs.all  <- vector('list', nrow(featTable))
  confidence.lower  <- p.val
  confidence.higher  <- p.val
  est  <- p.val
  ##############################################################################
  # calculate wilcoxon test and effect size for each feature
  for (f in row.names(featTable)) {
    
    x <- as.numeric(featTable[f, metaTable %>% filter(status==case) %>% pull(ID)])
    y <- as.numeric(featTable[f, metaTable %>% filter(status %in% ctr) %>% pull(ID)])
    
    # Wilcoxon
    res <- wilcox.test(x, y, exact=FALSE)
    p.val[f,1] <- res$p.value
    
    # AUC
    aucs.all[[f]][[1]]  <- c(pROC::roc(controls=y, cases=x, 
                                 direction='<', ci=TRUE, auc=TRUE)$ci)
    aucs.mat[f,1] <- c(pROC::roc(controls=y, cases=x, 
                           direction='<', ci=TRUE, auc=TRUE)$ci)[2]
    
    
    # FC 
    q.p <- quantile(x, probs=seq(.1, .9, .05))
    q.n <- quantile(y, probs=seq(.1, .9, .05))
    fc[f,1] <- sum(q.p - q.n)/length(q.p)
  }
  ##############################################################################
  # fdr correction
  p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method="fdr"),
                      check.names = FALSE)
  colnames(p.adj) <- "adj"
  
  # add fc and auc
  p.adj$p.val <- p.val
  p.adj$fc <- fc
  p.adj$auc.mat <- aucs.mat
  
  p.adj$log10p = -log10(as.numeric(p.adj$adj))
  p.adj <- as.data.frame(p.adj)
  rownames(p.adj) <- make.names(rownames(df), unique=TRUE)
  # add average relative abundance for each taxa 
  p.adj$rel <- rowSums(df)/nrow(df)
  #print(head(p.adj))
  p.adj <- p.adj %>% na.omit(p.adj)
  p.adj$significant <- ifelse(p.adj$adj < 0.01, "p.adj < 0.01", "not sig")
  p.adj$species <-rownames(p.adj)

  # save file
  write.table(p.adj, file=paste0(resultfolder, fileName, 'p.adj.tsv'), 
              sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
  return(p.adj)
}

