# SSGSEA_In_R


#In this manner, ssGSEA transforms a single sample's gene expression profile to a gene set enrichment profile. A gene set's enrichment score represents the activity        #level of the biological process in which the gene set's members are coordinately up- or down-regulated.


setwd("E:/Samiksha/M.Sc_Bioinformatics/SEM_3/Cancer_Genomics/Cancer_genomics_practical")

data<-readRDS('log2cpm.rds')

G_list <- getBM(filters = "ensembl_peptide_id", 
                attributes = c("ensembl_peptide_id", "entrezgene", "description"),
                values = data, mart = mart)
                
data[1:5,1:5]
data<-as.matrix(data)
data

gene_set = read.csv("markers2Sep.csv")
head(gene_set)

gene_sets = as.list(as.data.frame(gene_set))
print("genes set ready")
gene_set


ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

res=ssgsea(data,gene_set,norm=F,scale = T)
res
