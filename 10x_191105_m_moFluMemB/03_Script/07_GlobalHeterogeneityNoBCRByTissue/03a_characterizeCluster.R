# ##############################################################################
# This script aims to produce analysis to caracterise clusters fourn in
# previous step
# ##############################################################################


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# PAIRWISE MARKER GENES - IDENTIFICATION
# ---------------------------------------------------------------
# ---------------------------------------------------------------

## @knitr characterizeCluster_pairwiseMarkers

Idents( sc10x.rna.seurat) <- "seurat_clusters"

# Create the pairwise comparison by combinations of clusters
seurat_clusters <- sort( unique( Idents( sc10x.rna.seurat)))
pairwise <- combn( seurat_clusters, 2)

# Accumulate the comparisons results in a single dataframe
pairwise_marker_genes_df = data.frame()
for(i in 1:ncol(pairwise)) {
  # find the markers genes if first direction
  pairwise_marker_genes = FindMarkers( sc10x.rna.seurat, 
                                       ident.1 = pairwise[1, i], 
                                       ident.2 = pairwise[2, i],
                                       test.use        = FINDMARKERS_METHOD,
                                       only.pos        = FINDMARKERS_ONLYPOS,
                                       min.pct         = FINDMARKERS_MINPCT,
                                       logfc.threshold = FINDMARKERS_LOGFC_THR,
                                       verbose         = .VERBOSE);
  
  # find the markers genes if second direction
  pairwise_marker_genes_switch = FindMarkers( sc10x.rna.seurat, 
                                              ident.1 = pairwise[2, i], 
                                              ident.2 = pairwise[1, i],
                                              test.use        = FINDMARKERS_METHOD,
                                              only.pos        = FINDMARKERS_ONLYPOS,
                                              min.pct         = FINDMARKERS_MINPCT,
                                              logfc.threshold = FINDMARKERS_LOGFC_THR,
                                              verbose         = .VERBOSE);
  
  # For direct test, FoEliminate the non-significant genes and keep only 100 best
  pairwise_marker_genes = pairwise_marker_genes[ which( pairwise_marker_genes$p_val_adj < FINDMARKERS_PVAL_THR), ]
  if( nrow( pairwise_marker_genes) > 0){
    # Reorder results by adjusted p-value
    pairwise_marker_genes = pairwise_marker_genes[ order( pairwise_marker_genes[ ,"p_val_adj"], decreasing = FALSE), ]
    # Eliminate row names to allow accumulations in a single dataframe
    pairwise_marker_genes$gene.name = row.names( pairwise_marker_genes)
    row.names( pairwise_marker_genes) = seq( 1, nrow( pairwise_marker_genes), 1)
    # Add the information about compared clusters
    pairwise_marker_genes$cluster1 = pairwise[1, i]
    pairwise_marker_genes$cluster2 = pairwise[2, i]
    # Acculumulate the results
    pairwise_marker_genes_df = rbind( pairwise_marker_genes_df, pairwise_marker_genes)
  }
  
  # For switched test, Eliminate the non-significant genes and keep only 100 best
  pairwise_marker_genes_switch = pairwise_marker_genes_switch[ which( pairwise_marker_genes_switch$p_val_adj < FINDMARKERS_PVAL_THR), ]
  if( nrow( pairwise_marker_genes_switch) > 0){
    # Reorder results by adjusted p-value
    pairwise_marker_genes_switch = pairwise_marker_genes_switch[ order( pairwise_marker_genes_switch[ ,"p_val_adj"], decreasing = FALSE), ]
    # Eliminate row names to allow accumulations in a single dataframe
    pairwise_marker_genes_switch$gene.name = row.names( pairwise_marker_genes_switch)
    row.names( pairwise_marker_genes_switch) = seq( 1, nrow( pairwise_marker_genes_switch), 1)
    # Add the information about compared clusters
    pairwise_marker_genes_switch$cluster1 = pairwise[2, i]
    pairwise_marker_genes_switch$cluster2 = pairwise[1, i]
    # Acculumulate the results
    pairwise_marker_genes_df = rbind( pairwise_marker_genes_df, pairwise_marker_genes_switch)
  }
  
}

# Reorder columns,rename them and convert numbers to simplify reading of results
pairwise_marker_genes_df = pairwise_marker_genes_df[ , c( "cluster1", "cluster2", "gene.name", "avg_logFC", "pct.1", "pct.2",  "p_val", "p_val_adj")]
pairwise_marker_genes_df$p_val = signif( pairwise_marker_genes_df$p_val, 3)
pairwise_marker_genes_df$avg_logFC = signif( pairwise_marker_genes_df$avg_logFC, 3)
pairwise_marker_genes_df$p_val_adj = signif( pairwise_marker_genes_df$p_val_adj, 3)
names( pairwise_marker_genes_df) = c( "cluster", "vs.cluster", "gene.name", "avg_logFC", "pct.1", "pct.2", "p_val", "p_val_adj")


# Display the table of pairwise marker genes
datatable( pairwise_marker_genes_df, caption = "Pairwise marker genes", rownames = FALSE)

# Write to file the dataframe of pairwise marker genes
pairwise_marker_genes_file_path = file.path( PATH_ANALYSIS_OUTPUT, paste0( "pairwise_marker_genes_", CURRENT_TISSUE, ".tsv"))
write.table( pairwise_marker_genes_df, 
             file = pairwise_marker_genes_file_path,
             quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

# Summarize information per cluster and per gene by associating each gene to the cluster it is marker of
# and the list of clusters against which it has been detected
pairwise_topMarkers = data.frame()
pairwise_markers_summary_df = data.frame()
for( current_cluster in seurat_clusters){
  current_df = pairwise_marker_genes_df[ which( pairwise_marker_genes_df$cluster == current_cluster), ]
  new_df = data.frame()
  for( current_gene in unique( current_df$gene.name)){
    marked_cluster_set = as.character( sort( current_df[ which( current_df$gene.name == current_gene), "vs.cluster"]), decreasing = FALSE)
    new_df = rbind( new_df, data.frame( gene.name = current_gene,
                                        cluster = current_cluster,
                                        vs.clusters = paste( marked_cluster_set, collapse = ", "),
                                        vs.clusters.nb = length( marked_cluster_set)))
  }
  if( nrow( new_df) > 0){
    new_df = new_df[ order( new_df[ ,"vs.clusters.nb"], decreasing = TRUE), ]
    pairwise_markers_summary_df = rbind( pairwise_markers_summary_df, new_df)
    pairwise_topMarkers = rbind( pairwise_topMarkers, head( new_df, n = FINDMARKERS_SHOWTOP))
  }
}

# Display the summary table of pairwise marker genes
datatable( pairwise_markers_summary_df, caption = "Summary of pairwise marker genes", rownames = FALSE)

# Export the summary table to file
pairwise_markers_summary_file_path = file.path( PATH_ANALYSIS_OUTPUT, paste0( "pairwise_marker_genes_summary_", CURRENT_TISSUE, ".tsv"))
write.table( pairwise_markers_summary_df, 
             file = pairwise_markers_summary_file_path,
             quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

# Export the topMarkers summary table to file
pairwise_markers_summary_topMarkers_file_path = file.path( PATH_ANALYSIS_OUTPUT, paste0( "pairwise_marker_genes_summary_topMarkers_", CURRENT_TISSUE, ".tsv"))
write.table( pairwise_topMarkers, 
             file = pairwise_markers_summary_topMarkers_file_path,
             quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

