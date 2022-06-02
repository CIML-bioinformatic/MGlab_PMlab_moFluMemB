# #################################################################################################
# This script aims to prepare some figures of the final article
# #################################################################################################



## @knitr prepare_paper_figure_heatmap

# Build a dataframe of mean of expression of cluster of each tissue for the selected marker genes
# and a dataframe of percentage of cell with positive expression for the selected marker genes
marker_genes_cluster_means_df = data.frame()
marker_genes_cluster_pct_df = data.frame()
marker_genes_cluster_df_colnames = vector()
marker_genes_cluster_df_coltissues = vector()
for( current_tissue in names( PATH_VARIABLE_GENES_EXPRESSION_FILE)){
  
  cluster_mapping = cluster_mapping_list[[ current_tissue]]
  for( current_cluster in levels( cluster_mapping$cluster)){
    cluster_cell_set = cluster_mapping[ which( cluster_mapping$cluster == current_cluster), "cell.id"]
    expression_df = variable_genes_expression_list[[ current_tissue]][ MARKER_GENES_LIST, cluster_cell_set]
    means_set = rowMeans( expression_df)
    pct_set = apply( expression_df, 1, function( row){
      return( length( which( row > 0))/length( row))
    })
    marker_genes_cluster_means_df = rbind( marker_genes_cluster_means_df, t( data.frame( means_set)))
    marker_genes_cluster_pct_df = rbind( marker_genes_cluster_pct_df, t( data.frame( pct_set)))
    marker_genes_cluster_df_colnames = append( marker_genes_cluster_df_colnames, paste0( current_tissue, " C", current_cluster))
    marker_genes_cluster_df_coltissues = append( marker_genes_cluster_df_coltissues, current_tissue)
  }
}

# transpose the dataframe to have the genes in rows
marker_genes_cluster_means_df = data.frame( t( marker_genes_cluster_means_df))
names( marker_genes_cluster_means_df) = marker_genes_cluster_df_colnames

marker_genes_cluster_pct_df = data.frame( t( marker_genes_cluster_pct_df))
names( marker_genes_cluster_pct_df) = marker_genes_cluster_df_colnames

# Remove the cluster in Lung that was affected by dissociation process (LG cluster 2)
marker_genes_cluster_means_df =  marker_genes_cluster_means_df[ , - which( names( marker_genes_cluster_means_df) == "LG C2")]
marker_genes_cluster_pct_df =  marker_genes_cluster_pct_df[ , - which( names( marker_genes_cluster_pct_df) == "LG C2")]

# Plot the heatmaps in SVG files
annotation_col = data.frame( Tissue = marker_genes_cluster_df_coltissues)
row.names( annotation_col) = marker_genes_cluster_df_colnames

svg( file = file.path( PATH_ANALYSIS_OUTPUT, "Heatmap_MarkerGenes_MeanExpression2.svg") )
pheatmap::pheatmap( marker_genes_cluster_means_df,
                    cluster_rows = FALSE, cluster_cols= FALSE,
                    annotation_col = annotation_col,
                    cellwidth = 10,
                    gaps_col = c( 3,6),
                    color = viridis(20, option = "B"),
                    scale = "row")
dev.off()

svg( file = file.path( PATH_ANALYSIS_OUTPUT, "Heatmap_MarkerGenes_MeanExpression.svg") )
pheatmap::pheatmap( marker_genes_cluster_means_df,
                    cluster_rows = FALSE, cluster_cols= FALSE,
                    annotation_col = annotation_col,
                    cellwidth = 10,
                    gaps_col = c( 3,6),
                    color = viridis(20, option = "B"),
                    scale = "row")
dev.off()

svg( file = file.path( PATH_ANALYSIS_OUTPUT, "Heatmap_MarkerGenes_PctExpression.svg") )
pheatmap::pheatmap( marker_genes_cluster_pct_df,
                    cluster_rows = FALSE, cluster_cols= FALSE,
                    annotation_col = annotation_col,
                    cellwidth = 10,
                    gaps_col = c( 3,6),
                    color = viridis(20, option = "B"))
dev.off()

# test_cell_list = cluster_mapping_list[[ "SP"]][ which( cluster_mapping_list[[ "SP"]]$cluster == 2), "cell.id"]
# length( which( variable_genes_expression_list[[ "SP"]][ "Cd22", test_cell_list] > 0))
# length( variable_genes_expression_list[[ "SP"]][ "Cd22", test_cell_list])
