# ##############################################################################
# This script aims to analyze the clonotypes groups and their dispersion
# over the tissues/clusters in relation with data from Seurat analysis
# ##############################################################################

## @knitr analyze_clonotypes

# ................................................................................................
## CLONOTYPES versus SAMPLES
# ................................................................................................

cat("<HR><H5>CLONOTYPES versus SAMPLES</H5>")

# Show the dispersion of largest clonotypes on UMAP embedding
large.clonotype.similarity.symbol_plot_list = list()
for( large.clonotype.similarity.symbol in unique( clonotype_seurat_df$large.clonotype.similarity.symbol)){
  large.clonotype.similarity.symbol_plot_list[[ large.clonotype.similarity.symbol]]= 
    ggplot() + 
             geom_point( data = clonotype_seurat_df, aes( x=UMAP_1, y=UMAP_2), col="lightgrey") +
             geom_point( data = clonotype_seurat_df[ which( clonotype_seurat_df$large.clonotype.similarity.symbol == large.clonotype.similarity.symbol),], aes( x=UMAP_1, y=UMAP_2, col=seurat_clusters)) +
             scale_color_manual( values = seurat_cluster_palette) +
             ggtitle( large.clonotype.similarity.symbol) + 
             theme_minimal() + 
             theme( legend.position = "None", strip.text.x = element_text(size = 5)) +
             theme( axis.title = element_text(size = 5)) +
             theme(plot.title = element_text(size = 8))
}

grid.arrange( grobs = large.clonotype.similarity.symbol_plot_list, ncol = 3)

# Show the dispersion of all clonotypes over samples in a datatable
clonotype_sample_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "full.clonotype.similarity.symbol", "sample")]))
clonotype_sample_table$total = apply( clonotype_sample_table, 1, sum)
print( htmltools::tagList( datatable( clonotype_sample_table, caption = "Dispersion of all clonotypes over samples")))


# ................................................................................................
## CLONOTYPES versus TISSUE (when several tissues)
# ................................................................................................

if( length( unique( clonotype_seurat_df$tissue)) > 1){
  
  cat("<HR><H5>LARGEST CLONOTYPES versus TISSUE</H5>")
  
  # Show the dispersion of cells in large clonotype over samples in a datatable
  large_clonotype_tissue_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "tissue", "large.clonotype.similarity.symbol")]))
  large_clonotype_tissue_table_line_total = apply( large_clonotype_tissue_table, 1, sum)
  large_clonotype_tissue_table_column_total = apply( large_clonotype_tissue_table, 2, sum)
  large_clonotype_tissue_table = large_clonotype_tissue_table[ large_clonotype_tissue_table_line_total > 0, large_clonotype_tissue_table_column_total > 0]
  print( htmltools::tagList( datatable( t(large_clonotype_tissue_table), caption="Dispersion of largest clonotypes over tissues")))
  
  # Show the dispersion of samples over largest clonotypes in percentage in a cumulative barplot
  large_clonotype_tissue_wide = large_clonotype_tissue_table
  large_clonotype_tissue_wide = 100*large_clonotype_tissue_wide / large_clonotype_tissue_table_line_total
  large_clonotype_tissue_wide$large.clonotype.similarity.symbol = row.names( large_clonotype_tissue_wide)
  large_clonotype_tissue_wide = reshape2::melt( large_clonotype_tissue_wide, id.vars = c( "large.clonotype.similarity.symbol"))
  print( ggplot( large_clonotype_tissue_wide) +
           geom_bar( aes( x = large.clonotype.similarity.symbol, fill = variable, y=value), stat="identity") + 
           scale_fill_manual( values = large_clonotype_palette) +
           labs( fill = "Clonotype", x= "Tissue", y = "Percentage") + 
           theme_minimal() +
           ggtitle( "Percentage of cell in tissue for each clonotype")
  )
  
}


# ................................................................................................
## LARGEST CLONOTYPES versus SAMPLES
# ................................................................................................

cat("<HR><H5>LARGEST CLONOTYPES versus SAMPLES</H5>")

# Show the dispersion of cells in large clonotype over samples in a datatable
large_clonotype_sample_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "sample", "large.clonotype.similarity.symbol")]))
large_clonotype_sample_table_line_total = apply( large_clonotype_sample_table, 1, sum)
large_clonotype_sample_table_column_total = apply( large_clonotype_sample_table, 2, sum)
large_clonotype_sample_table = large_clonotype_sample_table[ large_clonotype_sample_table_line_total > 0, large_clonotype_sample_table_column_total > 0]
print( htmltools::tagList(datatable( t(large_clonotype_sample_table), caption="Dispersion of largest clonotypes over samples")))

# Show the dispersion of samples over largest clonotypes in percentage in a cumulative barplot
large_clonotype_sample_wide = large_clonotype_sample_table
large_clonotype_sample_wide = 100*large_clonotype_sample_wide / large_clonotype_sample_table_line_total
large_clonotype_sample_wide$large.clonotype.similarity.symbol = row.names( large_clonotype_sample_wide)
large_clonotype_sample_wide = reshape2::melt( large_clonotype_sample_wide, id.vars = c( "large.clonotype.similarity.symbol"))
print( ggplot( large_clonotype_sample_wide) +
         geom_bar( aes( x = large.clonotype.similarity.symbol, fill = variable, y=value), stat="identity") + 
         scale_fill_manual( values = large_clonotype_palette) +
         labs( fill = "Clonotype", x= "Sample", y = "Percentage") + 
         theme_minimal() +
         ggtitle( "Percentage of cell in sample for each large clonotype")
)

# ................................................................................................
## LARGEST CLONOTYPES versus CLUSTERS
# ................................................................................................

cat("<HR><H5>LARGEST CLONOTYPES versus CLUSTERS</H5>")

# Show the dispersion of cells in large clonotype over clusters in a datatable
large_clonotype_cluster_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "seurat_clusters", "large.clonotype.similarity.symbol")]))
large_clonotype_cluster_table_line_total = apply( large_clonotype_cluster_table, 1, sum)
large_clonotype_cluster_table_column_total = apply( large_clonotype_cluster_table, 2, sum)
large_clonotype_cluster_table = large_clonotype_cluster_table[ large_clonotype_cluster_table_line_total > 0, large_clonotype_cluster_table_column_total > 0]
print( htmltools::tagList( datatable( t( large_clonotype_cluster_table), caption = "Dispersion of largest clonotypes over clusters")))

# Show the dispersion of clusters over largest clonotypes in percentage in a cumulative barplot
large_clonotype_cluster_wide = large_clonotype_cluster_table
large_clonotype_cluster_wide = 100*large_clonotype_cluster_wide / large_clonotype_cluster_table_line_total
large_clonotype_cluster_wide$large.clonotype.similarity.symbol = row.names( large_clonotype_cluster_wide)
large_clonotype_cluster_wide = reshape2::melt( large_clonotype_cluster_wide, id.vars = c( "large.clonotype.similarity.symbol"))
print( ggplot( large_clonotype_cluster_wide) +
         geom_bar( aes( x = large.clonotype.similarity.symbol, fill = variable, y=value), stat="identity") + 
         scale_fill_manual( values = large_clonotype_palette) +
         labs( fill = "Clonotype", x= "Cluster", y = "Percentage") + 
         theme_minimal() + 
         ggtitle( "Percentage of cell in cluster for each large clonotype")
)

# #  For LG tissue plot the expression of cells from Large Clonotype on marker genes of cluster 0 and 1
# if( exists( "current_tissue") && current_tissue == "LG"){
# 
#   # Start by looking at classical marker genes
#   # ###########################################
#   
#   # Filter markers by cluster
#   topMarkers = by( CLUSTER_MARKER_GENES_DF, CLUSTER_MARKER_GENES_DF[["cluster"]], function(x)
#   {
#     # Filter markers based on adjusted PValue
#     x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
#     # Sort by decreasing logFC
#     x = x[ order(abs(x[["avg_logFC"]]), decreasing = TRUE), , drop = FALSE ]
#     # Return top ones
#     return( if(is.null( FINDMARKERS_SHOWTOP)) x else head( x, n = FINDMARKERS_SHOWTOP));
#   });
#   
#   # Get the gene markers of cluster 0 and 1
#   cluster_marker_genes_df = do.call( rbind, topMarkers);
#   cluster_marker_genes_df = rbind( cluster_marker_genes_df, CLUSTER_MARKER_GENES_DF[ which( CLUSTER_MARKER_GENES_DF$gene == "Cxcr3"),])
#   cluster_marker_genes = as.character( cluster_marker_genes_df[ which( cluster_marker_genes_df$cluster %in% c(0,1)), "gene"])
#   cluster_marker_genes = c( cluster_marker_genes, "Cxcr3")
#   
#   # Plot a heatmap of expression on marker genes of all cells
#   # ..............................................................................................
#   
#   # Define the annotation of genes
#   annotation_row = data.frame( Marker.cluster = as.factor( cluster_marker_genes_df[ which( cluster_marker_genes_df$gene %in% cluster_marker_genes), "cluster"]))
#   row.names( annotation_row) = cluster_marker_genes
# 
#   # Get the expression of selected cells on marker genes
#   cluster_markers_clonotype_expression_df = GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_marker_genes, ]
#   
#   # Get the annotation of cells in clusters and clonotype (for verification)
#   annotation_cell = data.frame( Cell.cluster = clonotype_seurat_df[ colnames( cluster_markers_clonotype_expression_df), "seurat_clusters"],
#                                 Cell.clonotype = clonotype_seurat_df[ colnames( cluster_markers_clonotype_expression_df), "large.clonotype.similarity.symbol"])
#   row.names( annotation_cell) = colnames( cluster_markers_clonotype_expression_df)
# 
#   # Define the colors of annotations
#   annotation_colors = list( Cell.cluster=seurat_cluster_palette,
#                             Marker.cluster=seurat_cluster_palette)
# 
#   # Plot a heatmap of the expression of selected cells on marker genes
#   pheatmap( cluster_markers_clonotype_expression_df,
#             annotation_col = annotation_cell,
#             annotation_row = annotation_row,
#             annotation_colors = annotation_colors,
#             cluster_rows = FALSE,
#             cluster_cols = TRUE,
#             show_colnames = FALSE,
#             main = paste( "Expression of cells on marker genes of clusters 0 and 1"))
#   
#   # Compare the distribution of some genes between clusters and through clonotype and non clonotype
#   # ..............................................................................................
#   for( current_gene in c( "Plac8")){
#     
#     plot_df = data.frame( Gene = cluster_markers_clonotype_expression_df[ current_gene, ],
#                           Cell.cluster = annotation_cell[ colnames( cluster_markers_clonotype_expression_df), "Cell.cluster"],
#                           In.clonotype = !is.na( annotation_cell[ colnames( cluster_markers_clonotype_expression_df), "Cell.clonotype"])
#                           )
#     print( ggplot(plot_df, aes(x=Cell.cluster, y=Gene)) +
#             geom_violin(aes(fill=In.clonotype)) +
#             xlab( "Cluster") + ylab( current_gene) + theme_classic()
#     )
#   }
#   
#   # For all the cells in large clonotypes, show a heatmap of their expression
#   # ................................................................................
#   
#   # -- get the cells associated to the large clonotype
#   current_clonotype_seurat_df = clonotype_seurat_df[ which( !is.na( clonotype_seurat_df$large.clonotype.similarity.symbol)), ]
#   current_clonotype_seurat_df = current_clonotype_seurat_df[ order( current_clonotype_seurat_df$seurat_cluster), ]
#   current_clonotype_cells = row.names( current_clonotype_seurat_df)
#   
#   # -- get the annotation of cells in clusters and clonotype (for verification)
#   annotation_cell = data.frame( Cell.cluster = clonotype_seurat_df[ current_clonotype_cells, "seurat_clusters"],
#                                 Cell.clonotype = clonotype_seurat_df[ current_clonotype_cells, "full.clonotype.similarity.symbol"])
#   row.names( annotation_cell) = current_clonotype_cells
#   
#   # -- get the expression of selected cells on marker genes
#   cluster_markers_clonotype_expression_df = GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_marker_genes, current_clonotype_cells]
#   
#   # -- define the colors of annotations
#   annotation_colors = list( Cell.cluster=seurat_cluster_palette,
#                             Marker.cluster=seurat_cluster_palette)
#   
#   # -- plot a heatmap of the expression of selected cells on marker genes
#   pheatmap( cluster_markers_clonotype_expression_df,
#             annotation_col = annotation_cell,
#             annotation_row = annotation_row,
#             annotation_colors = annotation_colors,
#             cluster_rows = FALSE,
#             cluster_cols = FALSE,
#             show_colnames = FALSE,
#             main = paste( "Expression of cells in all large clonotypes\non marker genes of clusters 0 and 1"))
# 
#   # For each large clonotype, plot the heatmap of gene expression on each clonotype cells
#   # ................................................................................
#   
#   for( current_clonotype in levels( clonotype_seurat_df$large.clonotype.similarity.symbol)){
#     if( !is.na( current_clonotype)){
#       # Get the cells associated to the clonotype
#       current_clonotype_seurat_df = clonotype_seurat_df[ which( clonotype_seurat_df$large.clonotype.similarity.symbol == current_clonotype), ]
#       current_clonotype_seurat_df = current_clonotype_seurat_df[ order( current_clonotype_seurat_df$seurat_cluster), ]
#       current_clonotype_cells = row.names( current_clonotype_seurat_df)
#       
#       # Get the annotation of cells in clusters and clonotype (for verification)
#       annotation_cell = data.frame( Cell.cluster = clonotype_seurat_df[ current_clonotype_cells, "seurat_clusters"],
#                                     Cell.clonotype = clonotype_seurat_df[ current_clonotype_cells, "full.clonotype.similarity.symbol"])
#       row.names( annotation_cell) = current_clonotype_cells
#       
#       # Get the expression of selected cells on marker genes
#       cluster_markers_clonotype_expression_df = GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_marker_genes, current_clonotype_cells]
#       
#       # Define the colors of annotations
#       annotation_colors = list( Cell.cluster=seurat_cluster_palette,
#                                 Marker.cluster=seurat_cluster_palette)
#       
#       # Plot a heatmap of the expression of selected cells on marker genes
#       pheatmap( cluster_markers_clonotype_expression_df,
#                 annotation_col = annotation_cell,
#                 annotation_row = annotation_row,
#                 annotation_colors = annotation_colors,
#                 cluster_rows = FALSE,
#                 cluster_cols = FALSE,
#                 show_colnames = FALSE,
#                 main = paste( "Expression of cells in clonotype", current_clonotype, "\non marker genes of clusters 0 and 1"))
#       
#       # Show in datatable the expression of cells
#       cell_in_cluster_0 = row.names( annotation_cell)[ which( annotation_cell$Cell.cluster == 0)]
#       cell_in_cluster_1 = row.names( annotation_cell)[ which( annotation_cell$Cell.cluster == 1)]
#       if( length( cell_in_cluster_0) > 0 && length( cell_in_cluster_1) > 0){
#         cluster_markers_clonotype_expression_cluster_0_df = as.data.frame( GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_marker_genes, cell_in_cluster_0])
#         cluster_markers_clonotype_expression_cluster_0_df$mean = apply( cluster_markers_clonotype_expression_cluster_0_df, 1, mean)
#         cluster_markers_clonotype_expression_cluster_1_df = as.data.frame( GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_marker_genes, cell_in_cluster_1])
#         cluster_markers_clonotype_expression_cluster_1_df$mean = apply( cluster_markers_clonotype_expression_cluster_1_df, 1, mean)
#         print( htmltools::tagList( datatable( cluster_markers_clonotype_expression_cluster_0_df, caption="Expression of cells in cluster 0")))
#         print( htmltools::tagList( datatable( cluster_markers_clonotype_expression_cluster_1_df, caption="Expression of cells in cluster 1")))
#       }
#     }
#   }
#   
#   # Next, look at pairwise marker genes
#   # ###########################################
#   
#   # Filter markers by cluster
#   topPairwiseMarkers = by( CLUSTER_PAIRWISE_MARKER_GENES_DF, CLUSTER_PAIRWISE_MARKER_GENES_DF[["cluster1"]], function(x)
#   {
#     # Filter markers based on adjusted PValue
#     x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
#     # Sort by decreasing logFC
#     x = x[ order(abs(x[["avg_logFC"]]), decreasing = TRUE), , drop = FALSE ]
#     # Return top ones
#     return( x);
#   });
#   
#   # Get the gene markers of cluster 0 and 1
#   cluster_pairwise_marker_genes_df = do.call( rbind, topPairwiseMarkers);
#   cluster_pairwise_marker_genes_df = rbind( cluster_pairwise_marker_genes_df, CLUSTER_PAIRWISE_MARKER_GENES_DF[ which( CLUSTER_PAIRWISE_MARKER_GENES_DF$gene == "Cxcr3"),])
#   cluster_pairwise_marker_genes = as.character( cluster_pairwise_marker_genes_df[ which( cluster_pairwise_marker_genes_df$cluster1 %in% c(0,1) & cluster_pairwise_marker_genes_df$cluster2 %in% c(0,1)), "gene"])
#   cluster_pairwise_marker_genes = c( cluster_pairwise_marker_genes, "Cxcr3")
#   
#   # Get the expression of selected cells on marker genes
#   cluster_pairwise_markers_clonotype_expression_df = GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_pairwise_marker_genes, ]
#   
#   # Get the annotation of cells in clusters and clonotype (for verification)
#   annotation_cell = data.frame( Cell.cluster = clonotype_seurat_df[ colnames( cluster_pairwise_markers_clonotype_expression_df), "seurat_clusters"])
#   row.names( annotation_cell) = colnames( cluster_pairwise_markers_clonotype_expression_df)
#   
#   # Define the colors of annotations
#   annotation_colors = list( Cell.cluster=seurat_cluster_palette,
#                             Marker.cluster=seurat_cluster_palette)
#   
#   # Plot a heatmap of the expression of selected cells on marker genes
#   pheatmap( cluster_pairwise_markers_clonotype_expression_df,
#             annotation_col = annotation_cell,
#             annotation_colors = annotation_colors,
#             cluster_rows = FALSE,
#             cluster_cols = TRUE,
#             show_colnames = FALSE,
#             main = paste( "Expression of cells on pariwise marker genes of clusters 0 and 1"))
#   
#   # For each clonotype,
#   for( current_clonotype in levels( clonotype_seurat_df$large.clonotype.similarity.symbol)){
#     if( !is.na( current_clonotype)){
#       # Get the cells associated to the clonotype
#       current_clonotype_seurat_df = clonotype_seurat_df[ which( clonotype_seurat_df$large.clonotype.similarity.symbol == current_clonotype), ]
#       current_clonotype_seurat_df = current_clonotype_seurat_df[ order( current_clonotype_seurat_df$seurat_cluster), ]
#       clonotype_cells = row.names( current_clonotype_seurat_df)
#       
#       # Get the annotation of cells in clusters and clonotype (for verification)
#       annotation_cell = data.frame( Cell.cluster = clonotype_seurat_df[ clonotype_cells, "seurat_clusters"],
#                                     Clonotype = clonotype_seurat_df[ clonotype_cells, "full.clonotype.similarity.symbol"])
#       row.names( annotation_cell) = clonotype_cells
#       
#       # Define the colors of annotations
#       annotation_colors = list( Cell.cluster=seurat_cluster_palette,
#                                 Marker.cluster=seurat_cluster_palette)
#       
#       # Get the expression of selected cells on marker genes
#       cluster_pairwise_markers_clonotype_expression_df = GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_pairwise_marker_genes, clonotype_cells]
#       
#       # Plot a heatmap of the expression of selected cells on marker genes
#       pheatmap( cluster_pairwise_markers_clonotype_expression_df,
#                 annotation_col = annotation_cell,
#                 annotation_colors = annotation_colors,
#                 cluster_rows = TRUE,
#                 cluster_cols = FALSE,
#                 show_colnames = FALSE,
#                 main = paste( "Expression of cells in clonotype", current_clonotype, "\non pairwise marker genes of clusters 0 and 1"))
#       
#       # Show in datatable the expression of cells
#       cluster_pairwise_markers_clonotype_expression_cluster_0_df = GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_pairwise_marker_genes, row.names( annotation_cell)[ which( annotation_cell$Cell.cluster == 0)]]
#       cluster_pairwise_markers_clonotype_expression_cluster_1_df = GetAssayData( custom.rna.seurat, assay="RNA", slot="data")[ cluster_pairwise_marker_genes, row.names( annotation_cell)[ which( annotation_cell$Cell.cluster == 1)]]
#       print( htmltools::tagList( datatable( as.matrix( cluster_pairwise_markers_clonotype_expression_cluster_0_df), caption="Expression of cells in cluster 0")))
#       print( htmltools::tagList( datatable( as.matrix( cluster_pairwise_markers_clonotype_expression_cluster_1_df), caption="Expression of cells in cluster 1")))
#     }
#   }
# }

# ................................................................................................
## ISOTYPE HEAVY-CHAIN versus CLUSTERS
# ................................................................................................

cat("<HR><H5>ISOTYPE HEAVY-CHAIN versus CLUSTERS</H5>")

# Show the dispersion of isotypes of heavy chain over clusters in a datatable
isotype_heavy_cluster_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "seurat_clusters", "isotype.heavy")]))
isotype_heavy_cluster_table_line_total = apply( isotype_heavy_cluster_table, 1, sum)
isotype_heavy_cluster_table_column_total = apply( isotype_heavy_cluster_table, 2, sum)
isotype_heavy_cluster_table = isotype_heavy_cluster_table[ isotype_heavy_cluster_table_line_total > 0, isotype_heavy_cluster_table_column_total > 0]
print( htmltools::tagList( datatable( t(isotype_heavy_cluster_table), caption = "Dispersion of heavy chain isotypes over clusters (all cells)")))

# Show the dispersion of clusters over isotype of heavy chain in percentage in a cumulative barplot
isotype_heavy_cluster_wide = isotype_heavy_cluster_table
isotype_heavy_cluster_wide = 100*isotype_heavy_cluster_wide / isotype_heavy_cluster_table_line_total
isotype_heavy_cluster_wide$isotype.heavy = row.names( isotype_heavy_cluster_wide)
isotype_heavy_cluster_wide = reshape2::melt( isotype_heavy_cluster_wide, id.vars = c( "isotype.heavy"))
print( ggplot( isotype_heavy_cluster_wide) +
         geom_bar( aes( x = isotype.heavy, fill = variable, y=value), stat="identity") + 
         scale_fill_manual( values = isotype_heavy_palette) +
         labs( fill = "Isotype", x= "Cluster", y = "Percentage") + 
         theme_minimal() +
         ggtitle( "Percentage of cell in cluster for each heavy chain isotype")
)

### plot the distribution of the Heavy chains for each tissue and each mouse
  
  # Get the count table of entries with same tissue, mouse, isotype.heavy
  tissue_mouse_isotype_df = as.data.frame( table( clonotype_seurat_df[, c( "tissue", "mouse", "isotype.heavy")]))
  # Reorder the isotypes
  tissue_mouse_isotype_df$isotype.heavy = factor( tissue_mouse_isotype_df$isotype.heavy, levels = c( "IGHD", "IGHM", "IGHG1", "IGHG2B", "IGHG2C", "IGHG3", "IGHA"))
  # Compute the percentage of cells in each isotype respect to the same tissue and mouse
  tissue_mouse_isotype_df$pct.tissue.mouse = apply( tissue_mouse_isotype_df, 1, function( row){
    return( 100 * as.numeric( row[ "Freq"])/sum( tissue_mouse_isotype_df[ tissue_mouse_isotype_df$tissue == row[ "tissue"] & tissue_mouse_isotype_df$mouse == row[ "mouse"], "Freq"]))
  })
  # Plot the distribution of percentage of the cells along the various isotypes
  print( ggplot( tissue_mouse_isotype_df) +
           stat_summary( aes(x=isotype.heavy, y=pct.tissue.mouse, col=tissue),
                         fun.y = "mean",
                         fun.ymin = function(x) mean(x) - sd(x),
                         fun.ymax = function(x) mean(x) + sd(x),
                         geom = "pointrange", group = 1, size = 1) +
           stat_summary( aes(x=isotype.heavy, y=pct.tissue.mouse, group = tissue, col=tissue), fun.y = "mean", geom = "line") +
           geom_jitter( aes(x=isotype.heavy, y=pct.tissue.mouse, col=tissue), width = 0.2) +
           xlab( "Isotype") + ylab( "Percentage of cells") + ggtitle( "Distribution of cells along isotypes and tissues") +
           theme_classic()
  )
  
  print( htmltools::tagList( datatable( as.matrix( tissue_mouse_isotype_df), caption="Distribution of cells along isotypes and tissues")))
  write.table( tissue_mouse_isotype_df,
               file = file.path( PATH_ANALYSIS_OUTPUT, paste0( SAMPLE_NAME, "_distribution_cells_HEAVYchain_isotypes_Alltissues.csv")),
               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
### plot the distribution of the Heavy chains for each cluster and each mouse
  
  # Get the count table of entries with same cluster, mouse, isotype.heavy
  cluster_mouse_isotype_df = as.data.frame( table( clonotype_seurat_df[ , c( "seurat_clusters", "mouse", "isotype.heavy")]))
  # Reorder the isotypes
  cluster_mouse_isotype_df$isotype.heavy = factor( cluster_mouse_isotype_df$isotype.heavy, levels = c( "IGHD", "IGHM", "IGHG1", "IGHG2B", "IGHG2C", "IGHG3", "IGHA"))
  # Compute the percentage of cells in each isotype respect to the same cluster and mouse
  cluster_mouse_isotype_df$pct.cluster.mouse = apply( cluster_mouse_isotype_df, 1, function( row){
    current_sum = sum( cluster_mouse_isotype_df[ cluster_mouse_isotype_df$seurat_clusters == row[ "seurat_clusters"] & cluster_mouse_isotype_df$mouse == row[ "mouse"], "Freq"])
    if( current_sum == 0){
      if( as.numeric( row[ "Freq"]) == 0){
        return( 0)
      }else{
        stop( "ERROR in table: sum is zero but nominator is not:", row)
      }
    }else{
      return( 100 * as.numeric( row[ "Freq"]) / current_sum)
    }
  })
  # Plot the distribution of percentage of the cells along the various isotypes
  print( ggplot( cluster_mouse_isotype_df) +
           stat_summary( aes(x=isotype.heavy, y=pct.cluster.mouse, col=seurat_clusters),
                         fun.y = "mean",
                         fun.ymin = function(x) mean(x) - sd(x),
                         fun.ymax = function(x) mean(x) + sd(x),
                         geom = "pointrange", group = 1, size = 1) +
           stat_summary( aes(x=isotype.heavy, y=pct.cluster.mouse, group = seurat_clusters, col=seurat_clusters), fun.y = "mean", geom = "line") +
           geom_jitter( aes(x=isotype.heavy, y=pct.cluster.mouse, col=seurat_clusters), width = 0.2) +
           xlab( "Isotype") + ylab( "Percentage of cells") + ggtitle( paste( "Distribution of cells along isotypes and clusters")) +
           theme_classic()
  )
  
  print( htmltools::tagList( datatable( as.matrix( cluster_mouse_isotype_df), caption=paste( "Distribution of cells along isotypes and clusters on tissue"))))
  write.table( cluster_mouse_isotype_df,
               file = file.path( PATH_ANALYSIS_OUTPUT, paste0( SAMPLE_NAME, "_distribution_cells_HEAVYchain_isotypes.csv")),
               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ................................................................................................
## ISOTYPE LIGHT-CHAIN versus CLUSTERS
# ................................................................................................

cat("<HR><H5>ISOTYPE LIGHT-CHAIN versus CLUSTERS</H5>")

# Show the dispersion of isotypes of light chain over clusters in a datatable
isotype_light_cluster_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "seurat_clusters", "isotype.light" )]))
isotype_light_cluster_table_line_total = apply( isotype_light_cluster_table, 1, sum)
isotype_light_cluster_table_column_total = apply( isotype_light_cluster_table, 2, sum)
isotype_light_cluster_table = isotype_light_cluster_table[ isotype_light_cluster_table_line_total > 0, isotype_light_cluster_table_column_total > 0]
print( htmltools::tagList( datatable( t(isotype_light_cluster_table), caption = "Dispersion of light chain isotypes over clusters (all cells)")))

# Show the dispersion of clusters over isotype of light chain in percentage in a cumulative barplot
isotype_light_cluster_wide = isotype_light_cluster_table
isotype_light_cluster_wide = 100*isotype_light_cluster_wide / isotype_light_cluster_table_line_total
isotype_light_cluster_wide$isotype.light = row.names( isotype_light_cluster_wide)
isotype_light_cluster_wide = reshape2::melt( isotype_light_cluster_wide, id.vars = c( "isotype.light"))
print( ggplot( isotype_light_cluster_wide) +
         geom_bar( aes( x = isotype.light, fill = variable, y=value), stat="identity") + 
         scale_fill_manual( values = isotype_light_palette) +
         labs( fill = "Isotype", x= "Cluster", y = "Percentage") + 
         theme_minimal() +
         ggtitle( "Percentage of cell in cluster for each light chain isotype")
)

### plot the distribution of the Light chains for each tissue and each mouse

  
  # Get the count table of entries with same tissue, mouse, isotype.heavy
  tissue_mouse_isotype_df = as.data.frame( table( clonotype_seurat_df[, c( "tissue", "mouse", "isotype.light")]))
  # Reorder the isotypes
  tissue_mouse_isotype_df$isotype.light = factor( tissue_mouse_isotype_df$isotype.light, levels = c( "IGKC", "IGLC1", "IGLC2", "IGLC3"))
  # Compute the percentage of cells in each isotype respect to the same tissue and mouse
  tissue_mouse_isotype_df$pct.tissue.mouse = apply( tissue_mouse_isotype_df, 1, function( row){
    return( 100 * as.numeric( row[ "Freq"])/sum( tissue_mouse_isotype_df[ tissue_mouse_isotype_df$tissue == row[ "tissue"] & tissue_mouse_isotype_df$mouse == row[ "mouse"], "Freq"]))
  })
  # Plot the distribution of percentage of the cells along the various isotypes
  print( ggplot( tissue_mouse_isotype_df) +
           stat_summary( aes(x=isotype.light, y=pct.tissue.mouse, col=tissue),
                         fun.y = "mean",
                         fun.ymin = function(x) mean(x) - sd(x),
                         fun.ymax = function(x) mean(x) + sd(x),
                         geom = "pointrange", group = 1, size = 1) +
           stat_summary( aes(x=isotype.light, y=pct.tissue.mouse, group = tissue, col=tissue), fun.y = "mean", geom = "line") +
           geom_jitter( aes(x=isotype.light, y=pct.tissue.mouse, col=tissue), width = 0.2) +
           xlab( "Isotype") + ylab( "Percentage of cells") + ggtitle( "Distribution of cells along isotypes and tissues") +
           theme_classic()
  )
  
  print( htmltools::tagList( datatable( as.matrix( tissue_mouse_isotype_df), caption="Distribution of cells along isotypes and tissues")))
  write.table( tissue_mouse_isotype_df,
               file = file.path( PATH_ANALYSIS_OUTPUT, paste0( SAMPLE_NAME, "_distribution_cells_LIGHTchain_isotypes_Alltissues.csv")),
               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
### If we are looking at data with a single tissue, plot the distribution of the Light chains for each cluster and each mouse
  
  # Get the count table of entries with same cluster, mouse, isotype.heavy
  cluster_mouse_isotype_df = as.data.frame( table( clonotype_seurat_df[ , c( "seurat_clusters", "mouse", "isotype.light")]))
  # Reorder the isotypes
  cluster_mouse_isotype_df$isotype.light = factor( cluster_mouse_isotype_df$isotype.light, levels = c( "IGKC", "IGLC1", "IGLC2", "IGLC3"))
  # Compute the percentage of cells in each isotype respect to the same cluster and mouse
  cluster_mouse_isotype_df$pct.cluster.mouse = apply( cluster_mouse_isotype_df, 1, function( row){
    current_sum = sum( cluster_mouse_isotype_df[ cluster_mouse_isotype_df$seurat_clusters == row[ "seurat_clusters"] & cluster_mouse_isotype_df$mouse == row[ "mouse"], "Freq"])
    if( current_sum == 0){
      if( as.numeric( row[ "Freq"]) == 0){
        return( 0)
      }else{
        stop( "ERROR in table: sum is zero but nominator is not:", row)
      }
    }else{
      return( 100 * as.numeric( row[ "Freq"]) / current_sum)
    }
  })
  # Plot the distribution of percentage of the cells along the various isotypes
  print( ggplot( cluster_mouse_isotype_df) +
           stat_summary( aes(x=isotype.light, y=pct.cluster.mouse, col=seurat_clusters),
                         fun.y = "mean",
                         fun.ymin = function(x) mean(x) - sd(x),
                         fun.ymax = function(x) mean(x) + sd(x),
                         geom = "pointrange", group = 1, size = 1) +
           stat_summary( aes(x=isotype.light, y=pct.cluster.mouse, group = seurat_clusters, col=seurat_clusters), fun.y = "mean", geom = "line") +
           geom_jitter( aes(x=isotype.light, y=pct.cluster.mouse, col=seurat_clusters), width = 0.2) +
           xlab( "Isotype") + ylab( "Percentage of cells") + ggtitle( paste( "Distribution of cells along isotypes and clusters")) +
           theme_classic()
  )
  
  print( htmltools::tagList( datatable( as.matrix( cluster_mouse_isotype_df), caption=paste( "Distribution of cells along isotypes and clusters on tissue"))))
  write.table( cluster_mouse_isotype_df,
               file = file.path( PATH_ANALYSIS_OUTPUT, paste0( SAMPLE_NAME, "_distribution_cells_LIGHTchain_isotypes.csv")),
               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# ................................................................................................
## V-SEGMENT HEAVY-CHAIN versus CLUSTERS
# ................................................................................................

cat("<HR><H5>V-SEGMENT HEAVY-CHAIN versus CLUSTERS</H5>")

# Show the dispersion of V-segment of heavy chain over clusters in a datatable
v.segment_heavy_cluster_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "seurat_clusters", "v.segment.heavy")]))
v.segment_heavy_cluster_table_line_total = apply( v.segment_heavy_cluster_table, 1, sum)
v.segment_heavy_cluster_table_column_total = apply( v.segment_heavy_cluster_table, 2, sum)
v.segment_heavy_cluster_table = v.segment_heavy_cluster_table[ v.segment_heavy_cluster_table_line_total > 0, v.segment_heavy_cluster_table_column_total > 0]
print( htmltools::tagList( datatable( t(v.segment_heavy_cluster_table), caption = "Dispersion of heavy chain V-Segments over clusters (all cells)")))

# Show the dispersion of clusters over V-segment of heavy chain in percentage in a cumulative barplot
v.segment_heavy_cluster_wide = v.segment_heavy_cluster_table
v.segment_heavy_cluster_wide = 100*v.segment_heavy_cluster_wide / v.segment_heavy_cluster_table_line_total
v.segment_heavy_cluster_wide$v.segment.heavy = row.names( v.segment_heavy_cluster_wide)
v.segment_heavy_cluster_wide = reshape2::melt( v.segment_heavy_cluster_wide, id.vars = c( "v.segment.heavy"))
print( ggplot( v.segment_heavy_cluster_wide) +
         geom_bar( aes( x = v.segment.heavy, fill = variable, y=value), stat="identity") + 
         scale_fill_manual( values = vsegment_heavy_palette) +
         labs( fill = "V Segment", x= "Cluster", y = "Percentage") + 
         theme_minimal() +
         ggtitle( "Percentage of cell in cluster for each heavy chain V-segment")
)

# ................................................................................................
## V-SEGMENT LIGHT-CHAIN versus CLUSTERS
# ................................................................................................

cat("<HR><H5>V-SEGMENT LIGHT-CHAIN versus CLUSTERS</H5>")

# Show the dispersion of V-segment of light chain over clusters in a datatable
v.segment_light_cluster_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "seurat_clusters", "v.segment.light")]))
v.segment_light_cluster_table_line_total = apply( v.segment_light_cluster_table, 1, sum)
v.segment_light_cluster_table_column_total = apply( v.segment_light_cluster_table, 2, sum)
v.segment_light_cluster_table = v.segment_light_cluster_table[ v.segment_light_cluster_table_line_total > 0, v.segment_light_cluster_table_column_total > 0]
print( htmltools::tagList( datatable( t(v.segment_light_cluster_table), caption = "Dispersion of light chain V-Segments over clusters (all cells)")))

# Show the dispersion of clusters over V-segment of light chain in percentage in a cumulative barplot
v.segment_light_cluster_wide = v.segment_light_cluster_table
v.segment_light_cluster_wide = 100*v.segment_light_cluster_wide / v.segment_light_cluster_table_line_total
v.segment_light_cluster_wide$v.segment.light = row.names( v.segment_light_cluster_wide)
v.segment_light_cluster_wide = reshape2::melt( v.segment_light_cluster_wide, id.vars = c( "v.segment.light"))
print( ggplot( v.segment_light_cluster_wide) +
          geom_bar( aes( x = v.segment.light, fill = variable, y=value), stat="identity") + 
          scale_fill_manual( values = vsegment_light_palette) +
          labs( fill = "V Segment", x= "Cluster", y = "Percentage") + 
          theme_minimal() +
          ggtitle( "Percentage of cell in cluster for each light chain V-segment")
)

# ................................................................................................
## HEAVY-CHAIN CDR3 LENGTH versus CLUSTERS
# ................................................................................................

cat("<HR><H5>HEAVY-CHAIN CDR3 LENGTH versus CLUSTERS</H5>")

# Show the dispersion of CDR3 Length of heavy chain over clusters in a datatable
cdr3.length_heavy_cluster_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "seurat_clusters", "cdr3.length.heavy")]))
cdr3.length_heavy_cluster_table_line_total = apply( cdr3.length_heavy_cluster_table, 1, sum)
cdr3.length_heavy_cluster_table_column_total = apply( cdr3.length_heavy_cluster_table, 2, sum)
cdr3.length_heavy_cluster_table = cdr3.length_heavy_cluster_table[ cdr3.length_heavy_cluster_table_line_total > 0, cdr3.length_heavy_cluster_table_column_total > 0]
print( htmltools::tagList( datatable( t( cdr3.length_heavy_cluster_table), caption = "Dispersion of CDR3 length over clusters (all cells)")))

# Show the dispersion of clusters CDR3 length of heavy chain in percentage in a cumulative barplot
cdr3.length_heavy_cluster_wide = cdr3.length_heavy_cluster_table
cdr3.length_heavy_cluster_wide = 100*cdr3.length_heavy_cluster_wide / cdr3.length_heavy_cluster_table_line_total
cdr3.length_heavy_cluster_wide$cluster = row.names( cdr3.length_heavy_cluster_wide)
cdr3.length_heavy_cluster_wide = reshape2::melt( cdr3.length_heavy_cluster_wide, id.vars = c( "cluster"))
print( ggplot( cdr3.length_heavy_cluster_wide) +
         geom_bar( aes( x = cluster, fill = variable, y=value), stat="identity") + 
         labs( fill = "CDR3 length", x= "Cluster", y = "Percentage") + 
         theme_minimal()  +
         ggtitle( "Percentage of cell in cluster for each heavy chain CDR3 length")
)

### plot the distribution of the Heavy chains for each tissue and each mouse

  
  # Get the count table of entries with same tissue, mouse, isotype.heavy
  tissue_mouse_cdr3length_df = as.data.frame( table( clonotype_seurat_df[, c( "tissue", "mouse", "cdr3.length.heavy")]))
  # Reorder the isotypes
  tissue_mouse_cdr3length_df$cdr3.length.heavy = as.numeric( as.character( tissue_mouse_cdr3length_df$cdr3.length.heavy))
  # Compute the percentage of cells in each isotype respect to the same tissue and mouse
  tissue_mouse_cdr3length_df$pct.tissue.mouse = apply( tissue_mouse_cdr3length_df, 1, function( row){
    return( 100 * as.numeric( row[ "Freq"])/sum( tissue_mouse_cdr3length_df[ tissue_mouse_cdr3length_df$tissue == row[ "tissue"] & tissue_mouse_cdr3length_df$mouse == row[ "mouse"], "Freq"]))
  })
  # Plot the distribution of percentage of the cells along the various isotypes
  print( ggplot( tissue_mouse_cdr3length_df) +
           stat_summary( aes(x=cdr3.length.heavy, y=pct.tissue.mouse, col=tissue),
                         fun.y = "mean",
                         fun.ymin = function(x) mean(x) - sd(x),
                         fun.ymax = function(x) mean(x) + sd(x),
                         geom = "pointrange", group = 1, size = 1) +
           stat_summary( aes(x=cdr3.length.heavy, y=pct.tissue.mouse, group = tissue, col=tissue), fun.y = "mean", geom = "line") +
           geom_jitter( aes(x=cdr3.length.heavy, y=pct.tissue.mouse, col=tissue), width = 0.2) +
           xlab( "CDR3 length") + ylab( "Percentage of cells") + ggtitle( "Distribution of cells along CDR3 length and tissues") +
           theme_classic()
  )
  
  # Print dataframe as datatable
  print( htmltools::tagList( datatable( as.matrix( tissue_mouse_cdr3length_df), caption="Distribution of cells along CDR3 length and tissues")))
  
  # Export dataframe to file
  cdr3.length.heavy.file_path=file.path( PATH_ANALYSIS_OUTPUT, paste0( SAMPLE_NAME, "_CDR3_length_heavy_stats_mouse_alltissues.csv"))
  write.table( tissue_mouse_cdr3length_df, file = cdr3.length.heavy.file_path,
               sep=",", col.names=TRUE, row.names=FALSE, quote = FALSE)
  
### plot the distribution of the Heavy chains for each cluster and each mouse  
  
  # Get the count table of entries with same cluster, mouse, isotype.heavy
  cluster_mouse_cdr3length_df = as.data.frame( table( clonotype_seurat_df[ , c( "seurat_clusters", "mouse", "cdr3.length.heavy")]))
  # Reorder the isotypes
  cluster_mouse_cdr3length_df$cdr3.length.heavy = as.numeric( as.character( cluster_mouse_cdr3length_df$cdr3.length.heavy))
  # Compute the percentage of cells in each isotype respect to the same cluster and mouse
  cluster_mouse_cdr3length_df$pct.cluster.mouse = apply( cluster_mouse_cdr3length_df, 1, function( row){
    current_sum = sum( cluster_mouse_cdr3length_df[ cluster_mouse_cdr3length_df$seurat_clusters == row[ "seurat_clusters"] & cluster_mouse_cdr3length_df$mouse == row[ "mouse"], "Freq"])
    if( current_sum == 0){
      if( as.numeric( row[ "Freq"]) == 0){
        return( 0)
      }else{
        stop( "ERROR in table: sum is zero but nominator is not:", row)
      }
    }else{
      return( 100 * as.numeric( row[ "Freq"]) / current_sum)
    }
  })
  # Plot the distribution of percentage of the cells along the various isotypes
  print( ggplot( cluster_mouse_cdr3length_df) +
           stat_summary( aes(x=cdr3.length.heavy, y=pct.cluster.mouse, col=seurat_clusters),
                         fun.y = "mean",
                         fun.ymin = function(x) mean(x) - sd(x),
                         fun.ymax = function(x) mean(x) + sd(x),
                         geom = "pointrange", group = 1, size = 1) +
           stat_summary( aes(x=cdr3.length.heavy, y=pct.cluster.mouse, group = seurat_clusters, col=seurat_clusters), fun.y = "mean", geom = "line") +
           geom_jitter( aes(x=cdr3.length.heavy, y=pct.cluster.mouse, col=seurat_clusters), width = 0.2) +
           xlab( "Isotype") + ylab( "Percentage of cells") + ggtitle( paste( "Distribution of cells along CDR3 length")) +
           theme_classic()
  )
  
  # Print dataframe as datatable
  print( htmltools::tagList( datatable( as.matrix( cluster_mouse_cdr3length_df), caption=paste( "Distribution of cells along CDR3 length on tissue"))))
  
  # Export dataframe to file
  cdr3.length.heavy.file_path=file.path( PATH_ANALYSIS_OUTPUT, paste0( SAMPLE_NAME, "_CDR3_length_heavy_stats_mouse.csv"))
  write.table( cluster_mouse_cdr3length_df, file = cdr3.length.heavy.file_path,
               sep=",", col.names=TRUE, row.names=FALSE, quote = FALSE)

# ................................................................................................
## LIGHT-CHAIN CDR3 LENGTH versus CLUSTERS
# ................................................................................................

cat("<HR><H5>LIGHT-CHAIN CDR3 LENGTH versus CLUSTERS</H5>")

# Show the dispersion of CDR3 Length of light chain over clusters in a datatable
cdr3.length_light_cluster_table = as.data.frame.matrix( table( clonotype_seurat_df[, c( "seurat_clusters", "cdr3.length.light")]))
cdr3.length_light_cluster_table_line_total = apply( cdr3.length_light_cluster_table, 1, sum)
cdr3.length_light_cluster_table_column_total = apply( cdr3.length_light_cluster_table, 2, sum)
cdr3.length_light_cluster_table = cdr3.length_light_cluster_table[ cdr3.length_light_cluster_table_line_total > 0, cdr3.length_light_cluster_table_column_total > 0]
print( htmltools::tagList( datatable( t(cdr3.length_light_cluster_table), caption = "Dispersion of CDR3 length over clusters (all cells)")))

# Show the dispersion of clusters over CDR3 length of light chain in percentage in a cumulative barplot
cdr3.length_light_cluster_wide = cdr3.length_light_cluster_table
cdr3.length_light_cluster_wide = 100*cdr3.length_light_cluster_wide / cdr3.length_light_cluster_table_line_total
cdr3.length_light_cluster_wide$cluster = row.names( cdr3.length_light_cluster_wide)
cdr3.length_light_cluster_wide = reshape2::melt( cdr3.length_light_cluster_wide, id.vars = c( "cluster"))
cdr3.length_light_cluster_wide$variable = factor( cdr3.length_light_cluster_wide$variable, levels = sort( levels( cdr3.length_light_cluster_wide$variable)))
print( ggplot( cdr3.length_light_cluster_wide) +
         geom_bar( aes( x = cluster, fill = variable, y=value), stat="identity") + 
         labs( fill = "CDR3 length", x= "Cluster", y = "Percentage") + 
         theme_minimal()  +
         ggtitle( "Percentage of cell in cluster for each light chain CDR3 length")
)

# ................................................................................................
## HEAVY-CHAIN NUMBER OF MUTATION versus CLUSTERS
# ................................................................................................

cat("<HR><H5>HEAVY-CHAIN NUMBER OF MUTATION versus CLUSTERS</H5>")

# Compose groups of mices to analyze : all mice together and then one mouse at a time
mice_groups = list( as.character( levels( clonotype_seurat_df$mouse)))
for( mice in levels( clonotype_seurat_df$mouse)){
  mice_groups = c( mice_groups, mice)
}

# Compute the percentage table at display the plots
for( current_mices in mice_groups){

  # Show the dispersion of number of mutation of heavy chain over clusters in a datatable
  mutation.qual_heavy_cluster_df = clonotype_seurat_df[ which( clonotype_seurat_df$mouse %in% current_mices), c( "seurat_clusters", "mutation.qual.heavy")]
  mutation.qual_heavy_cluster_df$mutation.qual.heavy = as.numeric( as.character( mutation.qual_heavy_cluster_df$mutation.qual.heavy))
  mutation.qual_heavy_cluster_table = as.data.frame.matrix( table( mutation.qual_heavy_cluster_df))
  mutation.qual_heavy_cluster_table_line_total = apply( mutation.qual_heavy_cluster_table, 1, sum)
  mutation.qual_heavy_cluster_table_column_total = apply( mutation.qual_heavy_cluster_table, 2, sum)
  mutation.qual_heavy_cluster_table = mutation.qual_heavy_cluster_table[ mutation.qual_heavy_cluster_table_line_total > 0, mutation.qual_heavy_cluster_table_column_total > 0]
  print( htmltools::tagList( datatable( t( mutation.qual_heavy_cluster_table), caption = "Dispersion of mutation number over clusters (all cells)")))
  
  # Show the dispersion of clusters CDR3 length of heavy chain in percentage in a cumulative barplot
  mutation.qual_heavy_cluster_wide = mutation.qual_heavy_cluster_table
  mutation.qual_heavy_cluster_wide = 100*mutation.qual_heavy_cluster_wide / mutation.qual_heavy_cluster_table_line_total
  mutation.qual_heavy_cluster_wide$cluster = row.names( mutation.qual_heavy_cluster_wide)
  mutation.qual_heavy_cluster_wide = reshape2::melt( mutation.qual_heavy_cluster_wide, id.vars = c( "cluster"))
  mutation.qual_heavy_cluster_wide$variable = as.numeric( as.character( mutation.qual_heavy_cluster_wide$variable))
  # print( ggplot( mutation.qual_heavy_cluster_wide) +
  #          geom_bar( aes( x = cluster, fill = variable, y=value), stat="identity") + 
  #          labs( fill = "# mutation", x= "Cluster", y = "Percentage") + 
  #          theme_minimal()  +
  #          ggtitle( paste( "Percentage of cell in cluster for each heavy chain number of mutation", 
  #                     "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_heavy_cluster_df),")"))
  # )
  
  print( ggplot( mutation.qual_heavy_cluster_wide) +
           geom_bar( aes( x = variable, fill = cluster, y=value), stat="identity", position="dodge") + 
           labs( fill = "Cluster", x= "# mutation", y = "Percentage") + 
           theme_minimal()  +
           theme(axis.text.x = element_text(angle = 45)) +
           ggtitle( paste( "Percentage of cell in cluster for each heavy chain number of mutation", 
                      "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_heavy_cluster_df),")"))
  )
  
  print( ggplot( mutation.qual_heavy_cluster_wide) +
           geom_bar( aes( x = variable, fill = cluster, y=value), stat="identity", position="dodge") + 
           labs( fill = "Cluster", x= "# mutation", y = "Percentage") + 
           facet_wrap( . ~ cluster) +
           theme_minimal()  +
           theme(axis.text.x = element_text(angle = 45)) +
           ggtitle( paste( "Percentage of cell in cluster for each heavy chain number of mutation", 
                           "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_heavy_cluster_df),")"))
  )

}

# ................................................................................................
cat("<HR><H5>HEAVY-CHAIN NUMBER OF MUTATION versus CLUSTERS between CLONOTYPE</H5>")
# ................................................................................................

mutation.qual_heavy_cluster_clonotype_df = clonotype_seurat_df[ , c( "seurat_clusters", "mutation.qual.heavy", "large.clonotype.similarity.symbol", "mouse")]
mutation.qual_heavy_cluster_clonotype_df$mutation.qual.heavy = as.numeric( as.character( mutation.qual_heavy_cluster_clonotype_df$mutation.qual.heavy))
mutation.qual_heavy_cluster_clonotype_df = mutation.qual_heavy_cluster_clonotype_df[ -which( is.na( mutation.qual_heavy_cluster_clonotype_df$large.clonotype.similarity.symbol)), ]
print( ggplot( mutation.qual_heavy_cluster_clonotype_df) +
         geom_boxplot( aes( x=seurat_clusters, y=mutation.qual.heavy, fill= large.clonotype.similarity.symbol)) +
         geom_jitter( aes( x=seurat_clusters, y=mutation.qual.heavy, col = mouse), width = 0.2) + 
         labs( fill = "Clonotype", x= "Cluster", y = "# mutation", col = "Mouse") + 
         facet_wrap( . ~ large.clonotype.similarity.symbol) + 
         theme_classic()
       
)


# ................................................................................................
## LIGHT-CHAIN NUMBER OF MUTATION versus CLUSTERS
# ................................................................................................

cat("<HR><H5>LIGHT-CHAIN NUMBER OF MUTATION versus CLUSTERS</H5>")

# Compose groups of mices to analyze : all mice together and then one mouse at a time
mice_groups = list( as.character( levels( clonotype_seurat_df$mouse)))
for( mice in levels( clonotype_seurat_df$mouse)){
  mice_groups = c( mice_groups, mice)
}

# Compute the percentage table at display the plots
for( current_mices in mice_groups){
  
  # Show the dispersion of number of mutation of light chain over clusters in a datatable
  mutation.qual_light_cluster_df = clonotype_seurat_df[ which( clonotype_seurat_df$mouse %in% current_mices), c( "seurat_clusters", "mutation.qual.light")]
  mutation.qual_light_cluster_df$mutation.qual.light = as.numeric( as.character( mutation.qual_light_cluster_df$mutation.qual.light))
  mutation.qual_light_cluster_table = as.data.frame.matrix( table( mutation.qual_light_cluster_df))
  mutation.qual_light_cluster_table_line_total = apply( mutation.qual_light_cluster_table, 1, sum)
  mutation.qual_light_cluster_table_column_total = apply( mutation.qual_light_cluster_table, 2, sum)
  mutation.qual_light_cluster_table = mutation.qual_light_cluster_table[ mutation.qual_light_cluster_table_line_total > 0, mutation.qual_light_cluster_table_column_total > 0]
  print( htmltools::tagList( datatable( t( mutation.qual_light_cluster_table), caption = "Dispersion of mutation number over clusters (all cells)")))
  
  # Show the dispersion of clusters CDR3 length of light chain in percentage in a cumulative barplot
  mutation.qual_light_cluster_wide = mutation.qual_light_cluster_table
  mutation.qual_light_cluster_wide = 100*mutation.qual_light_cluster_wide / mutation.qual_light_cluster_table_line_total
  mutation.qual_light_cluster_wide$cluster = row.names( mutation.qual_light_cluster_wide)
  mutation.qual_light_cluster_wide = reshape2::melt( mutation.qual_light_cluster_wide, id.vars = c( "cluster"))
  mutation.qual_light_cluster_wide$variable = as.numeric( as.character( mutation.qual_light_cluster_wide$variable))
  # print( ggplot( mutation.qual_light_cluster_wide) +
  #          geom_bar( aes( x = cluster, fill = variable, y=value), stat="identity") + 
  #          labs( fill = "# mutation", x= "Cluster", y = "Percentage") + 
  #          theme_minimal()  +
  #          ggtitle( paste( "Percentage of cell in cluster for each light chain number of mutation", 
  #                     "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_light_cluster_df),")"))
  # )
  
  print( ggplot( mutation.qual_light_cluster_wide) +
           geom_bar( aes( x = variable, fill = cluster, y=value), stat="identity", position="dodge") + 
           labs( fill = "Cluster", x= "# mutation", y = "Percentage") + 
           theme_minimal()  +
           theme(axis.text.x = element_text(angle = 45)) +
           ggtitle( paste( "Percentage of cell in cluster for each light chain number of mutation", 
                           "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_light_cluster_df),")"))
  )
  
  print( ggplot( mutation.qual_light_cluster_wide) +
           geom_bar( aes( x = variable, fill = cluster, y=value), stat="identity", position="dodge") + 
           labs( fill = "Cluster", x= "# mutation", y = "Percentage") + 
           facet_wrap( . ~ cluster) +
           theme_minimal()  +
           theme(axis.text.x = element_text(angle = 45)) +
           ggtitle( paste( "Percentage of cell in cluster for each light chain number of mutation", 
                           "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_light_cluster_df),")"))
  )
  
}

# ................................................................................................
cat("<HR><H5>LIGHT-CHAIN NUMBER OF MUTATION versus CLUSTERS between CLONOTYPE</H5>")
# ................................................................................................

mutation.qual_light_cluster_clonotype_df = clonotype_seurat_df[ , c( "seurat_clusters", "mutation.qual.light", "large.clonotype.similarity.symbol", "mouse")]
mutation.qual_light_cluster_clonotype_df$mutation.qual.light = as.numeric( as.character( mutation.qual_light_cluster_clonotype_df$mutation.qual.light))
mutation.qual_light_cluster_clonotype_df = mutation.qual_light_cluster_clonotype_df[ -which( is.na( mutation.qual_light_cluster_clonotype_df$large.clonotype.similarity.symbol)), ]
print( ggplot( mutation.qual_light_cluster_clonotype_df) +
         geom_boxplot( aes( x=seurat_clusters, y=mutation.qual.light, fill= large.clonotype.similarity.symbol)) +
         geom_jitter( aes( x=seurat_clusters, y=mutation.qual.light, col = mouse), width = 0.2) + 
         labs( fill = "Clonotype", x= "Cluster", y = "# mutation", col = "Mouse") + 
         facet_wrap( . ~ large.clonotype.similarity.symbol) + 
         theme_classic()
       
)
