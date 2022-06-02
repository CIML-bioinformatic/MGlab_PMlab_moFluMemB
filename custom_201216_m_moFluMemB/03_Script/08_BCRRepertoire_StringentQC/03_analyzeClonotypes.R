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
## ALL CLONOTYPES versus CLUSTERS
# ................................................................................................

cat("<HR><H5>ALL CLONOTYPES versus CLUSTERS</H5>")

# Compute the number of clonotypes shared by clusters. Compute also the overlap index
intersect_df = data.frame()
overlap_df = data.frame()
for( current_cluster_1 in levels( clonotype_seurat_df$seurat_clusters)){
  clonotype_cluster_1 = clonotype_seurat_df[ which( clonotype_seurat_df$seurat_clusters == current_cluster_1), "full.clonotype.similarity.symbol"]
  intersect_set = vector()
  overlap_set = vector()
  for( current_cluster_2 in levels( clonotype_seurat_df$seurat_clusters)){  
    if( current_cluster_1 != current_cluster_2){
      clonotype_cluster_2 = clonotype_seurat_df[ which( clonotype_seurat_df$seurat_clusters == current_cluster_2), "full.clonotype.similarity.symbol"]
      intersect_set = append( intersect_set, length( intersect( clonotype_cluster_1, clonotype_cluster_2)))
      overlap_set = append( overlap_set, length( intersect( clonotype_cluster_1, clonotype_cluster_2))/min( length( clonotype_cluster_1), length( clonotype_cluster_2)))
    }else{
      intersect_set = append( intersect_set, NA)
      overlap_set = append( overlap_set, NA)
    }
  }

  intersect_df = rbind( intersect_df, intersect_set)
  overlap_df = rbind( overlap_df, overlap_set)
}

# Plot the heatmap of clonotypes counts
names( intersect_df) = as.character( levels( clonotype_seurat_df$seurat_clusters))
intersect_df[ upper.tri(intersect_df)] = NA
intersect_df$Cluster = as.character( levels( clonotype_seurat_df$seurat_clusters))
intersect_melt_df = reshape2::melt( intersect_df, id=c( 'Cluster'))

ggplot( data = intersect_melt_df, aes(x=Cluster, y=variable, fill = value)) +
geom_tile( color = "white") +
scale_fill_gradientn( colours = c( "blue", "white", "red"),
                      space = "Lab",
                      name="Clonotype count\nbetween clusters",
                      na.value = "white") +
scale_y_discrete(position = "right") +
theme_classic2() + 
theme(axis.text.x = element_text( angle = 45, vjust = 1, size = 12, hjust = 1),
      axis.text.y = element_text( vjust = 1, size = 12, hjust = 1),
      axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = c( 0.2, 0.8)) +
coord_fixed()
  
ggsave( filename = file.path( PATH_ANALYSIS_OUTPUT, "SharedClonotypeAcrossClusters_Heatmap_Counts.svg"))

# Plot the heatmap of clonotypes overlap index
names( overlap_df) = as.character( levels( clonotype_seurat_df$seurat_clusters))
overlap_df[ upper.tri(overlap_df)] = NA
overlap_df$Cluster = as.character( levels( clonotype_seurat_df$seurat_clusters))
overlap_melt_df = reshape2::melt( overlap_df, id=c( 'Cluster'))

ggplot( data = overlap_melt_df, aes(x=Cluster, y=variable, fill = value)) +
  geom_tile( color = "white") +
  scale_fill_gradientn( colours = c( "blue", "white", "red"),
                        space = "Lab",
                        name="Clonotype overlap index\nbetween clusters",
                        na.value = "white") +
  scale_y_discrete(position = "right") +
  theme_classic2() + 
  theme(axis.text.x = element_text( angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text( vjust = 1, size = 12, hjust = 1),
        axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = c( 0.2, 0.8)) +
  coord_fixed()
ggsave( filename = file.path( PATH_ANALYSIS_OUTPUT, "SharedClonotypeAcrossClusters_Heatmap_OverlapIndex.svg"))


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
large_clonotype_cluster_wide$cluster = row.names( large_clonotype_cluster_wide)
large_clonotype_cluster_wide = reshape2::melt( large_clonotype_cluster_wide, id.vars = c( "cluster"))
print( ggplot( large_clonotype_cluster_wide) +
         geom_bar( aes( x = cluster, fill = variable, y=value), stat="identity") + 
         scale_fill_manual( values = large_clonotype_palette) +
         labs( fill = "Clonotype", x= "Cluster", y = "Percentage") + 
         theme_minimal() + 
         ggtitle( "Percentage of cell in cluster for each large clonotype")
)


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
isotype_heavy_cluster_wide$cluster = row.names( isotype_heavy_cluster_wide)
isotype_heavy_cluster_wide = reshape2::melt( isotype_heavy_cluster_wide, id.vars = c( "cluster"))
print( ggplot( isotype_heavy_cluster_wide) +
         geom_bar( aes( x = cluster, fill = variable, y=value), stat="identity") + 
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
