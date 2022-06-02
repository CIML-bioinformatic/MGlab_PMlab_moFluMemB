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


# Plot to file the heatmap of marker gene lists similarities

for( tissue1 in names( markers_heatmap_jaccard_similarity_list)){
  for( tissue2 in names( markers_heatmap_jaccard_similarity_list[[ tissue1]])){
    svg( file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Heatmap_MarkerGenes_JaccardComparison_", tissue1, "_vs_", tissue2, ".svg")))
      print( markers_heatmap_jaccard_similarity_list[[ tissue1]][[ tissue2]])
    dev.off()
  }
}

for( tissue1 in names( markers_heatmap_szymkiewicz_simpson_similarity_list)){
  for( tissue2 in names( markers_heatmap_szymkiewicz_simpson_similarity_list[[ tissue1]])){
    svg( file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Heatmap_MarkerGenes_Szymkiewicz-SimpsonComparison_", tissue1, "_vs_", tissue2, ".svg")))
    print( markers_heatmap_szymkiewicz_simpson_similarity_list[[ tissue1]][[ tissue2]])
    dev.off()
  }
}

ggsave( file = file.path( PATH_ANALYSIS_OUTPUT, "Heatmap_MarkerGenes_Jaccard_Comparison_AllTissues.svg"),
          plot = all_tissues_all_marker_jaccard_index_heatmap)

ggsave( file = file.path( PATH_ANALYSIS_OUTPUT, "Heatmap_MarkerGenes_Szymkiewicz-Simpson_Comparison_AllTissues.svg"),
        plot = all_tissues_all_marker_szymkiewicz_simpson_index_heatmap)


# ##################################################################################################
# Is there clonal overlap in memory B cells from different clusters in the different organs ?
#  
# For that analysis you can use the first dataset 190712, and potentially the second 191105 (but less clones detected). 
# The idea is to visualize and quantify the clonotypes which are found in multiple organ.subsets. One possible visualization are 
# circos plots as used in the powerpoint slide 1. An alternative visualization could be heatmaps backed up with pie charts 
# quantifying the numbers of clones in different categories, as shown in powerpoint slide 2.

# ##################################################################################################

# First, try a representation of the data as a Chord Diagram
# help : https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html#chord-diagram-colors
# .............................................................................................


# List the cells in the filtered set that have also clonotype info
filtered_cells_with_clonotype_info = intersect( merge_umap_cluster_alltissue$cell.id, row.names( clonotype_and_seurat_data_df))

# Prepare a dataframe with only the filtered cells having Clonotype info, merging the data from the two sources
chord_diagram_df = clonotype_and_seurat_data_df
chord_diagram_df$cell.id = row.names( chord_diagram_df)
chord_diagram_df = merge( chord_diagram_df[ filtered_cells_with_clonotype_info, c( "cell.id", "full.clonotype.similarity.symbol", "mouse") ], 
                          merge_umap_cluster_alltissue[ which( merge_umap_cluster_alltissue$cell.id %in% filtered_cells_with_clonotype_info), c( "cell.id", "cluster", "sample", "tissue" )],
                          by = "cell.id")

# Add a tissue.cluster column
chord_diagram_df$tissue.cluster = paste( chord_diagram_df$tissue,
                                         chord_diagram_df$cluster, sep = ".")

clono_count = table( chord_diagram_df$full.clonotype.similarity.symbol)
kept_clono = names( clono_count)[ which( clono_count >= 1)]
big_clono_chord_diagram_df = chord_diagram_df[ which( chord_diagram_df$full.clonotype.similarity.symbol %in% kept_clono),]
big_clono_chord_diagram_table_df = data.frame( table( as.character( big_clono_chord_diagram_df$full.clonotype.similarity.symbol), as.character( big_clono_chord_diagram_df$tissue.cluster)))

# Plot the complete clonotype versus tissue.cluster heatmap
ggplot( big_clono_chord_diagram_table_df) +
  geom_tile( aes( x= Var1, y = Var2, fill=Freq)) + theme( axis.text.x = element_text( angle = 90, size = 2))

# Parse the mouse to create plot for each of them
for( current_mouse in levels( chord_diagram_df$mouse)){
  
  cat("<H5>Analysis for mouse", current_mouse, "</H5>")
  
  # Select the data of the mouse
  mouse_chord_diagram_df = chord_diagram_df[ which( chord_diagram_df$mouse == current_mouse), ]

  
  # Create the dataframe required for the heatmaps and the ChordDiagram
  circle_plot_df = data.frame( stringsAsFactors = FALSE)
  interaction_df = data.frame( )
  overlap_df = data.frame( )
  
  # Define the list of Tissue.Cluster to parse
  tissu_cluster_set = sort( unique( mouse_chord_diagram_df$tissue.cluster))
  
  # Parse the pairs of Tissue.cluster to look for their overlap/interaction in clonotypes
  for( current_tissue_cluster_index_1 in 1:(length(tissu_cluster_set))){
    
    # -- Get the first tissue.cluster to look at
    current_tissue_cluster_1 = tissu_cluster_set[ current_tissue_cluster_index_1]
    
    # -- Initialize some variables used to build the heatmaps dataframes
    current_overlap_set = vector()
    current_interaction_set = vector()
    
    # -- Inner parse the Tissue.Cluster list
    for( current_tissue_cluster_index_2 in 1:length(tissu_cluster_set)){
      
      # -- Get the second tissue.cluster to look at
      current_tissue_cluster_2 = tissu_cluster_set[ current_tissue_cluster_index_2]
      
      # -- Get the useful data from global dataframe
      current_df_1 = mouse_chord_diagram_df[ mouse_chord_diagram_df$tissue.cluster == current_tissue_cluster_1, ]
      current_df_2 = mouse_chord_diagram_df[ mouse_chord_diagram_df$tissue.cluster == current_tissue_cluster_2, ]
      
      # -- Get the list of Clonotype in common to the two Tissue.Clusters
      common_clonotypes_number = length( intersect( current_df_1$full.clonotype.similarity.symbol, current_df_2$full.clonotype.similarity.symbol))
      
      # -- If we did not compare the two Tissue.Clusters already, compare them
      # -- Add the information for the ChordDiagram dataframe
      if( current_tissue_cluster_1 > current_tissue_cluster_2){   
        circle_plot_df = rbind( circle_plot_df, data.frame( from = current_tissue_cluster_1, to = current_tissue_cluster_2,
                                                          value = common_clonotypes_number, stringsAsFactors = FALSE))
      }

      # -- If we did not compare the two Tissue.Clusters already, compare them
      # if( current_tissue_cluster_1 > current_tissue_cluster_2){   
        # -- Accumulate the information for the interaction heatmap dataframe
        current_interaction_set = append( current_interaction_set, common_clonotypes_number)
        
        # -- Accumulate the information for the overlap heatmap dataframe
        current_overlap_set = append( current_overlap_set, common_clonotypes_number/min( length( unique( current_df_1$full.clonotype.similarity.symbol)),
                                                                                         length( unique( current_df_2$full.clonotype.similarity.symbol))))
      # }else{
      #  # If we did compare already the same tissue.clusters, just add NA to avoid replication of results due to symmetry
      #  current_overlap_set = append( current_overlap_set, NA)
      #  current_interaction_set = append( current_interaction_set, NA)
      # }
                              
    }
    
    # Accumulate the information to build the heatmap dataframes
    overlap_df = rbind( overlap_df, current_overlap_set)
    interaction_df = rbind( interaction_df, current_interaction_set)
  }
  
  # Rename the columns and row of the heatmap dataframes
  names( overlap_df) = tissu_cluster_set
  row.names( overlap_df) = tissu_cluster_set
  names(interaction_df) = tissu_cluster_set
  row.names( interaction_df) = tissu_cluster_set
  
  # Declare the tissue of origin of each tissue.cluster as Group
  group_df = mouse_chord_diagram_df[ !duplicated( mouse_chord_diagram_df$tissue.cluster), c( "tissue", "tissue.cluster")]
  group = structure( as.character( group_df$tissue)[ order( group_df$tissue.cluster)], names = group_df$tissue.cluster[ order( group_df$tissue.cluster)])
  
  # Provide a color for each concerned tissue.cluster
  grid.colors = c(LN.0 = "chartreuse2",  LN.1 = "chartreuse3", LN.2 = "chartreuse4",
                  SP.0 = "firebrick", SP.1 = "firebrick1", SP.2 = "firebrick2", SP.3 = "firebrick3", SP.4 = "firebrick4",
                  LG.0 = "blue", LG.1 = "blue2", LG.2 = "blue4")
  
  # Plot the Chord Diagram
  svg( file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "ChordDiagram_ClonotypeInteraction_Count_Mouse", current_mouse, ".svg")) )
    circlize::chordDiagram( circle_plot_df, 
                            group = group,  
                            annotationTrack = c( "grid", "axis"),
                            annotationTrackHeight = c(0.08, 0.1),
                            grid.col = grid.colors, 
                            preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(circle_plot_df))))))
      title( paste( "Mouse", current_mouse), cex = 0.6)
    
      # Add the tissue.cluster labels
      circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
        circlize::circos.text(CELL_META$xcenter, CELL_META$ylim[1]+0.5, CELL_META$sector.index, 
                    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
        }, bg.border = NA)
  dev.off()
  
  # Second, try a representation of the data as a heatmap
  # help : https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html#chord-diagram-colors
  # .............................................................................................
  
  
  # INTERACTION TABLE
  # Compute clutering on interaction count table
  interaction_hclust = hclust( dist( interaction_df))
  interaction_clustered_df = interaction_df[ interaction_hclust$order, interaction_hclust$order]
  
  # Set to NA all values in the lower triangle part of the matrix
  interaction_clustered_df[upper.tri(interaction_clustered_df)] = NA
  interaction_clustered_df$cluster.tissue = factor( row.names( interaction_clustered_df), levels = row.names( interaction_clustered_df))
  
  # Convert matrix from wide format to long formt to use in ggplot
  melted_interaction_clustered_df <- reshape2::melt( interaction_clustered_df, id = "cluster.tissue", na.rm = TRUE)
  
  # Use ggplot to plot the heatmap (only triangle) with and adapted scale color
  melted_interaction_clustered_heatmap = ggplot( data = melted_interaction_clustered_df, aes( cluster.tissue, variable,fill = value)) +
    geom_tile( color = "white") +
    scale_fill_gradientn( colours = c( "blue", "white", "red"),
                          colors = c( "blue", "white", "red"),
                          values = c( 0, 
                                      0.5,
                                      1),
                          space = "Lab",
                          name="Clonotype count") +
    scale_y_discrete(position = "right") +
    theme_classic2() + 
    theme(axis.text.x = element_text( angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.text.y = element_text( vjust = 1, size = 12, hjust = 1),
          axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = c( 0.3, 0.8)) +
    coord_fixed()
  
  # Save the heatmap to SVG file
  ggsave( melted_interaction_clustered_heatmap, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Heatmap_ClonotypeInteraction_Count_Mouse", current_mouse, ".svg")) )
  write.table( melted_interaction_clustered_df, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Heatmap_ClonotypeInteraction_Count_Mouse", current_mouse, ".tsv")),
               sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  print( melted_interaction_clustered_heatmap)
  
  # OVERLAP TABLE
  # Compute clustering on Overlap index table
  overlap_hclust = hclust( dist( overlap_df))
  overlap_clustered_df = overlap_df[ overlap_hclust$order, overlap_hclust$order]
  
  # Set to NA all values in the lower triangle part of the matrix
  overlap_clustered_df[upper.tri(overlap_clustered_df)] = NA
  overlap_clustered_df$cluster.tissue = factor( row.names( overlap_clustered_df), levels = row.names( overlap_clustered_df))
  
  # Convert matrix from wide format to long formt to use in ggplot
  melted_overlap_clustered_df <- reshape2::melt( overlap_clustered_df, id = "cluster.tissue", na.rm = TRUE)
  
  # Use ggplot to plot the heatmap (only triangle) with and adapted scale color
  melted_overlap_clustered_heatmap = ggplot( data = melted_overlap_clustered_df, aes( cluster.tissue, variable,fill = value)) +
    geom_tile( color = "white") +
    scale_fill_gradientn( colours = c( "blue", "white", "red", "grey"),
                          colors = c( "blue", "white", "red", "grey"),
                          values = c( 0, 
                                      0.25,
                                      0.9999,
                                      1),
                          space = "Lab",
                          name="Overlap index") +
    scale_y_discrete(position = "right") +
    theme_classic2() + 
    theme(axis.text.x = element_text( angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.text.y = element_text( vjust = 1, size = 12, hjust = 1),
          axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = c( 0.3, 0.8)) +
    coord_fixed()
  
  # Save the heatmap to SVG file and tab file
  ggsave( melted_overlap_clustered_heatmap, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Heatmap_ClonotypeInteraction_OverlapIndex_Mouse", current_mouse, ".svg")) )
  write.table( melted_overlap_clustered_df, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "Heatmap_ClonotypeInteraction_OverlapIndex_Mouse", current_mouse, ".tsv")),
               sep ="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  print( melted_overlap_clustered_heatmap)
}


