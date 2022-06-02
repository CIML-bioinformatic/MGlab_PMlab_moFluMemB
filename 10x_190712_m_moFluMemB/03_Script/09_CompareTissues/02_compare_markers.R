# ##################################################################################
# This script aims to compare the marker genes between clusters of different
# tissue
# ##################################################################################


## @knitr compare_markers_between_clusters


# Parse over the pairs of tissues
for( index_tissue_1 in 1: (length( names( PATH_ALL_MARKER_GENES_TABLE_FILE))-1)){
  for( index_tissue_2 in (index_tissue_1+1): length( names( PATH_ALL_MARKER_GENES_TABLE_FILE))){
  
    # Get the names of the tissues
    current_tissue_1 = names( PATH_ALL_MARKER_GENES_TABLE_FILE)[ index_tissue_1]
    current_tissue_2 = names( PATH_ALL_MARKER_GENES_TABLE_FILE)[ index_tissue_2]
    cat("<H4>Comparison of markers genes between clusters of", current_tissue_1, " and", current_tissue_2, "</H4>")
    
    print( ggplot( merge_umap_cluster_list[[ current_tissue_1]]) + 
             geom_point( aes( x=UMAP_1, y=UMAP_2, col=cluster)) +
             ggtitle( paste( "UMAP enbedding of", current_tissue_1)) +
             theme_classic()
    )
    print( ggplot( merge_umap_cluster_list[[ current_tissue_2]]) + 
             geom_point( aes( x=UMAP_1, y=UMAP_2, col=cluster)) +
             ggtitle( paste( "UMAP enbedding of", current_tissue_2)) +
             theme_classic()
    )
    
    # Get the table of markers genes of each tissues
    all_markers_tissue1_df = all_marker_table_list[[ current_tissue_1]]
    all_markers_tissue2_df = all_marker_table_list[[ current_tissue_2]]
    
    # Compute the size of intersection between list of genes of pairs of different clusters of the tissues
    # and compute the Jaccard index for each pair of clusters. Store both information in two different dataframes
    jaccard_index_df = data.frame()
    intersection_df = data.frame()
    for( cluster_tissue_1 in levels( all_markers_tissue1_df$cluster)){
      jaccard_index_set = vector()
      details_set = vector()

      for( cluster_tissue_2 in levels( all_markers_tissue2_df$cluster)){
        markers_tissue_1 = all_markers_tissue1_df[ which( all_markers_tissue1_df$cluster == cluster_tissue_1), "gene"]
        markers_tissue_2 = all_markers_tissue2_df[ which( all_markers_tissue2_df$cluster == cluster_tissue_2), "gene"]
        jaccard_index = length( intersect( markers_tissue_1, markers_tissue_2)) / length( unique( c( markers_tissue_1, markers_tissue_2)))
        jaccard_index_set = append( jaccard_index_set, jaccard_index)
        details_set = append( details_set, length( intersect( markers_tissue_1, markers_tissue_2)))
      }
      jaccard_index_df = rbind( jaccard_index_df, jaccard_index_set)
      intersection_df = rbind( intersection_df, details_set)
    }
    names( jaccard_index_df) = paste0( current_tissue_2, levels( all_markers_tissue2_df$cluster))
    row.names( jaccard_index_df) = paste0( current_tissue_1, levels( all_markers_tissue1_df$cluster))
    
    names( intersection_df) = paste0( current_tissue_2, levels( all_markers_tissue2_df$cluster))
    row.names( intersection_df) = paste0( current_tissue_1, levels( all_markers_tissue1_df$cluster))
    
    # Display the data.frame of jaccard indexes as a heatmap
    break_list = seq(0, 1, by = 0.1)
    pheatmap::pheatmap( 
                  jaccard_index_df,
                  color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length( break_list)),
                  cellwidth =70, cellheight =70,
                  fontsize = 18,
                  show_rownames = T, show_colnames = T,
                  angle_col = 45,
                  breaks = break_list,
                  display_numbers = matrix( paste0( as.matrix( intersection_df), " (", 100*signif( as.matrix( jaccard_index_df), 2), "%)"), ncol = length( levels( all_markers_tissue2_df$cluster))),
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  treeheight_row = 0, treeheight_col=0
      )
    
    # For each cluster of tissue1, search for the best match cluster (best jaccard index) among the cluster of the tissue2
    # then report the list of markers genes in the intersection of the cluster markers
    for( tissue1_index in 1:nrow( jaccard_index_df)){
      best_index = 0
      tissue1_cluster_name = row.names( jaccard_index_df)[ tissue1_index]
      
      for( tissue2_index in 1:ncol( jaccard_index_df)){
        tissue2_cluster_name = names( jaccard_index_df)[ tissue2_index]
        
        jaccard_index = jaccard_index_df[ tissue1_index, tissue2_index]
        
        if( jaccard_index > best_index){
          best_index = jaccard_index
          best_cluster_name = tissue2_cluster_name
          best_cluster = tissue2_index - 1
          markers_tissue_1 = all_markers_tissue1_df[ which( all_markers_tissue1_df$cluster == (tissue1_index-1)), "gene"]
          best_markers_tissue_2 = all_markers_tissue2_df[ which( all_markers_tissue2_df$cluster == (tissue2_index-1)), "gene"]
          best_markers_intersection = intersect( markers_tissue_1, best_markers_tissue_2)
        }
      }
      cat("<br><br><b>Best match for cluster ", tissue1_cluster_name, " (", length( markers_tissue_1), "genes) = ",
                                            best_cluster_name, " (", length( best_markers_tissue_2), "genes)", "</b>", sep="")
      cat("<br>Genes in intersection (", length( best_markers_intersection), " genes) : <br>")
      cat( pander( best_markers_intersection))
    }
    
    # For each cluster of tissue2, search fo rthe best match cluster (best jaccard index) among the cluster of the tissue1
    # then report the list of markers genes in the intersection of the cluster markers
    for( tissue2_index in 1:ncol( jaccard_index_df)){
      best_index = 0
      tissue2_cluster_name = names( jaccard_index_df)[ tissue2_index]
      
      for( tissue1_index in 1:nrow( jaccard_index_df)){
        tissue1_cluster_name = row.names( jaccard_index_df)[ tissue1_index]
        
        jaccard_index = jaccard_index_df[ tissue1_index, tissue2_index]
        
        if( jaccard_index > best_index){
          best_index = jaccard_index
          best_cluster_name = tissue1_cluster_name
          best_cluster = tissue1_index - 1
          markers_tissue_2 = all_markers_tissue2_df[ which( all_markers_tissue2_df$cluster == (tissue2_index-1)), "gene"]
          best_markers_tissue_1 = all_markers_tissue1_df[ which( all_markers_tissue1_df$cluster == (tissue1_index-1)), "gene"]
          best_markers_intersection = intersect( markers_tissue_2, best_markers_tissue_1)
        }
      }
      cat("<br><br><b>Best match for cluster ", tissue2_cluster_name, " (", length( markers_tissue_2), "genes) = ",
                                            best_cluster_name, " (", length( best_markers_tissue_1), "genes)", "</b>", sep="")
      cat("<br>Genes in intersection (", length( best_markers_intersection), " genes) :<br>")
      cat( pander( best_markers_intersection))
    }
    
    cat("<br><br>")
  }
}




