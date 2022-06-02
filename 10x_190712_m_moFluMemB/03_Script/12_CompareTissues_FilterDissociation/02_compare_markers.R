# ##################################################################################
# This script aims to compare the marker genes between clusters of different
# tissue
# ##################################################################################


## @knitr compare_markers_between_clusters


# Parse over the pairs of tissues
markers_heatmap_jaccard_similarity_list = list()
markers_heatmap_szymkiewicz_simpson_similarity_list = list()
for( index_tissue_1 in 1: (length( names( PATH_ALL_MARKER_GENES_TABLE_FILE))-1)){
  
  # Get the names of the tissue 1
  current_tissue_1 = names( PATH_ALL_MARKER_GENES_TABLE_FILE)[ index_tissue_1]
  
  # Initialize the list of heatmap figures
  markers_heatmap_jaccard_similarity_list[[ current_tissue_1]] = list()
  markers_heatmap_szymkiewicz_simpson_similarity_list[[ current_tissue_1]] = list()
  
  for( index_tissue_2 in (index_tissue_1+1): length( names( PATH_ALL_MARKER_GENES_TABLE_FILE))){
  
    # Get the names of the tissue 2
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
    jaccard_index_df = data.frame( stringsAsFactors = FALSE)
    szymkiewicz_simpson_index_df = data.frame( stringsAsFactors = FALSE)
    intersection_df = data.frame()
    for( cluster_tissue_1 in levels( all_markers_tissue1_df$cluster)){
      jaccard_index_set = vector()
      szymkiewicz_simpson_index_set = vector()
      details_set = vector()

      for( cluster_tissue_2 in levels( all_markers_tissue2_df$cluster)){
        markers_tissue_1 = all_markers_tissue1_df[ which( all_markers_tissue1_df$cluster == cluster_tissue_1), "gene"]
        markers_tissue_2 = all_markers_tissue2_df[ which( all_markers_tissue2_df$cluster == cluster_tissue_2), "gene"]
        jaccard_index = length( intersect( markers_tissue_1, markers_tissue_2)) / length( unique( c( markers_tissue_1, markers_tissue_2)))
        szymkiewicz_simpson_index = length( intersect( markers_tissue_1, markers_tissue_2)) / min( length( markers_tissue_1), length( markers_tissue_2))
        jaccard_index_set = append( jaccard_index_set, jaccard_index)
        szymkiewicz_simpson_index_set = append( szymkiewicz_simpson_index_set, szymkiewicz_simpson_index)
        details_set = append( details_set, length( intersect( markers_tissue_1, markers_tissue_2)))
      }
      jaccard_index_df = rbind( jaccard_index_df, jaccard_index_set)
      szymkiewicz_simpson_index_df = rbind( szymkiewicz_simpson_index_df, szymkiewicz_simpson_index_set)
      intersection_df = rbind( intersection_df, details_set)
    }
    names( jaccard_index_df) = paste0( current_tissue_2, levels( all_markers_tissue2_df$cluster))
    row.names( jaccard_index_df) = paste0( current_tissue_1, levels( all_markers_tissue1_df$cluster))
    
    names( szymkiewicz_simpson_index_df) = paste0( current_tissue_2, levels( all_markers_tissue2_df$cluster))
    row.names( szymkiewicz_simpson_index_df) = paste0( current_tissue_1, levels( all_markers_tissue1_df$cluster))
    
    names( intersection_df) = paste0( current_tissue_2, levels( all_markers_tissue2_df$cluster))
    row.names( intersection_df) = paste0( current_tissue_1, levels( all_markers_tissue1_df$cluster))
    
    # Display the data.frame of jaccard indexes as a heatmap and create the same heatmap without numbers in a list
    # for SVG export at the end of the report
    break_list = seq(0, 1, by = 0.1)
    cat("<BR><H5>Jaccard Index Similarity</H5>")
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
    
    markers_heatmap_jaccard_similarity_list[[ current_tissue_1]][[ current_tissue_2]] = pheatmap::pheatmap( 
      jaccard_index_df,
      color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length( break_list)),
      cellwidth =70, cellheight =70,
      fontsize = 18,
      show_rownames = T, show_colnames = T,
      angle_col = 45,
      breaks = break_list,
      # display_numbers = matrix( paste0( as.matrix( intersection_df), " (", 100*signif( as.matrix( jaccard_index_df), 2), "%)"), ncol = length( levels( all_markers_tissue2_df$cluster))),
      cluster_rows = FALSE, cluster_cols = FALSE,
      treeheight_row = 0, treeheight_col=0
      )
    
    # Display the data.frame of Szymkiewicz-Simpson indexes as a heatmap and create the same heatmap without numbers in a list
    # for SVG export at the end of the report    
    cat("<BR><H5>Szymkiewicz-Simpson Index Similarity</H5>")
    pheatmap::pheatmap( 
      szymkiewicz_simpson_index_df,
      color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length( break_list)),
      cellwidth =70, cellheight =70,
      fontsize = 18,
      show_rownames = T, show_colnames = T,
      angle_col = 45,
      breaks = break_list,
      display_numbers = matrix( paste0( as.matrix( intersection_df), " (", 100*signif( as.matrix( szymkiewicz_simpson_index_df), 2), "%)"), ncol = length( levels( all_markers_tissue2_df$cluster))),
      cluster_rows = FALSE, cluster_cols = FALSE,
      treeheight_row = 0, treeheight_col=0
    )
    
    markers_heatmap_szymkiewicz_simpson_similarity_list[[ current_tissue_1]][[ current_tissue_2]] = pheatmap::pheatmap( 
      szymkiewicz_simpson_index_df,
      color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length( break_list)),
      cellwidth =70, cellheight =70,
      fontsize = 18,
      show_rownames = T, show_colnames = T,
      angle_col = 45,
      breaks = break_list,
      # display_numbers = matrix( paste0( as.matrix( intersection_df), " (", 100*signif( as.matrix( szymkiewicz_simpson_index_df), 2), "%)"), ncol = length( levels( all_markers_tissue2_df$cluster))),
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


# Compare all clusters of all tissues together
# .............................................

cat("<H4>Comparison of markers genes between clusters of all tissues</H4>")

# Combine all marker datafrae fo all tissues, creating a new column to paste tissue name and cluster number
all_tissues_all_marker_table_df = data.frame()
for( current_tissue in names( all_marker_table_list)){
  current_df = all_marker_table_list[[ current_tissue]]
  current_df$tissue.cluster = paste0( current_tissue, "_", current_df$cluster)
  all_tissues_all_marker_table_df = rbind( all_tissues_all_marker_table_df, current_df)
}

all_tissues_all_marker_table_df = all_tissues_all_marker_table_df[ order( all_tissues_all_marker_table_df$p_val_adj, decreasing = FALSE),]

# Initialize result datafranme
jaccard_index_df = data.frame( stringsAsFactors = FALSE )
szymkiewicz_simpson_index_df = data.frame(stringsAsFactors = FALSE )
intersection_df = data.frame( stringsAsFactors = FALSE)

# determine the list of tissue per cluster combinaisons
tissue_cluster_set = sort( unique( all_tissues_all_marker_table_df$tissue.cluster))

# loop over pairs of tissue+cluster to compare their list of markers
for( cluster_tissue_1 in tissue_cluster_set){
  
  # Initialize inner loop result sets
  jaccard_index_set = vector()
  szymkiewicz_simpson_index_set = vector()
  details_set = vector()
  
  for( cluster_tissue_2 in tissue_cluster_set){
    
    # Determine the list of maker genes for each tissue+cluster
    markers_tissue_1 = all_tissues_all_marker_table_df[ which( all_tissues_all_marker_table_df$tissue.cluster == cluster_tissue_1), "gene"]
    markers_tissue_2 = all_tissues_all_marker_table_df[ which( all_tissues_all_marker_table_df$tissue.cluster == cluster_tissue_2), "gene"]
    
    markers_tissue_1 = markers_tissue_1[ 1:min( 100, length( markers_tissue_1))]
    markers_tissue_2 = markers_tissue_2[ 1:min( 100, length( markers_tissue_2))]
    
    # Compute the similarity indexes for th elist of marker genes
    jaccard_index = length( intersect( markers_tissue_1, markers_tissue_2)) / length( unique( c( markers_tissue_1, markers_tissue_2)))
    szymkiewicz_simpson_index = length( intersect( markers_tissue_1, markers_tissue_2)) / min( length( markers_tissue_1), length( markers_tissue_2))
    
    # Combine the results of inner loop
    jaccard_index_set = append( jaccard_index_set, jaccard_index)
    szymkiewicz_simpson_index_set = append( szymkiewicz_simpson_index_set, szymkiewicz_simpson_index)
    details_set = append( details_set, length( intersect( markers_tissue_1, markers_tissue_2)))
  }
  
  # Concatenate the outer loop result to a single dataframe for each index
  jaccard_index_df = rbind( jaccard_index_df, jaccard_index_set)
  szymkiewicz_simpson_index_df = rbind( szymkiewicz_simpson_index_df, szymkiewicz_simpson_index_set)
  intersection_df = rbind( intersection_df, details_set)
}

# Add the names of columns and row for each resulting dataframe
names( jaccard_index_df) = tissue_cluster_set
row.names( jaccard_index_df) = tissue_cluster_set

names( szymkiewicz_simpson_index_df) = tissue_cluster_set
row.names( szymkiewicz_simpson_index_df) = tissue_cluster_set

names( intersection_df) = tissue_cluster_set
row.names( intersection_df) = tissue_cluster_set

# Plot the heatmaps of clusters genes markers similarities
# For the Jaccard index
# .........................................................

cat("<H5>Heatmaps with Jaccard Index</H5>")

cat("<H6>Comparing", paste( names( PATH_ALL_MARKER_GENES_TABLE_FILE), collapse =" and "), "<H6>")

# Convert the Jaccard similarities dataframe to matrix 
jaccard_index_mat = as.matrix( jaccard_index_df)

# Computing heatmap for all tissues
# **************************************

# Plot the direct heatmap
print( pheatmap::pheatmap( 
  jaccard_index_mat,
  color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length( break_list)),
  na_col="white",
  cellwidth =30, cellheight =30,
  fontsize = 18,
  show_rownames = T, show_colnames = T,
  angle_col = 45,
  fontsize_number = 8,
  cluster_rows = TRUE, cluster_cols = TRUE,
  treeheight_row = 0, treeheight_col=0
))

# Cluster the column of matrix with the Jaccard distance
hc = hclust(dist(1 - jaccard_index_mat))
plot( hc)

# Reorder the matrix with the computed dendrogram
jaccard_index_mat = jaccard_index_mat[hc$order, hc$order]

# Set to NA all values in the lower triangle part of the matrix
jaccard_index_mat[lower.tri(jaccard_index_mat)] = NA

# Convert matrix from wide format to long formt to use in ggplot
melted_cormat <- reshape2::melt( jaccard_index_mat, na.rm = TRUE)

# Use ggplot to plot the heatmap (only triangle) with and adapted scale color
all_tissues_all_marker_jaccard_index_heatmap = ggplot( data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile( color = "white") +
  scale_fill_gradientn( colours = c( "blue", "white", "red", "grey"),
                        colors = c( "blue", "white", "red", "grey"),
                        values = c( 0, 
                                    0.5,
                                    0.9999,
                                    1),
                        space = "Lab",
                        name="Jaccard Similarity") +
  scale_y_discrete(position = "right") +
  theme_classic2() + 
  theme(axis.text.x = element_text( angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text( vjust = 1, size = 12, hjust = 1),
        axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = c( 0.3, 0.8)) +
  coord_fixed()

print( all_tissues_all_marker_jaccard_index_heatmap)

# Computing heatmap removing one tissue
# **************************************
for( current_excluded_tissue in names( PATH_ALL_MARKER_GENES_TABLE_FILE)){

  cat("<H6>Comparing", paste( names( PATH_ALL_MARKER_GENES_TABLE_FILE)[ ! names( PATH_ALL_MARKER_GENES_TABLE_FILE) == current_excluded_tissue], collapse =" and "), "<H6>")
  
  # Convert the Jaccard similarities dataframe to matrix 
  jaccard_index_mat = as.matrix( jaccard_index_df)
  jaccard_index_mat = jaccard_index_mat[ -which( grepl( current_excluded_tissue, rownames( jaccard_index_mat))),  -which( grepl( current_excluded_tissue, colnames( jaccard_index_mat)))]
  
  # Plot the direct heatmap
  print( pheatmap::pheatmap( 
    jaccard_index_mat,
    color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length( break_list)),
    na_col="white",
    cellwidth =30, cellheight =30,
    fontsize = 18,
    show_rownames = T, show_colnames = T,
    angle_col = 45,
    fontsize_number = 8,
    cluster_rows = TRUE, cluster_cols = TRUE,
    treeheight_row = 0, treeheight_col=0
  ))
  
  # Cluster the column of matrix with the Jaccard distance
  hc = hclust(dist(1 - jaccard_index_mat))
  plot( hc)
  
  # Reorder the matrix with the computed dendrogram
  jaccard_index_mat = jaccard_index_mat[hc$order, hc$order]
  
  # Set to NA all values in the lower triangle part of the matrix
  jaccard_index_mat[lower.tri(jaccard_index_mat)] = NA
  
  # Convert matrix from wide format to long formt to use in ggplot
  melted_cormat <- reshape2::melt( jaccard_index_mat, na.rm = TRUE)
  
  # Use ggplot to plot the heatmap (only triangle) with and adapted scale color
  all_tissuesm1_all_marker_jaccard_index_heatmap = ggplot( data = melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile( color = "white") +
    scale_fill_gradientn( colours = c( "blue", "white", "red", "grey"),
                          colors = c( "blue", "white", "red", "grey"),
                          values = c( 0, 
                                      0.5,
                                      0.9999,
                                      1),
                          space = "Lab",
                          name="Jaccard Similarity") +
    scale_y_discrete(position = "right") +
    theme_classic2() + 
    theme(axis.text.x = element_text( angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.text.y = element_text( vjust = 1, size = 12, hjust = 1),
          axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = c( 0.3, 0.8)) +
    coord_fixed()
  
  print( all_tissuesm1_all_marker_jaccard_index_heatmap)

}

# Plot the heatmaps of clusters genes markers similarities
# For the Szymkiewicz-Simpson index
# .........................................................

cat("<H5>Heatmap with Szymkiewicz-Simpson Index</H5>")

# Computing heatmap for all tissues
# **************************************
cat("<H6>Comparing", paste( names( PATH_ALL_MARKER_GENES_TABLE_FILE), collapse =" and "), "<H6>")

# Convert the Szymkiewicz-Simpson similarities dataframe to matrix 
szymkiewicz_simpson_index_mat = as.matrix( szymkiewicz_simpson_index_df)

# Plot the direct heatmap
print( pheatmap::pheatmap( 
  szymkiewicz_simpson_index_mat,
  color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length( break_list)),
  na_col="white",
  cellwidth =30, cellheight =30,
  fontsize = 18,
  show_rownames = T, show_colnames = T,
  angle_col = 45,
  # breaks = break_list,
  # display_numbers = matrix( paste0( as.matrix( intersection_df), "\n", 100*signif( as.matrix( szymkiewicz_simpson_index_df), 2), "%"), ncol = length( tissue_cluster_set)),
  fontsize_number = 8,
  cluster_rows = TRUE, cluster_cols = TRUE,
  treeheight_row = 0, treeheight_col=0
))

# Cluster the column of matrix with the Jaccard distances
hc = hclust( dist(1 - szymkiewicz_simpson_index_mat))
plot( hc)

# Reorder the matrix with the computed dendrogram
szymkiewicz_simpson_index_mat = szymkiewicz_simpson_index_mat[hc$order, hc$order]

# Set to NA all values in the lower triangle part of the matrix
szymkiewicz_simpson_index_mat[lower.tri(szymkiewicz_simpson_index_mat)] = NA

# Convert matrix from wide format to long formt to use in ggplot
melted_cormat <- reshape2::melt( szymkiewicz_simpson_index_mat, na.rm = TRUE)

# Compute the max value that is not 1 (self cluster similarity)
max_non_one_value = max (melted_cormat$value[ -which( melted_cormat$value == 1)])

# Use ggplot to plot the heatmap (only triangle) with and adapted scale color
all_tissues_all_marker_szymkiewicz_simpson_index_heatmap = ggplot( data = melted_cormat, aes(Var2, Var1, fill = value)) +
  geom_tile( color = "white") +
  scale_fill_gradientn( colours = c( "blue", "white", "red", "grey"),
                        colors = c( "blue", "white", "red", "grey"),
                        values = c( 0, 
                                    0.5,
                                    0.9999,
                                    1),
                        space = "Lab",
                        name="Szymkiewicz-Simpson Similarity") +  scale_y_discrete(position = "right") +
  theme_classic2() + 
  theme(axis.text.x = element_text( angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text( vjust = 1, size = 12, hjust = 1),
        axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = c( 0.3, 0.8)) +
  coord_fixed()

print( all_tissues_all_marker_szymkiewicz_simpson_index_heatmap)


# Computing heatmap removing one tissue
# **************************************
for( current_excluded_tissue in names( PATH_ALL_MARKER_GENES_TABLE_FILE)){

  
  cat("<H6>Comparing", paste( names( PATH_ALL_MARKER_GENES_TABLE_FILE)[ ! names( PATH_ALL_MARKER_GENES_TABLE_FILE) == current_excluded_tissue], collapse =" and "), "<H6>")
  
  # Convert the Szymkiewicz-Simpson similarities dataframe to matrix 
  szymkiewicz_simpson_index_mat = as.matrix( szymkiewicz_simpson_index_df)
  
  szymkiewicz_simpson_index_mat = szymkiewicz_simpson_index_mat[ -which( grepl( current_excluded_tissue, rownames( szymkiewicz_simpson_index_mat))),  -which( grepl( current_excluded_tissue, colnames( szymkiewicz_simpson_index_mat)))]
  
  # Plot the direct heatmap
  print( pheatmap::pheatmap( 
    szymkiewicz_simpson_index_mat,
    color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length( break_list)),
    na_col="white",
    cellwidth =30, cellheight =30,
    fontsize = 18,
    show_rownames = T, show_colnames = T,
    angle_col = 45,
    # breaks = break_list,
    # display_numbers = matrix( paste0( as.matrix( intersection_df), "\n", 100*signif( as.matrix( szymkiewicz_simpson_index_df), 2), "%"), ncol = length( tissue_cluster_set)),
    fontsize_number = 8,
    cluster_rows = TRUE, cluster_cols = TRUE,
    treeheight_row = 0, treeheight_col=0
  ))
  
  # Cluster the column of matrix with the Jaccard distances
  hc = hclust( dist(1 - szymkiewicz_simpson_index_mat))
  plot( hc)
  
  # Reorder the matrix with the computed dendrogram
  szymkiewicz_simpson_index_mat = szymkiewicz_simpson_index_mat[hc$order, hc$order]
  
  # Set to NA all values in the lower triangle part of the matrix
  szymkiewicz_simpson_index_mat[lower.tri(szymkiewicz_simpson_index_mat)] = NA
  
  # Convert matrix from wide format to long formt to use in ggplot
  melted_cormat <- reshape2::melt( szymkiewicz_simpson_index_mat, na.rm = TRUE)
  
  # Compute the max value that is not 1 (self cluster similarity)
  max_non_one_value = max (melted_cormat$value[ -which( melted_cormat$value == 1)])
  
  # Use ggplot to plot the heatmap (only triangle) with and adapted scale color
  all_tissuesm1_all_marker_szymkiewicz_simpson_index_heatmap = ggplot( data = melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile( color = "white") +
    scale_fill_gradientn( colours = c( "blue", "white", "red", "grey"),
                          colors = c( "blue", "white", "red", "grey"),
                          values = c( 0, 
                                      0.5,
                                      0.9999,
                                      1),
                          space = "Lab",
                          name="Szymkiewicz-Simpson Similarity") +  scale_y_discrete(position = "right") +
    theme_classic2() + 
    theme(axis.text.x = element_text( angle = 45, vjust = 1, size = 12, hjust = 1),
          axis.text.y = element_text( vjust = 1, size = 12, hjust = 1),
          axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
          legend.position = c( 0.3, 0.8)) +
    coord_fixed()
  
  print( all_tissuesm1_all_marker_szymkiewicz_simpson_index_heatmap)
}