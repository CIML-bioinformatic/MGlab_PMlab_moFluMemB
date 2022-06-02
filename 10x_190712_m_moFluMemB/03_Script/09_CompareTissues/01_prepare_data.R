# ##################################################################################
# This script aims to get the data from the tissue specific analysis
# ##################################################################################


## @knitr load_data

# Load all the table of maker gene for each tissue
all_marker_table_list = list()
for( current_tissue in names( PATH_ALL_MARKER_GENES_TABLE_FILE)){

  all_marker_table_list[[ current_tissue]] = read.table( file = PATH_ALL_MARKER_GENES_TABLE_FILE[[ current_tissue]], header = TRUE, sep="\t", quote = NULL)
  all_marker_table_list[[ current_tissue]]$cluster = as.factor( all_marker_table_list[[ current_tissue]]$cluster)
}

# Load the UMAP embeddings for each tissue
umap_embedding_list = list()
for( current_tissue in names( PATH_UMAP_EMBEDDING_FILE)){
  
  umap_embedding_list[[ current_tissue]] = read.table( file = PATH_UMAP_EMBEDDING_FILE[[ current_tissue]], header = TRUE, sep="\t", quote = NULL)
}

# Load the cluster mapping for each tissue
cluster_mapping_list = list()
for( current_tissue in names( PATH_CLUSTER_MAPPING_FILE)){
  
  cluster_mapping_list[[ current_tissue]] = read.table( file = PATH_CLUSTER_MAPPING_FILE[[ current_tissue]], header = TRUE, sep="\t", quote = NULL)
  cluster_mapping_list[[ current_tissue]]$cluster = as.factor( cluster_mapping_list[[ current_tissue]]$cluster)
}

# Merge the umap information with the cluster information for each tissue and plot the UMAP embedding with cluster colors
merge_umap_cluster_list = list()
for( current_tissue in names( PATH_CLUSTER_MAPPING_FILE)){
  
  merge_umap_cluster_list[[ current_tissue]] = cluster_mapping_list[[ current_tissue]]
  merge_umap_cluster_list[[ current_tissue]]$UMAP_1 = umap_embedding_list[[ current_tissue]][ merge_umap_cluster_list[[ current_tissue]]$cell.id, "UMAP_1"]
  merge_umap_cluster_list[[ current_tissue]]$UMAP_2 = umap_embedding_list[[ current_tissue]][ merge_umap_cluster_list[[ current_tissue]]$cell.id, "UMAP_2"]
  print( ggplot( merge_umap_cluster_list[[ current_tissue]]) + 
    geom_point( aes( x=UMAP_1, y=UMAP_2, col=cluster)) +
    ggtitle( paste( "UMAP enbedding of", current_tissue)) +
    theme_classic()
  )
}

# Load the matrices of expression of the most variables genes for each tissue
variable_genes_expression_list = list()
for( current_tissue in names( PATH_VARIABLE_GENES_EXPRESSION_FILE)){
  variable_genes_expression_list[[ current_tissue]] = read.table( file = PATH_VARIABLE_GENES_EXPRESSION_FILE[[current_tissue]], 
                                                                  header=TRUE, row.names = 1, sep="\t")
}




