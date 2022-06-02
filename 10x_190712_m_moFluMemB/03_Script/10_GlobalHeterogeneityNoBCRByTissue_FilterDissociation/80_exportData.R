# ########################################################################################
# This script aims to export data to file in various format to be used
# in next step analysis
# ########################################################################################


## @knitr export_data

# Export the Seurat object as RDS
seurat_rds_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "sc10x.rna.seurat_", CURRENT_TISSUE, ".RDS"))
cat("<br>Exportgin seurat object to RDS:", seurat_rds_file)
saveRDS( sc10x.rna.seurat, file = seurat_rds_file)

# Export marker genes for each cluster
all_marker_genes_table_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "allMarkers_", CURRENT_TISSUE, ".tab"))
cat("<br>Exporting all marker genes table :", all_marker_genes_table_file)
write.table( markers, file = all_marker_genes_table_file,
             sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Export marker genes for each cluster
marker_genes_table_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "topMarkersDT_", CURRENT_TISSUE, ".tab"))
cat("<br>Exporting marker genes table :", marker_genes_table_file)
write.table( topMarkersDT, file = marker_genes_table_file,
             sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Export the t-SNE embedding
tsne_embedding_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "tSNE_Embedding_", CURRENT_TISSUE, ".tab"))
cat("<br>Exporting t-SNE embedding:", tsne_embedding_file)
write.table( sc10x.rna.seurat@reductions$tsne@cell.embeddings, file = tsne_embedding_file,
             sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Export the UMAP embedding
umap_embedding_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "UMAP_Embedding_", CURRENT_TISSUE, ".tab"))
cat("<br>Exporting UMAP embedding:", umap_embedding_file)
write.table( sc10x.rna.seurat@reductions$umap@cell.embeddings, file = umap_embedding_file,
             sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Export the cluster assignation
cluster_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "clusters_", CURRENT_TISSUE, ".tab"))
cat("<br>Export cluster assignation file:", cluster_file)
write.table( data.frame( cell.id = row.names( sc10x.rna.seurat@meta.data),
                         cluster = sc10x.rna.seurat@meta.data$seurat_clusters,
                         sample = sc10x.rna.seurat@meta.data$MULTI_Sample),
             file = cluster_file ,
             sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
