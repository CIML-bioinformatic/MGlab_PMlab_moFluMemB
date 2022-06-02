# ########################################################
# This script aims to prepare figures for the final
# article
# ########################################################


## @knitr prepare_article_figure_UMAPplot

# Prepare the data
Idents( sc10x.rna.seurat) = "seurat_clusters"
dimReducData = cbind(as.data.frame(Embeddings(sc10x.rna.seurat, reduction = "umap")),
                     Cluster = as.numeric(as.character(Idents( sc10x.rna.seurat))));

# Save UMAP with cluster to SVG file and tsv file for publication
ggplot( dimReducData) + 
  geom_point( aes( x = UMAP_1, y = UMAP_2, color = as.factor( Cluster))) + 
  scale_color_manual( values = hue_pal()( nlevels( as.factor( dimReducData$Cluster)))) + 
  labs(title = paste( CURRENT_TISSUE, "mem B cells"), x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 8))) + 
  theme_classic()

ggsave( file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "UMAP_Embedding_WithCluster_", CURRENT_TISSUE, ".svg")), width = 10, height=10)
write.table( dimReducData, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "UMAP_Embedding_WithCluster_", CURRENT_TISSUE, ".tsv")),
             row.names = TRUE, col.names = NA, sep="\t", quote = FALSE)


# Save UMAP feature plots of some interesting genes
for( current_gene in c(  "Ccr6", "Fcer2a", "Plac8", "Rgs1", "Cxcr3", "Cd55", "Ms4a6", "Ms4a6c", "Ms4a6b", "Ms4a6d", "Itga4", "Ms4a6c", "Cr2", "Itgax")){
  if( current_gene %in% rownames( sc10x.rna.seurat)){
    print( FeaturePlot( sc10x.rna.seurat, features = current_gene))
    ggsave( file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "UMAP_FeaturePlot_", current_gene, "_", CURRENT_TISSUE, ".svg")), width = 10, height=10)
    cat("<BR>")
  }else{
    cat("<BR>Missing feature in Seurat object :", current_gene)
  }
}


