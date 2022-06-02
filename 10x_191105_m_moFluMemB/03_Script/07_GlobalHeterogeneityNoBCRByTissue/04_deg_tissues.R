# ############################################################################################
# This script aims to look at response
# of cells to gene signature from previous analysis
# ############################################################################################

# ...........................................
# Get gene signature from previous analysis
# ...........................................

## @knitr marker_genes_from_other_dataset

# Marker genes from the first project experiment
CLUSTER_MARKERS_FOR_SIGNATURE_FILE = paste0( "/mnt/NAS7/Collaboration/moFluMemB/03_BIOINFO_ANALYSIS/10x_190712_m_moFluMemB/05_Output/07_GlobalHeterogeneityNoBCRByTissue/allMarkers_", CURRENT_TISSUE, "_wilcox.tab")

# Load the data
cluster_marker_signature_df = read.table( CLUSTER_MARKERS_FOR_SIGNATURE_FILE, sep="\t", header = TRUE)
cluster_marker_signature_df$cluster = as.factor( cluster_marker_signature_df$cluster)

# Show the marker genes of the previous experiment on the current data
DefaultAssay( sc10x.rna.seurat) = "RNA"
all_idents = vector()
for( cluster_id in levels( cluster_marker_signature_df$cluster)){
  
  gene_for_signature = list( c( as.character( cluster_marker_signature_df[ which( cluster_marker_signature_df$cluster == cluster_id), "gene"])))
  ident_name = paste0( "Markers_Cluster_", as.character( cluster_id), "_")
  sc10x.rna.seurat = AddModuleScore( sc10x.rna.seurat, gene_for_signature, name = ident_name)
  all_idents = c( all_idents, paste0( ident_name, "1"))
  
}

FeaturePlot(object = sc10x.rna.seurat, features = all_idents, ncol = 2)





