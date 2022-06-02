# #################################################################################################
# This script aims to identify markers genes between tissues
# #################################################################################################

## @knitr identify_tissue_marker_genes

### Load the Seurat object from the QC and demultiplexing analysis
# .................................................................

sc10x.rna.seurat = readRDS( file.path( PATH_SEURAT_FILTERED_OBJECT, "sc10x.rna.seurat.RDS"))

cat("<br>Number of cells in loaded data:", length( Cells( sc10x.rna.seurat)))

### Filter some cells
# ...................................................

# Keep only the cells that were not affected by the stress induced by the cell dissociation process
for( CURRENT_TISSUE in names( SUBSAMPLE_BY_TISSUE)){
  dissociation_stress_cluster_number = DISSOCIATION_AFFECTED_CELLS_CLUSTER[[ CURRENT_TISSUE]]
  if( !is.na( dissociation_stress_cluster_number)){
    dissociation_stress_cluster_df = read.table( PATH_CLUSTER_MAPPING_FILE[[ CURRENT_TISSUE]], sep = "\t", header = TRUE)
    cells_stress_affected = dissociation_stress_cluster_df[ which( dissociation_stress_cluster_df$cluster == dissociation_stress_cluster_number), "cell.id"]
    sc10x.rna.seurat = subset( sc10x.rna.seurat, cells = setdiff( Cells( sc10x.rna.seurat), cells_stress_affected))
    cat("<br>Number of cells remaining after filtering cells affected by dissociation stress in ", CURRENT_TISSUE, ":", length( Cells( sc10x.rna.seurat)), sep="")
  }
}

### Remove the BCR genes from the Seurat object
# ...................................................

# Identify the list of BCR genes
igh_genes = rownames( sc10x.rna.seurat)[ grepl( "^Igh", rownames( sc10x.rna.seurat))]
igk_genes = rownames( sc10x.rna.seurat)[ grepl( "^Igk", rownames( sc10x.rna.seurat))]
igl_genes = rownames( sc10x.rna.seurat)[ grepl( "^Igl", rownames( sc10x.rna.seurat))]
bcr_genes = c( igh_genes, igk_genes, igl_genes)

cat("<br>Removing BCR genes. Number of genes removed:", length( bcr_genes))

# Identify the list of all features excluded BCR genes and keep only them in the Seurat object
genes_to_keep = setdiff( rownames( sc10x.rna.seurat), bcr_genes)
sc10x.rna.seurat = subset( sc10x.rna.seurat, features = genes_to_keep)

### Search for the marker genes of tissues
# ...................................................

# Add the metadata on the tissue
cell_tissue = unlist( sapply( sc10x.rna.seurat@meta.data$MULTI_Sample, function( sample){
  indexes = which( MOUSE_TO_SAMPLE_DF == sample, arr.ind = TRUE)
  if( nrow( indexes) == 1){
    return( row.names( indexes)[1])
  }else{
    return( NA)
  }
}), use.names = FALSE)
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_tissue, 
                                col.name = "tissue")

# Identify marker genes between tissues
Idents( sc10x.rna.seurat) <- "tissue"
markers = FindAllMarkers( object          = sc10x.rna.seurat,
                          test.use        = FINDMARKERS_METHOD,
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          verbose         = .VERBOSE);

#Print the datatable with all markers
datatable( markers[ which( markers$p_val_adj < FINDMARKERS_PVAL_THR), ], caption = paste( "All markers by clusters with adj.pval <", FINDMARKERS_PVAL_THR))

