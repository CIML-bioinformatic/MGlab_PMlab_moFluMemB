# #########################################
# This script reads the sc10x data from
# QC and filtering analysis
# #########################################


# READ DATA
# ---------

## @knitr loadData

### Load the Seurat object from the QC and demultiplexing analysis
# .................................................................

sc10x.rna.seurat = readRDS( file.path( PATH_SEURAT_FILTERED_OBJECT, "sc10x.rna.seurat.RDS"))

cat("<br>Number of cells in loaded data:", length( Cells( sc10x.rna.seurat)))

### Filter some cells
# ...................................................

# Keep only the cells from the selected tissue
cells_to_keep = row.names( sc10x.rna.seurat@meta.data)[ sc10x.rna.seurat@meta.data[ , "MULTI_Sample"] %in% SUBSAMPLE_BY_TISSUE[[ CURRENT_TISSUE]] ]
sc10x.rna.seurat = subset( sc10x.rna.seurat, cells = cells_to_keep)
cat("<br>Number of cells after tissue (",  CURRENT_TISSUE, ") selection: ", length( Cells( sc10x.rna.seurat)), sep="")

# Keep only the cells that were not affected by the stress induced by the cell dissociation process
dissociation_stress_cluster_number = DISSOCIATION_AFFECTED_CELLS_CLUSTER[[ CURRENT_TISSUE]]
if( !is.na( dissociation_stress_cluster_number)){
  dissociation_stress_cluster_df = read.table( PATH_CLUSTER_MAPPING_FILE[[ CURRENT_TISSUE]], sep = "\t", header = TRUE)
  cells_not_stress_affected = dissociation_stress_cluster_df[ which( dissociation_stress_cluster_df$cluster != dissociation_stress_cluster_number), "cell.id"]
  sc10x.rna.seurat = subset( sc10x.rna.seurat, cells = cells_not_stress_affected)
  cat("<br>Number of cells remaining after filtering cells affected by dissociation stress : ", length( Cells( sc10x.rna.seurat)), sep="")
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


## Managing the list of monitored genes
# ...................................................

# Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x.rna.seurat))));
matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x.rna.seurat))));
monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));

# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MONITORED_GENES = relist( rownames( GetAssayData( sc10x.rna.seurat))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)

# Finally remove names that did not match
MONITORED_GENES = lapply( MONITORED_GENES, na.omit);

## Export normalized expression table
# .....................................
write.table( GetAssayData(object = sc10x.rna.seurat, assay = "RNA", slot = "data"),
             file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "normalizedExpressions_", CURRENT_TISSUE, ".tsv")),
             row.names = TRUE, col.names = NA, sep="\t", quote = FALSE)

