# #########################################
# This script reads the sc10x data from
# QC and filtering analysis
# #########################################


# READ DATA
# ---------

## @knitr loadData

### Load the Seurat object from the QC and demultiplexing analysis
# .................................................................

custom.rna.seurat = readRDS( file.path( PATH_SEURAT_FILTERED_OBJECT, "custom.rna.seurat.RDS"))

cat("<br>Number of cells in loaded data:", length( Cells( custom.rna.seurat)))

### Remove the BCR genes from the Seurat object
# ...................................................

# Identify the list of BCR genes

igh_genes = rownames( custom.rna.seurat)[ grepl( "^Igh", rownames( custom.rna.seurat))]
igk_genes = rownames( custom.rna.seurat)[ grepl( "^Igk", rownames( custom.rna.seurat))]
igl_genes = rownames( custom.rna.seurat)[ grepl( "^Igl", rownames( custom.rna.seurat))]
bcr_genes = c( igh_genes, igk_genes, igl_genes)

cat("<br>Removing BCR genes. Number of genes removed:", length( bcr_genes))

# Identify the list of all features excluded BCR genes and keep only them in the Seurat object
genes_to_keep = setdiff( rownames( custom.rna.seurat), bcr_genes)
custom.rna.seurat = subset( custom.rna.seurat, features = genes_to_keep)


## Managing the list of monitored genes
# ...................................................

# Remove eventual NULL (empty) list elements from list of genes to be monitored
monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];

# Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
#matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( custom.rna.seurat))));
matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( custom.rna.seurat))));
monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));

# Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
MONITORED_GENES = relist( rownames( GetAssayData( custom.rna.seurat))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)

# Finally remove names that did not match
MONITORED_GENES = lapply( MONITORED_GENES, na.omit);

