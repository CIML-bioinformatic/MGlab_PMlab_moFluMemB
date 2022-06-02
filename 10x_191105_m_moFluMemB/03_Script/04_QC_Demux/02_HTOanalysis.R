#########################################################################################################
# This script aims to demultiplex the samples (9 samples) using the HTO counts
#########################################################################################################

## @knitr hto_analysis

### Get the HTO data and manage them with RNA data
# ....................................................................

# Setup Seurat object with HTO data
sc10x.hto.data <- Read10X( data.dir = PATH_HTO_COUNT_TABLE, gene.column=1)

# Look at the distribution of the HTO values accross samples
hto_distribution_df = data.frame( sample = rownames( sc10x.hto.data),
                                  min.hto = apply( sc10x.hto.data, 1, min),
                                  median.hto = apply( sc10x.hto.data, 1, median),
                                  max.hto = apply( sc10x.hto.data, 1, max))

datatable( hto_distribution_df, rownames = FALSE)

# Check for cells with its maximum HTO assigned to unmapped
unmapped_cell_barcodes = apply( sc10x.hto.data, 2, function( barcode.hto){
  if( max( barcode.hto) == barcode.hto[ "unmapped"]){
    return( TRUE)
  }
  return( FALSE)
})
cat("<BR>Number of cellular barcodes with their maximum HTO unmapped:", sum( unmapped_cell_barcodes))

# Exclude the unmapped counts
filtered.sc10x.hto.data= as.matrix( sc10x.hto.data[ - which( row.names( sc10x.hto.data) %in% c( "unmapped")),])

# Looking at the cell barcodes present in RNA counts and HTO counts
common_cell_barcodes <- intersect( colnames( sc10x.rna.seurat), colnames( filtered.sc10x.hto.data))
cat("<br>Number of barcode in HTO:", length( colnames( filtered.sc10x.hto.data)))
cat("<br>Number of barcode in RNA:", length( colnames( sc10x.rna.seurat)))
cat("<br>Number of barcode in common between RNA and HTO:", length( common_cell_barcodes))

# Subset HTO counts by joint cell barcodes with RNA counts
filtered.sc10x.hto.data <- filtered.sc10x.hto.data[, common_cell_barcodes]

# Subset RNA seurat object by joint cell barcodes with HTO counts
sc10x.rna.seurat = subset( sc10x.rna.seurat, cells = common_cell_barcodes)

# Normalize RNA data with log normalization
sc10x.rna.seurat = NormalizeData( object = sc10x.rna.seurat,
                                  normalization.method = DATA_NORM_METHOD,
                                  scale.factor = DATA_NORM_SCALEFACTOR,
                                  verbose = .VERBOSE);

sc10x.rna.seurat = ScaleData( object    = sc10x.rna.seurat,
                              do.center = DATA_CENTER,
                              do.scale  = DATA_SCALE,
                              vars.to.regress = DATA_VARS_REGRESS,
                              verbose = .VERBOSE);

# Add HTO data as a new assay independent from RNA
sc10x.rna.seurat[[ "HTO"]] <- CreateAssayObject( counts = filtered.sc10x.hto.data)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
sc10x.rna.seurat <- NormalizeData( sc10x.rna.seurat, assay = "HTO", normalization.method = "CLR")

### Demultiplex the cells with MULTIseqDemux
# ....................................................................

cat("<H5>Demultiplex cells with MULTIseqdemux</H5>")

# Compute the demultiplexing
sc10x.rna.seurat <- MULTIseqDemux( sc10x.rna.seurat, assay = "HTO", autoThresh = TRUE)

# Provide a short name to samples for easy display
# ....................................................................

sc10x.rna.seurat[[ "MULTI_Sample"]] = apply( sc10x.rna.seurat[[ "MULTI_ID"]], 1, function( sample_name){
  
  dash_index = regexpr('-', sample_name)[1]
  if( dash_index > 0){
    return( substring( sample_name, 1, (dash_index-1)))
  }else{
    return( sample_name)
  }
})


Idents( sc10x.rna.seurat) <- "MULTI_Sample"

# Plot the Ridge plots per mouse
RidgePlot( sc10x.rna.seurat, features = rownames(sc10x.rna.seurat[["HTO"]])[1:3], assay = "HTO",group.by = "MULTI_Sample",ncol = 2)
RidgePlot( sc10x.rna.seurat, features = rownames(sc10x.rna.seurat[["HTO"]])[4:6], assay = "HTO",group.by = "MULTI_Sample",ncol = 2)
RidgePlot( sc10x.rna.seurat, features = rownames(sc10x.rna.seurat[["HTO"]])[7:9], assay = "HTO",group.by = "MULTI_Sample",ncol = 2)

# Plot the t-SNE embeddings with a distance matrix based on HTO counts
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = sc10x.rna.seurat, assay = "HTO"))))
sc10x.rna.seurat <- RunTSNE(sc10x.rna.seurat, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot( sc10x.rna.seurat, group.by = "MULTI_Sample", reduction="tsne")  + ggtitle("t-SNE of all cells based on MULTIseq demultiplexing")+ 
  theme(plot.title = element_text(lineheight=.8, face="bold"))

# Print the table of cells to sample association
datatable(data.frame((table(sc10x.rna.seurat$MULTI_Sample))), caption = "Sample association of barcodes")

# Plot the distribution of RNA counts along samples/doublet/Negative
VlnPlot( sc10x.rna.seurat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# Filter cells based on MULTIseqDemux results
# ...................................................................

cat("<H5>Filtering doublet and negative cells</H5>")

cat( "<br>Number of remaining cells before filtering:", length( Cells( sc10x.rna.seurat)))
sc10x.rna.seurat <- subset( sc10x.rna.seurat, idents = "Doublet", invert = TRUE)
cat( "<br>Number of remaining cells after filtering doublets :", length( Cells( sc10x.rna.seurat)))
sc10x.rna.seurat <- subset(x = sc10x.rna.seurat, idents = "Negative", invert = TRUE )
cat( "<br>Number of remaining cells after filtering Negative :", length( Cells( sc10x.rna.seurat)))

cat( "<br><br><b>Number of remaining cells in final result :", length( Cells( sc10x.rna.seurat)), "</b>")


# Save Seurat object to RDS file
# ...................................................................

saveRDS( sc10x.rna.seurat, file = file.path( PATH_ANALYSIS_OUTPUT, "sc10x.rna.seurat.RDS"))

