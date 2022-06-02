# ##############################################################################
# This script aims to merge data from the seurat object built in previous
# analysis witht the data from clonotype analysis
# ##############################################################################

## @knitr merge_data

cat("<BR><H5>INFORMATION ON LOADED DATA FROM SEURAT ANALYSIS</H5>")

cat("<br><b>Number of cells in loaded Seurat data:", length( Cells( sc10x.rna.seurat)),"</b>")
cat("<br><b>Number of cells in Seurat data with heavy and light chain information:", length( intersect( heavy_and_light_chains_df$cell.bc, Cells( sc10x.rna.seurat))),"</b>")
cat("<HR>")

# ................................................................................................
# Compute several meta-data in the Seurat object
# ................................................................................................

# Add the metadata on the tissue
cell_tissue = unlist( sapply( sc10x.rna.seurat@meta.data$MULTI_Sample, function( sample){
  indexes = which( MOUSE_TO_SAMPLE_DF == sample, arr.ind = TRUE)
  if( nrow( indexes) == 1){
    return( row.names( MOUSE_TO_SAMPLE_DF)[ indexes[ 1, "row"]])
  }else{
    return( NA)
  }
}), use.names = FALSE)
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_tissue, 
                                col.name = "tissue")

# Add the metadata on the mouse
cell_mouse = unlist( sapply( sc10x.rna.seurat@meta.data$MULTI_Sample, function( sample){
  indexes = which( MOUSE_TO_SAMPLE_DF == sample, arr.ind = TRUE)
  if( nrow( indexes) == 1){
    return( names( MOUSE_TO_SAMPLE_DF)[ indexes[ 1, "col"]])
  }else{
    return( NA)
  }
}), use.names = FALSE)
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_mouse, 
                                col.name = "mouse")

# Add the metadata on isotypes (heavy and light chains)
cell_isotype_heavy_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "isotype.heavy"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_isotype_heavy_names, 
                                col.name = "isotype.heavy")

cell_isotype_light_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "isotype.light"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_isotype_light_names, 
                                col.name = "isotype.light")

# Add metadata on V-segment (heavy and light chains)
cell_vsegment_heavy_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "v.segment.heavy"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_vsegment_heavy_names, 
                                col.name = "v.segment.heavy")
cell_vsegment_light_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "v.segment.light"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_vsegment_light_names, 
                                col.name = "v.segment.light")

# Add metadata on CDR3 length (heavy and light chains)
cell_cdr3length_heavy_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "cdr3.length.heavy"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_cdr3length_heavy_names, 
                                col.name = "cdr3.length.heavy")
cell_cdr3length_light_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "cdr3.length.light"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_cdr3length_light_names, 
                                col.name = "cdr3.length.light")

# Add metadata on mutation qual (heavy and light chains)
cell_mutation_qual_heavy_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "mutation.qual.heavy"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_mutation_qual_heavy_names, 
                                col.name = "mutation.qual.heavy")
cell_mutation_qual_light_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "mutation.qual.light"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_mutation_qual_light_names, 
                                col.name = "mutation.qual.light")


# Add the metadata on clones
cell_clonotype_names = heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "full.clonotype.similarity"]
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_clonotype_names, 
                                col.name = "full.clonotype.similarity")

cell_clonotype_symbols = unname( sapply( cell_clonotype_names, function( clonotype){
  return( CLONOTYPE_NAME_TO_SYMBOL_MAPPING[ clonotype, "symbol"])
}))

sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = cell_clonotype_symbols, 
                                col.name = "full.clonotype.similarity.symbol")

# Select the largest clone groups and add it as metadata in seurat object
large_clonotype_set = names( sort( CLONOTYPE_GROUP_SIZE_SET, decreasing = TRUE))[ 1:11]
large_clonotypes_names = sapply( heavy_and_light_chains_df[ Cells( sc10x.rna.seurat), "full.clonotype.similarity"], function( clono){
  if( !is.na( clono) && clono %in% large_clonotype_set){
    return( clono)
  }else{
    return( NA)
  }
})
sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = large_clonotypes_names, 
                                col.name = "large.clonotype.similarity")

large_clonotypes_symbols = unname( sapply( large_clonotypes_names, function( clonotype){
  return( CLONOTYPE_NAME_TO_SYMBOL_MAPPING[ clonotype, "symbol"])
}))

sc10x.rna.seurat = AddMetaData( object = sc10x.rna.seurat, 
                                metadata = large_clonotypes_symbols, 
                                col.name = "large.clonotype.similarity.symbol")


# ................................................................................................
# Prepare data for analysis
# ................................................................................................

# Build a dataframe to merge umap embedding with various information (clonotype, isotype, sample, v.sgement...)
clonotype_seurat_df = data.frame( sc10x.rna.seurat@reductions$umap@cell.embeddings)
Idents( sc10x.rna.seurat) = "tissue"
clonotype_seurat_df$tissue = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "mouse"
clonotype_seurat_df$mouse = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$mouse = factor( clonotype_seurat_df$mouse, levels = sort( levels( clonotype_seurat_df$mouse)))
Idents( sc10x.rna.seurat) = "full.clonotype.similarity.symbol"
clonotype_seurat_df$full.clonotype.similarity.symbol = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "large.clonotype.similarity.symbol"
clonotype_seurat_df$large.clonotype.similarity.symbol = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "MULTI_Sample"
clonotype_seurat_df$sample = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "seurat_clusters"
clonotype_seurat_df$seurat_clusters = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "isotype.heavy"
clonotype_seurat_df$isotype.heavy = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "isotype.light"
clonotype_seurat_df$isotype.light = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "v.segment.heavy"
clonotype_seurat_df$v.segment.heavy = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "v.segment.light"
clonotype_seurat_df$v.segment.light = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( sc10x.rna.seurat) = "cdr3.length.heavy"
clonotype_seurat_df$cdr3.length.heavy = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$cdr3.length.heavy = factor( clonotype_seurat_df$cdr3.length.heavy, levels = sort( as.numeric( levels( clonotype_seurat_df$cdr3.length.heavy))))
Idents( sc10x.rna.seurat) = "cdr3.length.light"
clonotype_seurat_df$cdr3.length.light = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$cdr3.length.light = factor( clonotype_seurat_df$cdr3.length.light, levels = sort( as.numeric( levels( clonotype_seurat_df$cdr3.length.light))))
Idents( sc10x.rna.seurat) = "mutation.qual.heavy"
clonotype_seurat_df$mutation.qual.heavy = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$mutation.qual.heavy = factor( clonotype_seurat_df$mutation.qual.heavy, levels = sort( as.numeric( levels( clonotype_seurat_df$mutation.qual.heavy))))
Idents( sc10x.rna.seurat) = "mutation.qual.light"
clonotype_seurat_df$mutation.qual.light = Idents( sc10x.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$mutation.qual.light = factor( clonotype_seurat_df$mutation.qual.light, levels = sort( as.numeric( levels( clonotype_seurat_df$mutation.qual.light))))

# ................................................................................................
# Fix colors for each feature for a better display
# ................................................................................................

# Fix a color palette for the seurat clusters
if( !exists( "tissue_palette")){
  tissue_palette = brewer.pal( length( levels( clonotype_seurat_df$tissue)), "Dark2")
  names( tissue_palette) = levels( clonotype_seurat_df$tissue)
}

# Fix a color palette for the seurat clusters
if( !exists( "seurat_cluster_palette")){
  seurat_cluster_palette = brewer.pal( length( levels( clonotype_seurat_df$seurat_clusters)), "Set1")
  names( seurat_cluster_palette) = levels( clonotype_seurat_df$seurat_clusters)
}

# Fix a color palette for the largest clonotypes
if( !exists( "large_clonotype_palette")){
  large_clonotype_palette = brewer.pal( length( levels( clonotype_seurat_df$large.clonotype.similarity.symbol)), "Paired")
  names( large_clonotype_palette) = levels( clonotype_seurat_df$large.clonotype.similarity.symbol)
}

# Fix a color palette for the isotypes of heavy chains
if( !exists( "isotype_heavy_palette")){
  isotype_heavy_palette = brewer.pal( length( levels( clonotype_seurat_df$isotype.heavy)), "Set2")
  names( isotype_heavy_palette) = levels( clonotype_seurat_df$isotype.heavy)
}

# Fix a color palette for the isotypes of light chains
if( !exists( "isotype_light_palette")){
  isotype_light_palette = brewer.pal( length( levels( clonotype_seurat_df$isotype.light)), "Set2")
  names( isotype_light_palette) = levels( clonotype_seurat_df$isotype.light)
}

# Fix a color palette for the V-segment of heavy chains
if( !exists( "vsegment_heavy_palette")){
  vsegment_heavy_palette = colorRampPalette( brewer.pal( 8, "Dark2"))( length( levels( clonotype_seurat_df$v.segment.heavy)))
  names( vsegment_heavy_palette) = levels( clonotype_seurat_df$v.segment.heavy)
}

# Fix a color palette for the V-segment of light chains
if( !exists( "vsegment_light_palette")){
  vsegment_light_palette = colorRampPalette( brewer.pal( 8, "Dark2"))( length( levels( clonotype_seurat_df$v.segment.light)))
  names( vsegment_light_palette) = levels( clonotype_seurat_df$v.segment.light)
}

# ................................................................................................
# Plot the clusters on a UMAP embedding (for reference)
# ................................................................................................
print( ggplot() + 
         geom_point( data = clonotype_seurat_df, aes( x=UMAP_1, y=UMAP_2, col=seurat_clusters)) +
         scale_color_manual( values = seurat_cluster_palette) +
         theme_minimal() +
         ggtitle( "Clusters on UMAP embedding")
)



