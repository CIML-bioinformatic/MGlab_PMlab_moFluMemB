# ##############################################################################
# This script aims to merge data from the seurat object built in previous
# analysis witht the data from clonotype analysis
# ##############################################################################

## @knitr merge_data

cat("<BR><H5>INFORMATION ON LOADED DATA FROM SEURAT ANALYSIS</H5>")

cat("<br><b>Number of cells in loaded Seurat data:", length( Cells( custom.rna.seurat)),"</b>")
cat("<br><b>Number of cells in Seurat data with heavy and light chain information:", length( intersect( heavy_and_light_chains_df$cell.plate.bcid, Cells( custom.rna.seurat))),"</b>")
cat("<HR>")


custom.rna.seurat = subset( custom.rna.seurat, cells = intersect( heavy_and_light_chains_df$cell.plate.bcid, Cells( custom.rna.seurat)))

# ...................................................
# Load information on cluster marker genes 
# as defined by previous analysis
# ...................................................
if( !exists( "current_tissue")){
  CLUSTER_MARKER_GENES_DF = data.frame()
  for( tissue in names( PATH_CLUSTER_MARKER_GENES_PER_TISSUE)){
    current_table = read.table( PATH_CLUSTER_MARKER_GENES_PER_TISSUE[[ tissue]], sep="\t", header=TRUE)
    current_table$tissue = rep( tissue, nrow( current_table))
    CLUSTER_MARKER_GENES_DF = rbind( CLUSTER_MARKER_GENES_DF, current_table)  
  }
}else{
  CLUSTER_MARKER_GENES_DF = read.table( PATH_CLUSTER_MARKER_GENES_PER_TISSUE[[ current_tissue]], sep="\t", header=TRUE)
  CLUSTER_MARKER_GENES_DF$tissue = rep( current_tissue, nrow( CLUSTER_MARKER_GENES_DF))
}



# ................................................................................................
# Compute several meta-data in the Seurat object
# ................................................................................................

# Add the metadata on isotypes (heavy and light chains)
cell_isotype_heavy_names = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "isotype.heavy"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_isotype_heavy_names, 
                                col.name = "isotype.heavy")

cell_isotype_light_names = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "isotype.light"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_isotype_light_names, 
                                col.name = "isotype.light")

# Add metadata on V-segment (heavy and light chains)
cell_vsegment_heavy_names = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "v.segment.heavy"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_vsegment_heavy_names, 
                                col.name = "v.segment.heavy")
cell_vsegment_light_names = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "v.segment.light"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_vsegment_light_names, 
                                col.name = "v.segment.light")

# Add metadata on CDR3 length (heavy and light chains)
cell_cdr3length_heavy_names = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "cdr3.length.heavy"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_cdr3length_heavy_names, 
                                col.name = "cdr3.length.heavy")
cell_cdr3length_light_names = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "cdr3.length.light"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_cdr3length_light_names, 
                                col.name = "cdr3.length.light")

# Add metadata on mutation qual (heavy and light chains)
cell_mutation_qual_heavy_names = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "mutation.qual.heavy"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_mutation_qual_heavy_names, 
                                col.name = "mutation.qual.heavy")
cell_mutation_qual_light_names = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "mutation.qual.light"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_mutation_qual_light_names, 
                                col.name = "mutation.qual.light")


# Add the metadata on clones

# -- Add the true clonotype name (composition of VDJ genes) NOT distinguishing the mouse of origin
cell_clonotype_names_mixed = heavy_and_light_chains_df[ Cells( custom.rna.seurat), "full.clonotype.similarity.mixed"]
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_clonotype_names_mixed, 
                                col.name = "full.clonotype.similarity.mixed")

# -- Convert the true clontype name to a clonotype symbol, NOT distinguishing the mouse of origin
cell_clonotype_symbols_mixed = unname( sapply( cell_clonotype_names_mixed, function( clonotype){
  return( CLONOTYPE_NAME_TO_SYMBOL_MAPPING[ clonotype, "symbol"])
}))
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = cell_clonotype_symbols_mixed, 
                                col.name = "full.clonotype.similarity.symbol.mixed")

# -- Convert the true clontype name to a clonotype symbol, distinguishing the mouse of origin
Idents( custom.rna.seurat) = "mouse"
cell_clonotype_names = paste( Idents( custom.rna.seurat), cell_clonotype_names_mixed, sep=".")
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                 metadata = cell_clonotype_names, 
                                 col.name = "full.clonotype.similarity")

cell_clonotype_symbols = paste( Idents( custom.rna.seurat), cell_clonotype_symbols_mixed, sep=".")
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                 metadata = cell_clonotype_symbols, 
                                 col.name = "full.clonotype.similarity.symbol")


# Select the largest clone groups and add it as metadata in Seurat object
# A large clonotype is clonotype that have, in the current tissue tissue (if so), at least 10 cells
# OR in all tissues, at least 20 cells

# -- retrieve the large clonotypes from all tissues
cat("<BR><b>NOTE : Defining LARGE CLONOTYPE as clonotypes with more than", GLOBAL_LARGE_CLUSTER_THRESHOLD, "cells in all tissues</b>")
if( !exists( "current_tissue")){
  global_large_clonotype_set = unique( names( which( table( cell_clonotype_names) >= GLOBAL_LARGE_CLUSTER_THRESHOLD)))
}

# --build the large clonotype factor for the Seurat object, based on distinguished mouse symbol

# -- Add the large clonotype name to Seurat object
Idents( custom.rna.seurat) = "full.clonotype.similarity"
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = sapply( Idents( custom.rna.seurat), function( clono){
                                  clono = as.character( clono)
                                  if( !is.na( clono) && clono %in% global_large_clonotype_set){
                                    return( clono)
                                  }else{
                                    return( NA)
                                  }
                                }), 
                                col.name = "large.clonotype.similarity")

# -- Add the large clonotype symbol to Seurat object
Idents( custom.rna.seurat) = "large.clonotype.similarity"
large_clonotype_indexes = which( !is.na( Idents( custom.rna.seurat)))
Idents( custom.rna.seurat) = "full.clonotype.similarity.symbol"
large.clonotype.similarity.symbol = rep( NA, length( Idents( custom.rna.seurat)))
large.clonotype.similarity.symbol[ large_clonotype_indexes] = as.character( Idents( custom.rna.seurat)[ large_clonotype_indexes])
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, 
                                metadata = large.clonotype.similarity.symbol, 
                                col.name = "large.clonotype.similarity.symbol")


# ................................................................................................
# Prepare data for analysis
# ................................................................................................

# Build a dataframe to merge umap embedding with various information (clonotype, isotype, sample, v.sgement...)
clonotype_seurat_df = data.frame( custom.rna.seurat@reductions$umap@cell.embeddings)
Idents( custom.rna.seurat) = "tissue"
clonotype_seurat_df$tissue = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "plate"
clonotype_seurat_df$plate = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "mouse"
clonotype_seurat_df$sample = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "mouse"
clonotype_seurat_df$mouse = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$mouse = factor( clonotype_seurat_df$mouse, levels = sort( levels( clonotype_seurat_df$mouse)))
Idents( custom.rna.seurat) = "full.clonotype.similarity.symbol"
clonotype_seurat_df$full.clonotype.similarity.symbol = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "large.clonotype.similarity.symbol"
clonotype_seurat_df$large.clonotype.similarity.symbol = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "seurat_clusters"
clonotype_seurat_df$seurat_clusters = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "isotype.heavy"
clonotype_seurat_df$isotype.heavy = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "isotype.light"
clonotype_seurat_df$isotype.light = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "v.segment.heavy"
clonotype_seurat_df$v.segment.heavy = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "v.segment.light"
clonotype_seurat_df$v.segment.light = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
Idents( custom.rna.seurat) = "cdr3.length.heavy"
clonotype_seurat_df$cdr3.length.heavy = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$cdr3.length.heavy = factor( clonotype_seurat_df$cdr3.length.heavy, levels = sort( as.numeric( levels( clonotype_seurat_df$cdr3.length.heavy))))
Idents( custom.rna.seurat) = "cdr3.length.light"
clonotype_seurat_df$cdr3.length.light = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$cdr3.length.light = factor( clonotype_seurat_df$cdr3.length.light, levels = sort( as.numeric( levels( clonotype_seurat_df$cdr3.length.light))))
Idents( custom.rna.seurat) = "mutation.qual.heavy"
clonotype_seurat_df$mutation.qual.heavy = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$mutation.qual.heavy = factor( clonotype_seurat_df$mutation.qual.heavy, levels = sort( as.numeric( levels( clonotype_seurat_df$mutation.qual.heavy))))
Idents( custom.rna.seurat) = "mutation.qual.light"
clonotype_seurat_df$mutation.qual.light = Idents( custom.rna.seurat)[ row.names( clonotype_seurat_df)]
clonotype_seurat_df$mutation.qual.light = factor( clonotype_seurat_df$mutation.qual.light, levels = sort( as.numeric( levels( clonotype_seurat_df$mutation.qual.light))))

# Write the merged information to file
if( !exists( "current_tissue")){
  write.table( clonotype_seurat_df, file = file.path( PATH_ANALYSIS_OUTPUT, "clonotype_and_seurat_merged_data.tsv"),
             sep="\t", col.names = NA, quote = FALSE)
}else{
  write.table( clonotype_seurat_df, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( "clonotype_and_seurat_merged_data_", current_tissue, ".tsv")),
               sep="\t", col.names = NA, quote = FALSE)
  
}


# ................................................................................................
# Fix colors for each feature for a better display
# ................................................................................................

# Define large color palette of colors
if( !exists( "color_vector")){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  color_vector = unlist( mapply( brewer.pal, qual_col_pals$maxcolors, rownames( qual_col_pals)))
}

# Fix a color palette for the seurat clusters
if( !exists( "tissue_palette")){
  tissue_palette = "lightblue"
  names( tissue_palette) = levels( clonotype_seurat_df$tissue)
}

# Fix a color palette for the seurat clusters
if( !exists( "seurat_cluster_palette")){
  seurat_cluster_palette = brewer.pal( length( levels( clonotype_seurat_df$seurat_clusters)), "Set1")
  names( seurat_cluster_palette) = levels( clonotype_seurat_df$seurat_clusters)
}

# Fix a color palette for the largest clonotypes
if( !exists( "large_clonotype_palette")){
  large_clonotype_palette = color_vector [1: length( levels( clonotype_seurat_df$large.clonotype.similarity.symbol))]
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



