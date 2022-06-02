# ##############################################################################
# This script aims at analyzing the hierarchy of cells in each clonotype
# looking at the mutation levels
# ##############################################################################

## @knitr hierarchy_analysis

# Get the list of clonotypes to analyze
Idents( sc10x.rna.seurat) = "large.clonotype.similarity.symbol"
clonotype_to_analyze = na.omit( unique( as.character( Idents( sc10x.rna.seurat))))

# Parse the list of clonotype to analyse to build a phylogenetic tree for each of them
for( current_clonotype in clonotype_to_analyze){
  
  # Get the sequence and cell infos specific od the clonotype
  clonotype_info_df = heavy_and_light_chains_df[ which( heavy_and_light_chains_df$full.clonotype.similarity.symbol == current_clonotype), ]
  
  # Add the cell symbol (human readble barcode)
  clonotype_info_df$cell.symbol = paste0( current_clonotype, "_CELL", seq( 1, nrow( clonotype_info_df),1))
  
  # Add the tissue of origin of each cell and the corresponding tissue color
  clonotype_info_df$cell.tissue = sapply( clonotype_info_df$cell.bc, function( cellbc){
    return( cell_tissue_set[ cellbc])
  })
  
  # Save as fatsa file the sequences of the cells of the clonotype
  clonotype_info_file = file.path( PATH_ANALYSIS_OUTPUT, paste0( current_clonotype, "_HeavyAndLightChains_sequence.fasta"))
  write.fasta( sequences = as.list( clonotype_info_df$full.trimmed.chain), names = clonotype_info_df$cell.symbol,
               file.out = clonotype_info_file, nbchar = max( unlist( sapply( clonotype_info_df$full.trimmed.chain, nchar)))+1)
  
  # Read the just saved fatsa file to get the sequences in the right format
  sequences = readAAStringSet( clonotype_info_file)
  
  # Build a multiple alignment from the sequences list using the MUSCLE method (better for DNA than ClustalW)
  alignment = msa( sequences, method = "Muscle")
  
  # Build the distance matrix between sequences
  alignment.seqinr <- msaConvert( alignment, type="seqinr::alignment")
  distmat = dist.alignment( alignment.seqinr, "similarity")
  
  # Build the phylogenetic tree of the cells
  phylotree <- nj( distmat)
  
  # Plot the phylogenetic tree with cells colored by tissue
  colors = unlist( sapply( phylotree$tip.label, function( tip){
    return( tissue_palette[ clonotype_info_df[ which( clonotype_info_df$cell.symbol == tip), "cell.tissue"]])
  }))
  plot( phylotree, cex=0.6, main= paste( "Phylogenetic Tree of", current_clonotype), tip.color = colors)
  
  phylotree_mrca = mrca( phylotree, full = FALSE)
  getMRCA( phylotree, tip = c( "CLONOTYPE31_CELL6", "CLONOTYPE31_CELL30"))
  phylo_igraph = as.igraph.phylo( phylotree)
  
  plot( phylo_igraph)
  plot( phylo_igraph, layout=layout_with_kk, vertex.color="green")
  write_graph( phylo_igraph, file = file.path( PATH_ANALYSIS_OUTPUT, paste0( current_clonotype, "_graph.tab")), format = "edgelist")
  
}

