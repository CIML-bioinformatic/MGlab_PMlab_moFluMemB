# ##################################################################################
# This script aims to get the data from BCR analysis
# ##################################################################################


# =================================================================================================
# READ BCR DATA ON CELLS FROM BCR ANALYSIS RESULTS
# =================================================================================================

## @knitr prepare_bcr_data

cat("<H4>Read the data from BCR reconstruction and analysis</H4>")

# Load the data from the file produced by migmap at the end of the PMlab BCR analysis pipeline
bcr_migmap_heavy_df = read.table( PATH_BCR_MIGMAP_HEAVYCHAIN_ANALYSIS_FILE, sep=",", header = TRUE, stringsAsFactors = FALSE, quote='"')
bcr_migmap_light_df = read.table( PATH_BCR_MIGMAP_LIGHTCHAIN_ANALYSIS_FILE, sep=",", header = TRUE, stringsAsFactors = FALSE, quote='"')

# For light chain and heavy chain dataframe, get the barcode of the cells from the read name and put it as column in the corresponding dataframe
# The read name is composed of
# '>' + cellbarcode + "_TRINITY_" + some trinity infos 
# ........................................................................................
bcr_migmap_heavy_df$cell.bc = unlist( sapply( bcr_migmap_heavy_df$read.header, function( code){
  
  # get the first "_" index since the cell barcode is before
  index_underscore = regexpr( "_", code)
  if( index_underscore > 0){
    code = substr( code , start = 1, stop = index_underscore-1)
  }
  
  # get the first ">" index since the cell barcode is after
  index_greater = regexpr( ">", code)
  if( index_underscore >= 0){
    code = substr( code , start = index_greater+1, stop = nchar( code))
  }
  
  return( code)
}), use.names = FALSE)


bcr_migmap_light_df$cell.bc = unlist( sapply( bcr_migmap_light_df$read.header, function( code){
  
  # get the first "_" index since the cell barcode is before
  index_underscore = regexpr( "_", code)
  if( index_underscore > 0){
    code = substr( code , start = 1, stop = index_underscore-1)
  }
  
  # get the first ">" index since the cell barcode is after
  index_greater = regexpr( ">", code)
  if( index_underscore >= 0){
    code = substr( code , start = index_greater+1, stop = nchar( code))
  }
  
  return( code)
}), use.names = FALSE)


# For light chain and heavy chain dataframe, get the ID of the plate of each cell from the UniqueCellID and put it as column in the corresponding dataframe
# The read name is composed of
# '>' + cellbarcode + "_TRINITY_" + some trinity infos 
# plateID + "_" + cell_barcodeID + "_>" + cell_barcode + "_" + some Trinity info
# exemple : plate1_BC10_>TCTCACAC_TRINITY_DN30_c0_g1_i1
# ........................................................................................
bcr_migmap_heavy_df$cell.plate = unlist( sapply( bcr_migmap_heavy_df$UniqueCellID_UniqueIsoform, function( code){
  
  # get the first "_" index since the plate ID is before
  index_underscore = regexpr( "_", code)
  if( index_underscore > 0){
    code = substr( code , start = 1, stop = index_underscore-1)
  }
  
  return( code)
}), use.names = FALSE)

bcr_migmap_light_df$cell.plate = unlist( sapply( bcr_migmap_light_df$UniqueCellID_UniqueIsoform, function( code){
  
  # get the first "_" index since the plate ID is before
  index_underscore = regexpr( "_", code)
  if( index_underscore > 0){
    code = substr( code , start = 1, stop = index_underscore-1)
  }
  
  return( code)
}), use.names = FALSE)


# For light chain and heavy chain dataframe, get the ID of the barcode of each cell from the UniqueCellID and put it as column in the corresponding dataframe
# The read name is composed of
# '>' + cellbarcode + "_TRINITY_" + some trinity infos 
# plateID + "_" + cell_barcodeID + "_>" + cell_barcode + "_" + some Trinity info
# exemple : plate1_BC10_>TCTCACAC_TRINITY_DN30_c0_g1_i1
# ........................................................................................
bcr_migmap_heavy_df$cell.bcid = unlist( sapply( bcr_migmap_heavy_df$UniqueCellID_UniqueIsoform, function( code){
  
  # get the first "_" index since the barcodeID is after
  index_underscore = regexpr( "_", code)
  if( index_underscore > 0){
    code = substr( code , start = index_underscore+1, stop = nchar( code))
  }
  
  # get the first "_>" index since the barcodeID is before
  index_greater = regexpr( "_>", code)
  if( index_greater > 0){
    code = substr( code , start = 1, stop = index_greater-1)
  }
  
  return( code)
}), use.names = FALSE)

bcr_migmap_light_df$cell.bcid = unlist( sapply( bcr_migmap_light_df$UniqueCellID_UniqueIsoform, function( code){
  
  # get the first "_" index since the barcodeID is after
  index_underscore = regexpr( "_", code)
  if( index_underscore > 0){
    code = substr( code , start = index_underscore+1, stop = nchar( code))
  }
  
  # get the first "_>" index since the barcodeID is after
  index_greater = regexpr( "_>", code)
  if( index_greater > 0){
    code = substr( code , start = 1, stop = index_greater-1)
  }
  
  return( code)
}), use.names = FALSE)

# For light chain and heavy chain dataframe, merge the plateID and the cell barcode to have a unique identifier of the cell
# For light chain and heavy chain dataframe, merge the plateID and the cell barcodeID to have an other unique identifier of the cell
# ........................................................................................
bcr_migmap_heavy_df$cell.plate.bc = paste( bcr_migmap_heavy_df$cell.plate , bcr_migmap_heavy_df$cell.bc, sep="_")
bcr_migmap_light_df$cell.plate.bc = paste( bcr_migmap_light_df$cell.plate , bcr_migmap_light_df$cell.bc, sep="_")
bcr_migmap_heavy_df$cell.plate.bcid = paste( bcr_migmap_heavy_df$cell.plate , bcr_migmap_heavy_df$cell.bcid, sep="_")
bcr_migmap_light_df$cell.plate.bcid = paste( bcr_migmap_light_df$cell.plate , bcr_migmap_light_df$cell.bcid, sep="_")

# Count at the cells with heavy / light chain information
# ........................................................................................
cat("<br>Number of cells with heavy chain information in Migmap results:", nrow( bcr_migmap_heavy_df))
cat("<br>Number of cells with light chain information in Migmap results:", nrow( bcr_migmap_light_df))

# Get the complete light chain and heavy chain information for each cells in separated dataframes (one for light chain, one for heavy chain)
# clonotype = composition of detected V+D+J segments
# sequence = composition of FR1 + CDR1 + FR2+ CDR2 + FR3 + CDR3 + FR4 (only the FR4_AA_LENGTH_CUT first AA)
# ........................................................................................

# Compose the various strings
light_chain_df = data.frame( stringsAsFactors = FALSE)
heavy_chain_df = data.frame( stringsAsFactors = FALSE)

# Parse the line to analyse each cells
for( line_index in 1:nrow( bcr_migmap_light_df)){
  
  # Get the nucleotidic sequence : FR1 + CDR1 + FR2+ CDR2 + FR3 + CDR3 + FR4 (only the FR4_AA_LENGTH_CUT first AA)
  sequence = paste0( bcr_migmap_light_df$fr1nt[ line_index], bcr_migmap_light_df$cdr1nt[ line_index],
                    bcr_migmap_light_df$fr2nt[ line_index], bcr_migmap_light_df$cdr2nt[ line_index],
                    bcr_migmap_light_df$fr3nt[ line_index], bcr_migmap_light_df$cdr3nt[ line_index],
                    substr( bcr_migmap_light_df$fr4nt[ line_index], start =1 , stop = FR4_AA_LENGTH_CUT))
  trim.sequence = gsub("^N*", "", sequence)
  
  restricted.sequence = 
  
  # Get the clonotype : detected V+D+J segments
  clonotype = paste( bcr_migmap_light_df$cdr3nt[ line_index], bcr_migmap_light_df$v.segment[ line_index], bcr_migmap_light_df$d.segment[ line_index], bcr_migmap_light_df$j.segment[ line_index], sep="/")
  
  # Consider that a similar clonotype come from (same V+ same J + same length of CDR3)
  # This will be used for clonotype definition later
  clonotype.similarity = paste( bcr_migmap_light_df$v.segment[ line_index], bcr_migmap_light_df$j.segment[ line_index], nchar( bcr_migmap_light_df$cdr3nt[ line_index]), sep="/")
  
  # Get the isotype from the constant region
  isotype = bcr_migmap_light_df$CstRegion[ line_index]
  
  # Get the V-segment information
  v.segment = bcr_migmap_light_df$v.segment[ line_index]
  tokens = strsplit( v.segment, split='-', fixed=TRUE)
  if( length( tokens) == 1 && length( tokens[[1]]) == 2){
    v.segment.subgroup = tokens[[1]][1]
  }else{
    v.segment.subgroup = NA
  }
  
  # Get the V-segment exact information
  v.segment = bcr_migmap_light_df$v.segment[ line_index]
  tokens = strsplit( v.segment, split='*', fixed=TRUE)
  if( length( tokens) == 1 && length( tokens[[1]]) == 2){
    v.segment.exact = tokens[[1]][1]
  }else{
    v.segment.exact = v.segment
  }
  
  # Get the CDR3 length
  cdr3.length = nchar( bcr_migmap_light_df$cdr3nt[ line_index])
  
  # mutation number in sequence
  mutation.qual = nchar( bcr_migmap_light_df$mutations.qual[ line_index])
  
  # Store the information to the correct dataframe according chain type
  light_chain_df = rbind( light_chain_df, data.frame( cell.bc = bcr_migmap_light_df$cell.bc[ line_index],
                                                      cell.plate.bc = bcr_migmap_light_df$cell.plate.bc[ line_index],
                                                      cell.plate.bcid = bcr_migmap_light_df$cell.plate.bcid[ line_index],
                                                      clonotype = clonotype,
                                                      clonotype.similarity = clonotype.similarity,
                                                      sequence = sequence,
                                                      trim.sequence = trim.sequence,
                                                      isotype = isotype,
                                                      cdr3.length = cdr3.length, 
                                                      v.segment.exact = v.segment.exact,
                                                      v.segment = v.segment.subgroup,
                                                      mutation.qual = mutation.qual,
                                                      stringsAsFactors = FALSE))
}

# Parse the line to analyse each cells
for( line_index in 1:nrow( bcr_migmap_heavy_df)){
  
  # Get the nucleotidic sequence : FR1 + CDR1 + FR2+ CDR2 + FR3 + CDR3 + FR4 (only the FR4_AA_LENGTH_CUT first AA)
  sequence = paste0( bcr_migmap_heavy_df$fr1nt[ line_index], bcr_migmap_heavy_df$cdr1nt[ line_index],
                     bcr_migmap_heavy_df$fr2nt[ line_index], bcr_migmap_heavy_df$cdr2nt[ line_index],
                     bcr_migmap_heavy_df$fr3nt[ line_index], bcr_migmap_heavy_df$cdr3nt[ line_index],
                     substr( bcr_migmap_heavy_df$fr4nt[ line_index], start =1 , stop = FR4_AA_LENGTH_CUT))
  trim.sequence = gsub("^N*", "", sequence)
  
  # Get the clonotype : detected V+D+J segments
  clonotype = paste( bcr_migmap_heavy_df$cdr3nt[ line_index], bcr_migmap_heavy_df$v.segment[ line_index], bcr_migmap_heavy_df$d.segment[ line_index], bcr_migmap_heavy_df$j.segment[ line_index], sep="/")
  
  # Consider that a similar clonotype come from (same V+ same J + same length of CDR3)
  # This will be used for clonotype definition later
  clonotype.similarity = paste( bcr_migmap_heavy_df$v.segment[ line_index], bcr_migmap_heavy_df$j.segment[ line_index], nchar( bcr_migmap_heavy_df$cdr3nt[ line_index]), sep="/")
  
  # Get the isotype from the constant region
  isotype = bcr_migmap_heavy_df$CstRegion[ line_index]
  
  # Get the V-segment information
  v.segment = bcr_migmap_heavy_df$v.segment[ line_index]
  tokens = strsplit( v.segment, split='-', fixed=TRUE)
  if( length( tokens) == 1 && length( tokens[[1]]) == 2){
    v.segment.subgroup = tokens[[1]][1]
  }else{
    v.segment.subgroup = NA
  }
  
  # Get the V-segment exact information
  v.segment = bcr_migmap_heavy_df$v.segment[ line_index]
  tokens = strsplit( v.segment, split='*', fixed=TRUE)
  if( length( tokens) == 1 && length( tokens[[1]]) == 2){
    v.segment.exact = tokens[[1]][1]
  }else{
    v.segment.exact = v.segment
  }
  
  # Get the CDR3 length
  cdr3.length = nchar( bcr_migmap_heavy_df$cdr3nt[ line_index])
  
  # mutation number in sequence
  mutation.qual = nchar( bcr_migmap_heavy_df$mutations.qual[ line_index])
  
  # Store the information to the correct dataframe according chain type
  heavy_chain_df = rbind( heavy_chain_df, data.frame( cell.bc = bcr_migmap_heavy_df$cell.bc[ line_index],
                                                      cell.plate.bc = bcr_migmap_heavy_df$cell.plate.bc[ line_index],
                                                      cell.plate.bcid = bcr_migmap_heavy_df$cell.plate.bcid[ line_index],
                                                      clonotype = clonotype,
                                                      clonotype.similarity = clonotype.similarity,
                                                      sequence = sequence,
                                                      trim.sequence = trim.sequence,
                                                      isotype = isotype,
                                                      v.segment.exact = v.segment.exact,
                                                      v.segment = v.segment.subgroup,
                                                      cdr3.length = cdr3.length,
                                                      mutation.qual = mutation.qual,
                                                      stringsAsFactors = FALSE))
}

cat("<br>Number of barcodes with heavy chain information after sequence merging:", nrow( heavy_chain_df))
cat("<br>Number of barcodes with light chain information after sequence merging:", nrow( light_chain_df))

# Remove duplicated information (same barcode + sequence) in each dataframe
duplicated_barcode_heavy_chain = duplicated( paste( heavy_chain_df$cell.bc, heavy_chain_df$sequence))
duplicated_barcode_light_chain = duplicated( paste( light_chain_df$cell.bc, light_chain_df$sequence))

cat("<br>Number of barcodes with duplicated heavy chain information (same barcode+sequence):", sum( duplicated_barcode_heavy_chain))
cat("<br>Number of barcodes with duplicated light chain information (same barcode+sequence):", sum( duplicated_barcode_light_chain))

heavy_chain_df = heavy_chain_df[ !duplicated_barcode_heavy_chain, ]
light_chain_df = light_chain_df[ !duplicated_barcode_light_chain, ]

cat("<br>Number of barcodes with heavy chain information after duplication removal:",nrow( heavy_chain_df))
cat("<br>Number of barcodes with light chain information after duplication removal:", nrow( light_chain_df))

# Merge the heavy and light chains information together, keeping only the cells with both heavy and light chains information
heavy_and_light_chains_df = merge( heavy_chain_df, light_chain_df, by="cell.plate.bcid", suffixes=c( ".heavy", ".light"))
row.names( heavy_and_light_chains_df) = heavy_and_light_chains_df$cell.plate.bcid

# .............................................................
# Analyse the complete clonotype groups i.e. the group of
# cells with the same light chain and heavy chain clonotypes
# .............................................................

cat("<H4>Compose complete clonotype groups</H4>")

cat("<br>Number of barcodes with light and heavy chain information after data merging:", nrow( heavy_and_light_chains_df))

# Compose the heavy and light chain sequence in a single sequence
heavy_and_light_chains_df$full.clonotype = paste0( heavy_and_light_chains_df$clonotype.heavy, "<<>>", heavy_and_light_chains_df$clonotype.light)
heavy_and_light_chains_df$full.clonotype.similarity.mixed = paste0( heavy_and_light_chains_df$clonotype.similarity.heavy, "<<>>", heavy_and_light_chains_df$clonotype.similarity.light)
heavy_and_light_chains_df$full.chain = paste0( heavy_and_light_chains_df$sequence.heavy, heavy_and_light_chains_df$sequence.light)
heavy_and_light_chains_df$full.trimmed.chain = paste0( heavy_and_light_chains_df$trim.sequence.heavy, heavy_and_light_chains_df$trim.sequence.light)

# Save the full sequences to file (heavy chain + light chain)
heavy_and_light_chain_fasta_file = file.path( PATH_ANALYSIS_OUTPUT, "HeavyAndLightChains_sequence.fasta")
write.fasta( sequences = as.list( heavy_and_light_chains_df$full.trimmed.chain), names = heavy_and_light_chains_df$cell.plate.bcid,
             file.out = heavy_and_light_chain_fasta_file, nbchar = max( unlist( sapply( heavy_and_light_chains_df$full.trimmed.chain, nchar)))+1)

# Determine size of groups of same clonotype similarity
clonotype_similarity_set = unique( heavy_and_light_chains_df$full.clonotype.similarity.mixed)
CLONOTYPE_GROUP_LIST = list()
CLONOTYPE_GROUP_SIZE_SET = vector()
for( clonotype in clonotype_similarity_set){
  CLONOTYPE_GROUP_LIST[[ clonotype]] = heavy_and_light_chains_df[ which( heavy_and_light_chains_df$full.clonotype.similarity.mixed == clonotype), "cell.plate.bc"]
  CLONOTYPE_GROUP_SIZE_SET = append( CLONOTYPE_GROUP_SIZE_SET, length( CLONOTYPE_GROUP_LIST[[ clonotype]]))
}
names( CLONOTYPE_GROUP_SIZE_SET) = clonotype_similarity_set

# Look at the distribution of clone group sizes
cat("<br>Number of complete clonotype groups found (same light chain and heavy chain clonotype):", length( CLONOTYPE_GROUP_LIST))
datatable( as.data.frame( table( CLONOTYPE_GROUP_SIZE_SET)), colnames = c( "Number of cell in clonotype", "Nb of clonotypes"), rownames = FALSE,
           caption = "Distribution of clonotype group sizes")

# ...................................................
# Execute clusteOmega to analyse the distance between
# sequences of light+heavy chains
# ...................................................

cat("<H4>Analyse complete (light+heavy) chain sequence distance between cells</H4>")

# Excute ClustOmega on the sequences of Heavy and Light chains to compute the distance matrix and the alignments
# clustalo --full --in= PATH_EXPERIMENT_OUTPUT/08_BCRrepertoire/HeavyAndLightChains_sequence.fasta --distmat-out=PATH_EXPERIMENT_OUTPUT/08_BCRrepertoire/HeavyAndLightChains_distMat.csv --out=PATH_EXPERIMENT_OUTPUT/08_BCRrepertoire/HeavyAndLightChains_alignment.csv
# Do not execute ClustalO if result data alreadu exists (apart if the force mode is on).
# If the execution must be done, the possibly existing previous result files must be removed or ClustalO refuse to run

# -- Define the ClustalO result files
heavy_and_light_chain_distmat_file = file.path( PATH_ANALYSIS_OUTPUT, "HeavyAndLightChains_distMat.csv")
heavy_and_light_chain_alignment_file = file.path( PATH_ANALYSIS_OUTPUT, "HeavyAndLightChains_alignment.csv")

# -- Try to execute ClustalO
if( FORCE_CLUSTALO || !file.exists( heavy_and_light_chain_distmat_file) || !file.exists( heavy_and_light_chain_alignment_file)){
  if( file.exists( heavy_and_light_chain_distmat_file)){
    file_delete( heavy_and_light_chain_distmat_file)
  }
  if( file.exists( heavy_and_light_chain_alignment_file)){
    file_delete( heavy_and_light_chain_alignment_file)
  }
  cat("<br>Executing ClustalOmega to get distance between sequences")
  system2( command = "clustalo",
           args = c( " --full",
                   paste0( " --in=", heavy_and_light_chain_fasta_file),
                   paste0( " --distmat-out=", heavy_and_light_chain_distmat_file),
                   paste0( " --out=", heavy_and_light_chain_alignment_file)),
            wait = TRUE)
}else{
  cat("<br>Do not execute ClustalOmega and get already existing ClustalOmega results")
}

# Read the ClustalOmega distance matrix
heavy_and_light_chain_distmat_df = read.table( heavy_and_light_chain_distmat_file, header = FALSE, row.names = 1, skip = 1)
names( heavy_and_light_chain_distmat_df) = row.names( heavy_and_light_chain_distmat_df)

# Plot the distribution of distances found in the distance matrix
ggplot( data.frame( distance = as.numeric( unlist( heavy_and_light_chain_distmat_df)))) + 
        geom_histogram( aes( x=distance), fill="lightgreen") +
        ggtitle( "Distribution of distance between full BCR sequences (Heavy+Light)") +
        theme_classic()

# Look at the distance scores between cells in the same clonotype group
distances_list = list()
all_distances = vector()
for( clonotype in names( CLONOTYPE_GROUP_LIST)){
  cell_bc_set = CLONOTYPE_GROUP_LIST[[ clonotype]]
  if( length( cell_bc_set) > 1){
    pairs = combn( cell_bc_set, 2)
    distances = apply( pairs, 2, function( col){
      return( heavy_and_light_chain_distmat_df[ col[1], col[2]])
    })
    distances_list[[ clonotype]] = distances
    all_distances = append( all_distances, distances)
  }
}

ggplot( data.frame( distance = all_distances)) + 
  geom_histogram( aes( x=distance), fill="lightblue") +
  ggtitle( "Distribution of distance between full BCR sequences (Heavy+Light)\n between same clonotype cells") +
  theme_classic()

# Build a mapping from clonotype name to clonotype symbol (human readable)
cat("<H4>Map clonotype name to human readable symbol</H4>")
CLONOTYPE_NAME_TO_SYMBOL_MAPPING = data.frame( name = names( CLONOTYPE_GROUP_SIZE_SET), symbol = paste0( "CLONOTYPE", seq(1, length( CLONOTYPE_GROUP_SIZE_SET), 1)))
row.names( CLONOTYPE_NAME_TO_SYMBOL_MAPPING) = CLONOTYPE_NAME_TO_SYMBOL_MAPPING$name
datatable( CLONOTYPE_NAME_TO_SYMBOL_MAPPING, rownames = FALSE, caption = "Mapping of clonotype name to human readable symbol")

# Add the symbol of clonotype to the global dataframe
heavy_and_light_chains_df$full.clonotype.similarity.symbol.mixed = unlist( sapply( heavy_and_light_chains_df$full.clonotype.similarity.mixed, function( name){
  return( CLONOTYPE_NAME_TO_SYMBOL_MAPPING[ name, "symbol"])
}), use.names = FALSE)

# Save the global dataframe to file
write.table( heavy_and_light_chains_df, file = file.path( PATH_ANALYSIS_OUTPUT, "heavy_and_light_chains_allinfos.tsv"),
             sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

