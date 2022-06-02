# ##############################################################################
# This script aims at analyzing the FACS phenotypes in relation with the
# BCR data
# ##############################################################################

## @knitr analyze_bcrrepertoire_to_phenotype

# Read the meta-data and phenotype data from single file

cat("<BR>Reading meta-data and phenotype data from", PATH_CELL_PREPROCESSING_METADATA)
metadata_and_phenotype_df = read.table( PATH_CELL_PREPROCESSING_METADATA, sep=",", quote = '"', header = TRUE)
metadata_and_phenotype_df$Plate = factor( metadata_and_phenotype_df$Plate, levels = paste0( "plate", seq( 1, length( unique( metadata_and_phenotype_df$Plate)))))

# Change the cell plate+barcodeid to match those of clonotype data
metadata_and_phenotype_df$plate.bcid = gsub( paste0( SAMPLE_NAME, "_"), "", row.names( metadata_and_phenotype_df))

# Change scale for asinh of HA phenotype
metadata_and_phenotype_df$asinh_50_CXCR3.BV421 = asinh( metadata_and_phenotype_df$Comp.460_50.PB....CXCR3.BV421 / 50)
metadata_and_phenotype_df$asinh_20_CXCR3.BV421 = asinh( metadata_and_phenotype_df$Comp.460_50.PB....CXCR3.BV421 / 20)

# Look at the common information between clonotype and phenotype
cat("<BR>Number of cells in the meta-data/phenotype data :", nrow( metadata_and_phenotype_df))
cat("<BR>Number of cells in the clonotype data :", nrow( clonotype_seurat_df))
cat("<BR>Number of cells in both data:", length( intersect( row.names( clonotype_seurat_df), metadata_and_phenotype_df$plate.bcid)))

# Limit the phenotype information to the cell with also clonotype information
metadata_and_phenotype_df = metadata_and_phenotype_df[ which( metadata_and_phenotype_df$plate.bcid %in% row.names( clonotype_seurat_df)), ]

# =========================================================================
# Looking at Inde sort data on all cells
# =========================================================================

cat("<H4>Relation betwen CXCR3, CCR6 and plates</H4>")

# Look at CXRCR3 x CCR6 and their crossrelation with plates
# ..........................................................

# -- Look at CXCR3
ggplot( metadata_and_phenotype_df) + geom_density( aes_string( x=PHENOTYPE_CXCR3_INDEXSORT_SOURCE)) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_CXCR3_THRESHOLD, col="Red") + theme_light()
ggplot( metadata_and_phenotype_df) + geom_density( aes_string( x=PHENOTYPE_CXCR3_INDEXSORT_SOURCE, col="Plate"))  + 
  geom_vline( xintercept=PHENOTYPE_ASINH_CXCR3_THRESHOLD, col="Red") +
  facet_wrap( .~Plate) + theme_light() + theme( legend.position = "None")

# -- Look at CCR6
ggplot( metadata_and_phenotype_df) + geom_density( aes_string( x=PHENOTYPE_CCR6_INDEXSORT_SOURCE))  + 
  geom_vline( xintercept=PHENOTYPE_ASINH_CCR6_THRESHOLD, col="Blue") + theme_light()
ggplot( metadata_and_phenotype_df) + geom_density( aes_string( x=PHENOTYPE_CCR6_INDEXSORT_SOURCE, col="Plate"))  + 
  geom_vline( xintercept=PHENOTYPE_ASINH_CCR6_THRESHOLD, col="Blue")  +
  facet_wrap( .~Plate) + theme_light() + theme( legend.position = "None")


cat("<H4>Relation betwen CXCR3 x CCR6 and plates</H4>")

# -- Look at CXCR3 x CCR6
ggplot( metadata_and_phenotype_df) + 
  geom_point( aes_string( x=PHENOTYPE_CXCR3_INDEXSORT_SOURCE, y = PHENOTYPE_CCR6_INDEXSORT_SOURCE, col="Plate")) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_CXCR3_THRESHOLD, col="Red") +
  geom_hline( yintercept=PHENOTYPE_ASINH_CCR6_THRESHOLD, col="Blue") + theme_light()
ggplot( metadata_and_phenotype_df) + 
  geom_point( aes_string( x=PHENOTYPE_CXCR3_INDEXSORT_SOURCE, y = PHENOTYPE_CCR6_INDEXSORT_SOURCE, col="Plate")) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_CXCR3_THRESHOLD, col="Red") +
  geom_hline( yintercept=PHENOTYPE_ASINH_CCR6_THRESHOLD, col="Blue") +
  facet_wrap( .~Plate) + theme_light() + theme( legend.position = "None")

# Look at HA x NP and their crossrelation with plates
# .....................................................

cat("<H4>Relation betwen HA, NP and plates</H4>")

# -- Look at HA
ggplot( metadata_and_phenotype_df) + geom_density( aes_string( x=PHENOTYPE_HA_INDEXSORT_SOURCE)) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_HA_THRESHOLD, col="Red") + theme_light()
ggplot( metadata_and_phenotype_df) + geom_density( aes_string( x=PHENOTYPE_HA_INDEXSORT_SOURCE, col="Plate")) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_HA_THRESHOLD, col="Red") +
  facet_wrap( .~Plate) + theme_light() + theme( legend.position = "None")

# -- Look at NP
ggplot( metadata_and_phenotype_df) + geom_density( aes_string( x=PHENOTYPE_NP_INDEXSORT_SOURCE)) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_NP_THRESHOLD, col="Blue") + theme_light()
ggplot( metadata_and_phenotype_df) + geom_density( aes_string( x=PHENOTYPE_NP_INDEXSORT_SOURCE, col="Plate")) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_NP_THRESHOLD, col="Blue") +
  facet_wrap( .~Plate) + theme_light() + theme( legend.position = "None")

cat("<H5>Look at relation betwen HA x NP and plates</H5>")

# -- Look at HP x NA
ggplot( metadata_and_phenotype_df) + 
  geom_point( aes_string( x=PHENOTYPE_HA_INDEXSORT_SOURCE, y = PHENOTYPE_NP_INDEXSORT_SOURCE, col="Plate")) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_HA_THRESHOLD, col="Red") +
  geom_hline( yintercept=PHENOTYPE_ASINH_NP_THRESHOLD, col="Blue") + theme_light()
ggplot( metadata_and_phenotype_df) + 
  geom_point( aes_string( x=PHENOTYPE_HA_INDEXSORT_SOURCE, y = PHENOTYPE_NP_INDEXSORT_SOURCE, col="Plate")) + 
  geom_vline( xintercept=PHENOTYPE_ASINH_HA_THRESHOLD, col="Red") +
  geom_hline( yintercept=PHENOTYPE_ASINH_NP_THRESHOLD, col="Blue") +
  facet_wrap( .~Plate) + theme_light() + theme( legend.position = "None")


# Place gates on  CXCR3xCCR6 and HAxNP data to generate phenotype names
# ......................................................................

cat("<H4>Apply gates in CXCR3xCCR6 and on HAxNP and look at their interactions</H4>")

# -- Give names to cell phenoytpes using CXCR3 data gates

metadata_and_phenotype_df$CXCR3.phenotype= unlist( apply( metadata_and_phenotype_df, 1, function( row){
  cxcr3.value = as.numeric( as.character( row[ PHENOTYPE_CXCR3_INDEXSORT_SOURCE]))
  
  if( cxcr3.value < PHENOTYPE_ASINH_CXCR3_THRESHOLD){
    return( PHENOTYPE_CXCR3_NAMES[ "Neg"])
  }else {
    return( PHENOTYPE_CXCR3_NAMES[ "Pos"])
  }
}), use.names = FALSE)

# -- Give names to cell phenoytpes using CCR6 data gates

metadata_and_phenotype_df$CCR6.phenotype= unlist( apply( metadata_and_phenotype_df, 1, function( row){
  ccr6.value = as.numeric( as.character( row[ PHENOTYPE_CCR6_INDEXSORT_SOURCE]))
  
  if( ccr6.value < PHENOTYPE_ASINH_CCR6_THRESHOLD){
    return( PHENOTYPE_CCR6_NAMES[ "Neg"])
  }else {
    return( PHENOTYPE_CCR6_NAMES[ "Pos"])
  }
}), use.names = FALSE)

# -- Give names to cell phenoytpes using CXCR3 x CCR6 data gates

metadata_and_phenotype_df$CXCR3xCCR6= unlist( apply( metadata_and_phenotype_df, 1, function( row){
  cxcr3.value = as.numeric( as.character( row[ PHENOTYPE_CXCR3_INDEXSORT_SOURCE]))
  ccr6.value = as.numeric( as.character( row[ PHENOTYPE_CCR6_INDEXSORT_SOURCE]))
  
  if( cxcr3.value < PHENOTYPE_ASINH_CXCR3_THRESHOLD && ccr6.value < PHENOTYPE_ASINH_CCR6_THRESHOLD){
    return( PHENOTYPE_CXCR3xCCR6_NAMES[ "NegNeg"])
  }else if( cxcr3.value >= PHENOTYPE_ASINH_CXCR3_THRESHOLD && ccr6.value < PHENOTYPE_ASINH_CCR6_THRESHOLD){
    return( PHENOTYPE_CXCR3xCCR6_NAMES[ "PosNeg"])
  }else if( cxcr3.value < PHENOTYPE_ASINH_CXCR3_THRESHOLD && ccr6.value >= PHENOTYPE_ASINH_CCR6_THRESHOLD){
    return( PHENOTYPE_CXCR3xCCR6_NAMES[ "NegPos"])
  }else{
  return( PHENOTYPE_CXCR3xCCR6_NAMES[ "PosPos"])
  }
}), use.names = FALSE)
metadata_and_phenotype_df$CXCR3xCCR6 = factor( metadata_and_phenotype_df$CXCR3xCCR6, levels= PHENOTYPE_CXCR3xCCR6_NAMES_ORDER)

# -- Give names to cell phenoytpes using HA data gates
metadata_and_phenotype_df$HA.phenotype= unlist( apply( metadata_and_phenotype_df, 1, function( row){
  ha.value = as.numeric( as.character( row[ PHENOTYPE_HA_INDEXSORT_SOURCE]))
  
  if( ha.value < PHENOTYPE_ASINH_HA_THRESHOLD){
    return( PHENOTYPE_HA_NAMES[ "Neg"])
  }else {
    return( PHENOTYPE_HA_NAMES[ "Pos"])
  }
}), use.names = FALSE)


# -- Give names to cell phenoytpes using NP data gates
metadata_and_phenotype_df$NP.phenotype= unlist( apply( metadata_and_phenotype_df, 1, function( row){
  np.value = as.numeric( as.character( row[ PHENOTYPE_NP_INDEXSORT_SOURCE]))
  
  if( np.value < PHENOTYPE_ASINH_NP_THRESHOLD){
    return( PHENOTYPE_NP_NAMES[ "Neg"])
  }else {
    return( PHENOTYPE_NP_NAMES[ "Pos"])
  }
}), use.names = FALSE)

# -- Give names to cell phenoytpes using HA x NP data gates

metadata_and_phenotype_df$HAxNP= unlist( apply( metadata_and_phenotype_df, 1, function( row){
  ha.value = as.numeric( as.character( row[ PHENOTYPE_HA_INDEXSORT_SOURCE]))
  np.value = as.numeric( as.character( row[ PHENOTYPE_NP_INDEXSORT_SOURCE]))
  
  if( ha.value < PHENOTYPE_ASINH_HA_THRESHOLD && np.value < PHENOTYPE_ASINH_NP_THRESHOLD){
    return( PHENOTYPE_HAxNP_NAMES[ "NegNeg"])
  }else if( ha.value >= PHENOTYPE_ASINH_HA_THRESHOLD && np.value < PHENOTYPE_ASINH_NP_THRESHOLD){
    return( PHENOTYPE_HAxNP_NAMES[ "PosNeg"])
  }else if( ha.value < PHENOTYPE_ASINH_HA_THRESHOLD && np.value >= PHENOTYPE_ASINH_NP_THRESHOLD){
    return( PHENOTYPE_HAxNP_NAMES[ "NegPos"])
  }else{
    return( PHENOTYPE_HAxNP_NAMES[ "PosPos"])
  }
}), use.names = FALSE)
metadata_and_phenotype_df$HAxNP = factor( metadata_and_phenotype_df$HAxNP, levels= PHENOTYPE_HAxNP_NAMES_ORDER)

# Look at the chi2 statistics between CXCR3xCCR6 phenotypes against plates
plate_CXCR3xCCR6_wide_df = as.data.frame.matrix( table( metadata_and_phenotype_df$CXCR3xCCR6, metadata_and_phenotype_df$Plate))
chisq_plate_CXCR3xCCR6 = chisq.test(  plate_CXCR3xCCR6_wide_df)
cat("<BR>")
datatable( t( plate_CXCR3xCCR6_wide_df), caption = "Number of cells of each CXCR3xCCR6 phenotype for each plate")
cat("<BR>")
pheatmap( plate_CXCR3xCCR6_wide_df, cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells of each\nHAxNP phenotype for each plate",
          display_numbers = TRUE, number_format = "%d")
pheatmap( chisq_plate_CXCR3xCCR6$residuals, cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nnumber of cells of each\nCXCR3xCCR6 phenotype for each plate")

# Look at the chi2 statistics between HAxNP phenotypes against plates
plate_HAxNP_wide_df = as.data.frame.matrix( table( metadata_and_phenotype_df$HAxNP, metadata_and_phenotype_df$Plate))
chisq_plate_HAxNP = chisq.test(  plate_HAxNP_wide_df)
cat("<BR>")
datatable( t( plate_HAxNP_wide_df), caption = "Number of cells of each HAxNP phenotype for each plate")
cat("<BR>")
pheatmap( plate_HAxNP_wide_df, cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells of each\nHAxNP phenotype for each plate",
          display_numbers = TRUE, number_format = "%d")
pheatmap( chisq_plate_HAxNP$residuals, cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nnumber of cells of each\nHAxNP phenotype for each plate")

# Look at the chi2 statistics between CXCR3xCCR6 phenotypes against HAxNP phenotypes
HAxNP_CXCR3xCCR6_wide_df = as.data.frame.matrix( table( metadata_and_phenotype_df$HAxNP, metadata_and_phenotype_df$CXCR3xCCR6))
chisq_HAxNP_CXCR3xCCR6 = chisq.test(  HAxNP_CXCR3xCCR6_wide_df)
cat("<BR>")
datatable( HAxNP_CXCR3xCCR6_wide_df, caption = "Number of cells of each CXCR3xCCR6 phenotype for each HAxNP phenotype")
cat("<BR>")
pheatmap( HAxNP_CXCR3xCCR6_wide_df, cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells of each\nCXCR3xCCR6 phenotype for each HAxNP phenotype",
          display_numbers = TRUE, number_format = "%d")
pheatmap( floor( chisq_HAxNP_CXCR3xCCR6$expected), cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Expected number of cells of chi2 test of each\nCXCR3xCCR6 phenotype for each HAxNP phenotype",
          display_numbers = TRUE, number_format = "%d")
pheatmap( chisq_HAxNP_CXCR3xCCR6$residuals, cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nnumber of cells of each\nCXCR3xCCR6 phenotype for each HAxNP phenotype",
          display_numbers = TRUE)


# ..................................................
# Merge phenotype information with clonotype information
# and analyse their relationship
# ..................................................
cat("<H4>Merging BCR data to Index sort data</H4>")

# Merge the phenotype and clonotype dataframes
phenotype_clonotype_seurat_df = clonotype_seurat_df
phenotype_clonotype_seurat_df$plate.bcid = row.names( phenotype_clonotype_seurat_df)
phenotype_clonotype_seurat_df = merge( phenotype_clonotype_seurat_df, metadata_and_phenotype_df, 
                                       by.x = "plate.bcid", by.y= "plate.bcid", 
                                       all.x = TRUE, all.y = FALSE,
                                       suffixes = c( ".clonotype", ".phenotype"))

cat("<BR>Number of cells in merged data:", nrow( phenotype_clonotype_seurat_df))

# Plot the phenotype information on CXCR3 and CCR6 on the UMAP

ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes_string( x="UMAP_1", y = "UMAP_2", col=PHENOTYPE_CXCR3_INDEXSORT_SOURCE)) +
  theme_classic() + labs( col = "CXCR3") +
ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes( x=UMAP_1, y = UMAP_2, col= CXCR3.phenotype)) +
  theme_classic() + labs( col = "CXCR3") +
ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes_string( x="UMAP_1", y = "UMAP_2", col=PHENOTYPE_CCR6_INDEXSORT_SOURCE)) +
  theme_classic() + labs( col = "CCR6") +
ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes( x=UMAP_1, y = UMAP_2, col= CCR6.phenotype)) +
  theme_classic() + labs( col = "CCR6")

ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes( x=UMAP_1, y = UMAP_2, col=CXCR3xCCR6)) +
  theme_classic() + scale_color_manual( values = c( "red", "blue", "springgreen3", "lavenderblush3") )

# Plot the phenotype information on NA and NP on the UMAP

ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes_string( x="UMAP_1", y = "UMAP_2", col=PHENOTYPE_HA_INDEXSORT_SOURCE)) +
  theme_classic() + labs( col = "HA") +
  ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes( x=UMAP_1, y = UMAP_2, col= HA.phenotype)) +
  theme_classic() + labs( col = "HA") +
  ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes_string( x="UMAP_1", y = "UMAP_2", col=PHENOTYPE_NP_INDEXSORT_SOURCE)) +
  theme_classic() + labs( col = "NP") +
  ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes( x=UMAP_1, y = UMAP_2, col= NP.phenotype)) +
  theme_classic() + labs( col = "NP")

ggplot( phenotype_clonotype_seurat_df) + 
  geom_point( aes( x=UMAP_1, y = UMAP_2, col=HAxNP)) +
  theme_classic() + scale_color_manual( values = c( "red", "blue", "springgreen3", "lavenderblush3") )


# Look at the chi2 statistics between CXCR3xCCR6 phenotypes against cluster
# .........................................................................................................

cat("<H4>CXCR3xCCR6 phenotypes against Seurat clusters</H4>")
cluster_CXCR3xCCR6_wide_df = as.data.frame.matrix( table( phenotype_clonotype_seurat_df$CXCR3xCCR6, phenotype_clonotype_seurat_df$seurat_clusters))
chisq_cluster_CXCR3xCCR6 = chisq.test(  cluster_CXCR3xCCR6_wide_df)

datatable( cluster_CXCR3xCCR6_wide_df, caption = "Number of cells of each CXCR3xCCR6 phenotype for each cluster")
cat("<BR>")
pheatmap( cluster_CXCR3xCCR6_wide_df, cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells of each\nCXCR3xCCR6 phenotype for each Seurat cluster",
          display_numbers = TRUE, number_format = "%d")
pheatmap( floor( chisq_cluster_CXCR3xCCR6$expected), cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Expected number of cells of chi2 test of each\nCXCR3xCCR6 phenotype for each Seurat cluster",
          display_numbers = TRUE, number_format = "%d")
pheatmap( chisq_cluster_CXCR3xCCR6$residuals, cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nnumber of cells of each\nCXCR3xCCR6 phenotype for each Seurat cluster",
          display_numbers = TRUE)

# Look at the chi2 statistics between HAxNP phenotypes against cluster
# .........................................................................................................

cat("<H4>HAxNP phenotypes against Seurat clusters</H4>")
cluster_HAxNP_wide_df = as.data.frame.matrix( table( phenotype_clonotype_seurat_df$HAxNP, phenotype_clonotype_seurat_df$seurat_clusters))
chisq_cluster_HAxNP = chisq.test(  cluster_HAxNP_wide_df)

datatable( cluster_HAxNP_wide_df, caption = "Number of cells of each HAxNP phenotype for each cluster")
cat("<BR>")
pheatmap( cluster_HAxNP_wide_df, cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells of each\nHAxNP phenotype for each Seurat cluster",
          display_numbers = TRUE, number_format = "%d")
pheatmap( floor( chisq_cluster_HAxNP$expected), cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Expected number of cells of chi2 test of each\nHAxNP phenotype for each Seurat cluster",
          display_numbers = TRUE, number_format = "%d")
pheatmap( chisq_cluster_HAxNP$residuals, cellwidth = 30, cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nnumber of cells of each\nHAxNP phenotype for each Seurat cluster",
          display_numbers = TRUE)


# Look at the chi2 statistics between Seurat clusters against large clonotypes
# .........................................................................................................

cat("<H4>Large clonotypes against Seurat clusters</H4>")
largeclonotype_plate_wide_df = as.data.frame.matrix( table( phenotype_clonotype_seurat_df$large.clonotype.similarity.symbol, phenotype_clonotype_seurat_df$Plate))
chisq_largeclonotype_plate = chisq.test(  largeclonotype_plate_wide_df)

datatable( largeclonotype_plate_wide_df, caption = "Number of cells of each plate for each large clonotype")
cat("<BR>")
pheatmap( largeclonotype_plate_wide_df, cellwidth = 15, cellheight = 15, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Number of cells of each\nplate for each large clonotype",
          display_numbers = TRUE, number_format = "%d", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap( floor( chisq_largeclonotype_plate$expected), cellwidth = 15, cellheight = 15, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Expected number of cells of chi2 test of each\nplate for each large clonotype",
          display_numbers = TRUE, number_format = "%d", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap( chisq_largeclonotype_plate$residuals, cellwidth = 15, cellheight = 15, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nnumber of cells of each\nplate for each large clonotype",
         cluster_rows = FALSE, cluster_cols = FALSE)

# Look at the chi2 statistics between CXCR3xCCR6 phenotypes against large clonotypes
# .........................................................................................................

cat("<H4>Large clonotypes against CXCR3xCCR6 phenotypes</H4>")
largeclonotype_CXCR3xCCR6_wide_df = as.data.frame.matrix( table( phenotype_clonotype_seurat_df$CXCR3xCCR6, phenotype_clonotype_seurat_df$large.clonotype.similarity.symbol))
chisq_largeclonotype_CXCR3xCCR6 = chisq.test(  largeclonotype_CXCR3xCCR6_wide_df)

datatable( largeclonotype_CXCR3xCCR6_wide_df, caption = "Number of cells of each CXCR3xCCR6 phenotype for each large clonotype")
cat("<BR>")
pheatmap( largeclonotype_CXCR3xCCR6_wide_df, cellwidth = 30, cellheight = 30, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Number of cells of each\nCXCR3xCCR6 phenotype for each large clonotype",
          display_numbers = TRUE, number_format = "%d", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap( floor( chisq_largeclonotype_CXCR3xCCR6$expected), cellwidth = 30, cellheight = 30, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Expected number of cells of chi2 test of each\nCXCR3xCCR6 phenotype for each large clonotype",
          display_numbers = TRUE, number_format = "%d", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap( chisq_largeclonotype_CXCR3xCCR6$residuals, cellwidth = 30, cellheight = 30, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nnumber of cells of each\nCXCR3xCCR6 phenotype for each large clonotype",
          display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)


# Look at the chi2 statistics between HAxNP phenotypes against large clonotypes
# .........................................................................................................

cat("<H4>Large clonotypes against HAxNP phenotypes</H4>")
largeclonotype_HAxNP_wide_df = as.data.frame.matrix( table( phenotype_clonotype_seurat_df$HAxNP, phenotype_clonotype_seurat_df$large.clonotype.similarity.symbol))
chisq_largeclonotype_HAxNP = chisq.test(  largeclonotype_HAxNP_wide_df)

datatable( largeclonotype_HAxNP_wide_df, caption = "Number of cells of each HAxNP phenotype for each large clonotype")
cat("<BR>")
pheatmap( largeclonotype_HAxNP_wide_df, cellwidth = 30, cellheight = 30, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Number of cells of each\nHAxNP phenotype for each large clonotype",
          display_numbers = TRUE, number_format = "%d", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap( floor( chisq_largeclonotype_HAxNP$expected), cellwidth = 30, cellheight = 30, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Expected number of cells of chi2 test of each\nHAxNP phenotype for each large clonotype",
          display_numbers = TRUE, number_format = "%d", cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap( chisq_largeclonotype_HAxNP$residuals, cellwidth = 30, cellheight = 30, largeclonotype_rows = FALSE, largeclonotype_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nnumber of cells of each\nHAxNP phenotype for each large clonotype",
          display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)

# Look at the chi2 statistics between CXCR3xCCR6 and HAxNP phenotypes against large clonotypes
# .........................................................................................................

cat("<H4>Large clonotypes against CXCR3xCCR6 and HAxNP phenotypes</H4>")
largeclonotype_CXCR3xCCR6_HAxNP_long_df = as.data.frame.table( table( phenotype_clonotype_seurat_df$CXCR3xCCR6, phenotype_clonotype_seurat_df$HAxNP, phenotype_clonotype_seurat_df$large.clonotype.similarity.symbol))
names( largeclonotype_CXCR3xCCR6_HAxNP_long_df) = c( "CXCR3xCCR6", "HAxNP", "Large Clonotype","Frequency")

for( CXCR3xCCR6_value in levels( largeclonotype_CXCR3xCCR6_HAxNP_long_df$CXCR3xCCR6)){
  current_df = largeclonotype_CXCR3xCCR6_HAxNP_long_df[ which( largeclonotype_CXCR3xCCR6_HAxNP_long_df$CXCR3xCCR6 == CXCR3xCCR6_value), ]
  current_df = current_df[ order( current_df$HAxNP),]
  print( htmltools::tagList( datatable( current_df, rownames = FALSE, options = list(dom = 't'))))
}


# Look at the chi2 statistics between CXCR3xCCR6 and HAxNP phenotypes against heavy chain mutation number
# .........................................................................................................

cat("<H4>Heavy chain mutation number against CXCR3xCCR and HAxNP phenotypes</H4>")
mutationqualheavy_CXCR3xCCR6_HAxNP_long_df = as.data.frame.table( table( phenotype_clonotype_seurat_df$CXCR3xCCR6, phenotype_clonotype_seurat_df$HAxNP, phenotype_clonotype_seurat_df$mutation.qual.heavy))
names( mutationqualheavy_CXCR3xCCR6_HAxNP_long_df) = c( "CXCR3xCCR6", "HAxNP", "Number.mutation","Frequency")

# -- Parse the list of CXCR3xCCR6 phenotypes
for( CXCR3xCCR6_value in levels( mutationqualheavy_CXCR3xCCR6_HAxNP_long_df$CXCR3xCCR6)){
  
  # -- Limit the data to the current CXCR3xCCR6 phenotype
  current_df = mutationqualheavy_CXCR3xCCR6_HAxNP_long_df[ which( mutationqualheavy_CXCR3xCCR6_HAxNP_long_df$CXCR3xCCR6 == CXCR3xCCR6_value), ]
  current_df = current_df[ order( current_df$HAxNP),]
  
  # -- Print the data on datatable
  print( htmltools::tagList( datatable( current_df, rownames = FALSE, options = list(dom = 't'))))
  
  # -- Compute the total number of cell for each HAxNP phenotype
  current_total_df = as.data.frame( current_df %>% group_by( HAxNP) %>% summarise( sum( Frequency)))
  row.names( current_total_df) = current_total_df$HAxNP
  
  # -- Compute the density of mutation number for each HAxNP phenotype
  current_df$Density = unlist( apply( current_df, 1, function( row){
    phenotype = as.character( row[ "HAxNP"])
    frequency = as.numeric( as.character( row[ "Frequency"]))
    total = as.numeric( as.character( current_total_df[ phenotype, 2]))
    if( total > 0){
      return( frequency / total)
    }else{
      return( 0)
    }
  }), use.names = FALSE)
  
  # -- Plot the number of mutations frequency and density
  print( ggplot( current_df) + 
    geom_bar( aes( x=Number.mutation, y=Frequency, fill = HAxNP), stat = "identity", position = "dodge")+
    theme_minimal() + ggtitle( "Distribution of heavy chain mutation frequency for", CXCR3xCCR6_value)
  )
  
  print( ggplot( current_df) + 
    geom_bar( aes( x=Number.mutation, y=Density, fill = HAxNP), stat = "identity", position = "dodge")+
    theme_minimal() + ggtitle( "Distribution of heavy chain mutation density for", CXCR3xCCR6_value)
  )
}


# Look at the chi2 statistics between CXCR3xCCR6 and HAxNP phenotypes against light chain mutation number
# .........................................................................................................

cat("<H4>Light chain mutation number against CXCR3xCCR and HAxNP phenotypes</H4>")
mutationquallight_CXCR3xCCR6_HAxNP_long_df = as.data.frame.table( table( phenotype_clonotype_seurat_df$CXCR3xCCR6, phenotype_clonotype_seurat_df$HAxNP, phenotype_clonotype_seurat_df$mutation.qual.light))
names( mutationquallight_CXCR3xCCR6_HAxNP_long_df) = c( "CXCR3xCCR6", "HAxNP", "Number.mutation","Frequency")

# -- Parse the list of CXCR3xCCR6 phenotypes
for( CXCR3xCCR6_value in levels( mutationquallight_CXCR3xCCR6_HAxNP_long_df$CXCR3xCCR6)){
  
  # -- Limit the data to the current CXCR3xCCR6 phenotype
  current_df = mutationquallight_CXCR3xCCR6_HAxNP_long_df[ which( mutationquallight_CXCR3xCCR6_HAxNP_long_df$CXCR3xCCR6 == CXCR3xCCR6_value), ]
  current_df = current_df[ order( current_df$HAxNP),]
  
  # -- Print the data on datatable
  print( htmltools::tagList( datatable( current_df, rownames = FALSE, options = list(dom = 't'))))
  
  # -- Compute the total number of cell for each HAxNP phenotype
  current_total_df = as.data.frame( current_df %>% group_by( HAxNP) %>% summarise( sum( Frequency)))
  row.names( current_total_df) = current_total_df$HAxNP
  
  # -- Compute the density of mutation number for each HAxNP phenotype
  current_df$Density = unlist( apply( current_df, 1, function( row){
    phenotype = as.character( row[ "HAxNP"])
    frequency = as.numeric( as.character( row[ "Frequency"]))
    total = as.numeric( as.character( current_total_df[ phenotype, 2]))
    if( total > 0){
      return( frequency / total)
    }else{
      return( 0)
    }
  }), use.names = FALSE)
  
  # -- Plot the number of mutations frequency and density
  print( ggplot( current_df) + 
           geom_bar( aes( x=Number.mutation, y=Frequency, fill = HAxNP), stat = "identity", position = "dodge")+
           theme_minimal() + ggtitle( "Distribution of light chain mutation frequency for", CXCR3xCCR6_value)
  )
  
  print( ggplot( current_df) + 
           geom_bar( aes( x=Number.mutation, y=Density, fill = HAxNP), stat = "identity", position = "dodge")+
           theme_minimal() + ggtitle( "Distribution of light chain mutation density for", CXCR3xCCR6_value)
  )
}
