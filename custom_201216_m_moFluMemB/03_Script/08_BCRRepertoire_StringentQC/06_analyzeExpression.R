# ##############################################################################
# This script aims at analyzing the FACS phenotypes in relation with the
# gene expression
# ##############################################################################

## @knitr analyze_phenotype_to_expression

# Plot a heatmap to see the expression of genes of interest on population of cell CXCR3+/CCR6+, CXCR3-/CCR6-, CXRCR3-/CCR6+
# #########################################################################################################################

# -- Get the cells ID for each phenotype groups
CXCR3p.CCR6p_cells = gsub( "201216_m_moFluMemB_", "", row.names( metadata_and_phenotype_df)[ which( metadata_and_phenotype_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "PosPos"])])
CXCR3n.CCR6n_cells = gsub( "201216_m_moFluMemB_", "", row.names( metadata_and_phenotype_df)[ which( metadata_and_phenotype_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "NegNeg"])])
CXCR3n.CCR6p_cells = gsub( "201216_m_moFluMemB_", "", row.names( metadata_and_phenotype_df)[ which( metadata_and_phenotype_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "NegPos"])])

# -- Get the expression of genes of interest for each phenotype group and make the mean over the cells
CXCR3p.CCR6p_expression_heatmap_df = as.data.frame( GetAssayData(object = custom.rna.seurat, assay = "RNA")[ GENES_OF_INTEREST, 
                                                                                 c( CXCR3p.CCR6p_cells)])
CXCR3p.CCR6p_expression_heatmap_df$CXCR3p.CCR6p = apply( CXCR3p.CCR6p_expression_heatmap_df, 1, mean)


CXCR3n.CCR6n_expression_heatmap_df = as.data.frame( GetAssayData(object = custom.rna.seurat, assay = "RNA")[ GENES_OF_INTEREST, 
                                                                                 c(CXCR3n.CCR6n_cells)])
CXCR3n.CCR6n_expression_heatmap_df$CXCR3n.CCR6n = apply( CXCR3n.CCR6n_expression_heatmap_df, 1, mean)


CXCR3n.CCR6p_expression_heatmap_df = as.data.frame( GetAssayData(object = custom.rna.seurat, assay = "RNA")[ GENES_OF_INTEREST, 
                                                                                 c( CXCR3n.CCR6p_cells)])
CXCR3n.CCR6p_expression_heatmap_df$CXCR3n.CCR6p = apply( CXCR3n.CCR6p_expression_heatmap_df, 1, mean)

# Group the expression mean information fo each phenotype group
expression_heatmap_df = data.frame( CXCR3p.CCR6p = CXCR3p.CCR6p_expression_heatmap_df$CXCR3p.CCR6p,
                                    CXCR3n.CCR6n = CXCR3n.CCR6n_expression_heatmap_df$CXCR3n.CCR6n,
                                    CXCR3n.CCR6p = CXCR3n.CCR6p_expression_heatmap_df$CXCR3n.CCR6p)
row.names( expression_heatmap_df) = GENES_OF_INTEREST

# Plot a heatmap of the mean expression of gene of interest for each phenotype group
cat("<BR>Heatmap of Genes of interest Mean Expression  on phenotype groups CXCR3xCCR6")
print( pheatmap::pheatmap( expression_heatmap_df,
                           cluster_rows = TRUE, cluster_cols= FALSE,
                           cellwidth = 10, cellheight = 10,
                           color = viridis::viridis(20, option = "B"),
                           scale = "row", labels_col = c("CXCR3+ CCR6+", "CXCR3- CCR6-", "CXCR3- CCR6+")
))

svg( file = file.path( PATH_ANALYSIS_OUTPUT, "Heatmap_InterestGenes_MeanExpression_phenotypeCXCR3xCCR6.svg") )
  print( pheatmap::pheatmap( expression_heatmap_df,
                      cluster_rows = TRUE, cluster_cols= FALSE,
                      cellwidth = 10, cellheight = 10,
                      color = viridis::viridis(20, option = "B"),
                      scale = "row", labels_col = c("CXCR3+ CCR6+", "CXCR3- CCR6-", "CXCR3- CCR6+")
                      )
  )
dev.off()


# Plot a heatmap to see the expression of genes of interest on population of cell CXCR3+/CCR6+, CXCR3-/CCR6-, CXRCR3-/CCR6+
# with the details on the mice
# #########################################################################################################################

# -- For each mouse, get the expression of genes of interest for each phenotype group and make the mean over the cells
allmouse_expression_heatmap_df = data.frame()
for( current_mouse in levels( metadata_and_phenotype_df$Mouse_ID)){
  
  # -- Get the cells ID for each phenotype groups
  CXCR3p.CCR6p_cells = gsub( "201216_m_moFluMemB_", "", row.names( metadata_and_phenotype_df)[ which( metadata_and_phenotype_df$Mouse_ID == current_mouse & metadata_and_phenotype_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "PosPos"])])
  CXCR3n.CCR6n_cells = gsub( "201216_m_moFluMemB_", "", row.names( metadata_and_phenotype_df)[ which( metadata_and_phenotype_df$Mouse_ID == current_mouse & metadata_and_phenotype_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "NegNeg"])])
  CXCR3n.CCR6p_cells = gsub( "201216_m_moFluMemB_", "", row.names( metadata_and_phenotype_df)[ which( metadata_and_phenotype_df$Mouse_ID == current_mouse & metadata_and_phenotype_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "NegPos"])])
  
  # -- Get the expression of genes of interest for each phenotype group and make the mean over the cells
  mouse_CXCR3p.CCR6p_expression_heatmap_df = as.data.frame( GetAssayData(object = custom.rna.seurat, assay = "RNA")[ GENES_OF_INTEREST, 
                                                                                                               c( CXCR3p.CCR6p_cells)])
  mouse_CXCR3p.CCR6p_expression_heatmap_df$CXCR3p.CCR6p = apply( mouse_CXCR3p.CCR6p_expression_heatmap_df, 1, mean)
  
  mouse_CXCR3n.CCR6n_expression_heatmap_df = as.data.frame( GetAssayData(object = custom.rna.seurat, assay = "RNA")[ GENES_OF_INTEREST, 
                                                                                                               c(CXCR3n.CCR6n_cells)])
  mouse_CXCR3n.CCR6n_expression_heatmap_df$CXCR3n.CCR6n = apply( mouse_CXCR3n.CCR6n_expression_heatmap_df, 1, mean)
  
  mouse_CXCR3n.CCR6p_expression_heatmap_df = as.data.frame( GetAssayData(object = custom.rna.seurat, assay = "RNA")[ GENES_OF_INTEREST, 
                                                                                                               c( CXCR3n.CCR6p_cells)])
  mouse_CXCR3n.CCR6p_expression_heatmap_df$CXCR3n.CCR6p = apply( mouse_CXCR3n.CCR6p_expression_heatmap_df, 1, mean)
  
  # Group the expression mean information fo each phenotype group
  mouse_expression_heatmap_df = data.frame( CXCR3p.CCR6p = mouse_CXCR3p.CCR6p_expression_heatmap_df$CXCR3p.CCR6p,
                                            CXCR3n.CCR6n = mouse_CXCR3n.CCR6n_expression_heatmap_df$CXCR3n.CCR6n,
                                            CXCR3n.CCR6p = mouse_CXCR3n.CCR6p_expression_heatmap_df$CXCR3n.CCR6p)
  row.names( mouse_expression_heatmap_df) = GENES_OF_INTEREST
  names( mouse_expression_heatmap_df) = paste0( names( mouse_expression_heatmap_df), "_", gsub( "moLung", "", current_mouse))
  
  # Accumulate over the mice 
  if( nrow( allmouse_expression_heatmap_df) == 0){
    allmouse_expression_heatmap_df = mouse_expression_heatmap_df
  }else{
    allmouse_expression_heatmap_df = cbind( allmouse_expression_heatmap_df, mouse_expression_heatmap_df)
  }
}

# Order the columns to group by phenotype
allmouse_expression_heatmap_df = allmouse_expression_heatmap_df[ , c(1,4,7,2,5,8,3,6,9)]

# Plot a heatmap of the mean expression of gene of interest for each phenotype group
cat("<BR>Heatmap of Genes of interest Mean Expression  on phenotype groups CXCR3xCCR6 with mice details")
dev.off()
print( pheatmap::pheatmap( allmouse_expression_heatmap_df,
                           cluster_rows = TRUE, cluster_cols= TRUE,
                           cellwidth = 10, cellheight = 10,
                           color = viridis::viridis(20, option = "B"),
                           scale = "row" )
)


svg( file = file.path( PATH_ANALYSIS_OUTPUT, "Heatmap_InterestGenes_MeanExpression_phenotypeCXCR3xCCR6_micedetail_notclustered.svg") )
print( pheatmap::pheatmap( allmouse_expression_heatmap_df,
                           cluster_rows = TRUE, cluster_cols= FALSE,
                           cellwidth = 10, cellheight = 10,
                           color = viridis::viridis(20, option = "B"),
                           scale = "row" )
)
dev.off()


svg( file = file.path( PATH_ANALYSIS_OUTPUT, "Heatmap_InterestGenes_MeanExpression_phenotypeCXCR3xCCR6_micedetail_clustered.svg") )
  print( pheatmap::pheatmap( allmouse_expression_heatmap_df,
                           cluster_rows = TRUE, cluster_cols= TRUE,
                           cellwidth = 10, cellheight = 10,
                           color = viridis::viridis(20, option = "B"),
                           scale = "row" )
        )
dev.off()

