# #################################################################################################
# This script aims to prepare some figures of the final article
# #################################################################################################


## @knitr analysis_specific_clonotype

# ##################################################################################################
#
# The NP-spe B cells have a restricted repertoire. If we select cells from the FB5P-seq dataset that
# are IGHV9-3/4+ IGKV5-48+ (100 cells), they are all NP-spe (100%), while in the total dataset, 
# there are 526/1694 NP-spe. Thus, we can use the IGHV9-3/4+ IGKV5-48+ filter to confidently 
# identify NP-spe cells from our previous 10x datasets. In the first 10x dataset, 
# there are 274 such cells (IGHV9-3/4+ IGKV5-48+) among 3733 total cells with paired IGH+IGKL
# sequences (from an old pim matrix file, not the official Lionel dataset J). 
# Many of those putative NP-spe cells contain large clones shared across organs. 
#“Gating” on those cells we may check which transcriptional clusters from Lung, LN and Spleen 
# those cells belong to !
#
# ##################################################################################################

cat("<H3>Dispersion Cells with BCR", paste( paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                         "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or ")), "</H3>")

# List the cells in the filtered set that have also clonotype info
filtered_cells_with_clonotype_info = intersect( merge_umap_cluster_alltissue$cell.id, row.names( clonotype_and_seurat_data_df))

# Merge clonotype and seurat info
clonotype_and_seurat_and_tissue_data = clonotype_and_seurat_data_df
clonotype_and_seurat_and_tissue_data$cell.id = row.names( clonotype_and_seurat_and_tissue_data)
clonotype_and_seurat_and_tissue_data = merge( 
                          clonotype_and_seurat_and_tissue_data[ filtered_cells_with_clonotype_info,], 
                          merge_umap_cluster_alltissue[ which( merge_umap_cluster_alltissue$cell.id %in% filtered_cells_with_clonotype_info), ],
                          by = "cell.id", suffixes = c( ".global", ".tissue"))

# Combine tissue and cluster of tissue in a single feature
clonotype_and_seurat_and_tissue_data$tissue.cluster = paste( clonotype_and_seurat_and_tissue_data$tissue.tissue,
                                                             clonotype_and_seurat_and_tissue_data$cluster, sep = ".")

row.names( clonotype_and_seurat_and_tissue_data) = clonotype_and_seurat_and_tissue_data$cell.id

# Get the cells using the top v-genes
cells_with_vgene_specific = row.names( clonotype_and_seurat_and_tissue_data)[ which( clonotype_and_seurat_and_tissue_data$v.segment.exact.heavy %in% SELECTED_V.SEGMENT.EXACT.HEAVY
                                                                             & clonotype_and_seurat_and_tissue_data$v.segment.exact.light %in% SELECTED_V.SEGMENT.EXACT.LIGHT)]

# Note the cells using the specific BCR in the dataframe for future plot
clonotype_and_seurat_and_tissue_data$specific.BCR = FALSE
clonotype_and_seurat_and_tissue_data[ cells_with_vgene_specific, "specific.BCR"] = TRUE

# ................................................................

cat("<H4>Dispersion over mice of cells with BCR", paste( paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
    "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or ")), "</H4>")

# Dispersion of cells in the mouse
mouse.spe_df = as.data.frame.matrix( t( table( clonotype_and_seurat_and_tissue_data[ cells_with_vgene_specific, "mouse"])))

# Compute the ratio of cell in each Ag-specificity group
nb_cell_mouse.spe = vector()
for( current_mouse.spe in names( mouse.spe_df)){
  nb_cell_mouse.spe = append( nb_cell_mouse.spe, length( which( clonotype_and_seurat_and_tissue_data$mouse == current_mouse.spe)))
}
mouse.spe_ratio_df = signif( 100* mouse.spe_df / nb_cell_mouse.spe,3)


datatable( t(mouse.spe_df),  
           colnames = c("Count"),
           caption= paste( "Dispersion of cells count in mouse using BCR\n", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                                 "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or ")))

datatable( t(mouse.spe_ratio_df),  
           colnames = c("Percentage"),
           caption= paste( "Dispersion of cells percentage in mouse using BCR\n", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                                    "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse="or ")))

# ................................................................

cat("<H4>Dispersion over tissue of cells with BCR", paste( paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                         "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or ")), "</H4>")

# Dispersion of cells in the tissue
tissue.spe_df = as.data.frame.matrix( t( table( clonotype_and_seurat_and_tissue_data[ cells_with_vgene_specific, "tissue.tissue"])))

# Compute the ratio of cell in each Ag-specificity group
nb_cell_tissue.spe = vector()
for( current_tissue.spe in names( tissue.spe_df)){
  nb_cell_tissue.spe = append( nb_cell_tissue.spe, length( which( clonotype_and_seurat_and_tissue_data$tissue.tissue == current_tissue.spe)))
}
tissue.spe_ratio_df = signif( 100* tissue.spe_df / nb_cell_tissue.spe,3)


datatable( t(tissue.spe_df),  
           colnames = c("Count"),
           caption= paste( "Dispersion of cells count in tissue using BCR\n", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                                    "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or ")))

datatable( t(tissue.spe_ratio_df), 
           colnames = c("Percentage"),
           caption= paste( "Dispersion of cells percentage in tissue using BCR\n", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                                          "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or ")))

# ................................................................

cat("<H4>Dispersion over tissue.cluster of cells with BCR", paste( paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                         "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or ")), "</H4>")
# Dispersion of cells in the tissue.cluster
tissue.cluster.spe_df = as.data.frame.matrix( t( table( clonotype_and_seurat_and_tissue_data[ cells_with_vgene_specific, "tissue.cluster"])))

# Compute the ratio of cell in each Ag-specificity group
nb_cell_tissue.cluster.spe = vector()
for( current_tissue.cluster.spe in names( tissue.cluster.spe_df)){
  nb_cell_tissue.cluster.spe = append( nb_cell_tissue.cluster.spe, length( which( clonotype_and_seurat_and_tissue_data$tissue.cluster == current_tissue.cluster.spe)))
}
tissue.cluster.spe_ratio_df = signif( 100* tissue.cluster.spe_df / nb_cell_tissue.cluster.spe,3)

# Display the result in table
datatable( t( tissue.cluster.spe_df), 
           colnames = c("Count"),
           caption= paste( "Dispersion of cells count using BCR\n", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                                             "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or "), "in clusters"))

datatable( t(tissue.cluster.spe_ratio_df), 
           colnames = c( "Percentage"),
           caption= paste( "Dispersion of cells percentage using BCR\n", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                                                           "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or "), "in clusters"))

# Plot the BCR specific cells on UMAP wihtout and with the clusters in color

all_tissue_mice_table_df = data.frame()
for( current_tissue in names( PATH_CLUSTER_MAPPING_FILE)){
  
  # get the tissue specific data
  tissue_df = clonotype_and_seurat_and_tissue_data[ which( clonotype_and_seurat_and_tissue_data$tissue.global == current_tissue), ]
  
  # Plot the BCR specific cells in black on a UMAP with other cells in grey
  print( ggplot() + 
           geom_point( data = tissue_df, aes( x=UMAP_1.tissue, y=UMAP_2.tissue), col="grey") +
           geom_point( data = tissue_df[ which( tissue_df$specific.BCR == TRUE), ], aes( x=UMAP_1.tissue, y=UMAP_2.tissue), col="black") +
           ggtitle( paste( "UMAP enbedding of", current_tissue, "with cells using\nBCR", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                           "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or "))) +
           theme_classic()
  )
  
  ggsave( file=file.path( PATH_ANALYSIS_OUTPUT, paste0( "UMAP_", current_tissue, "_BCRspecific_", paste( c( SELECTED_V.SEGMENT.EXACT.HEAVY, SELECTED_V.SEGMENT.EXACT.LIGHT), collapse="_"), ".svg" )))
  
  # Plot the BCR specific cells in black on a UMAP with other cells in cluster color
  print( ggplot() + 
           geom_point( data = tissue_df, aes( x=UMAP_1.tissue, y=UMAP_2.tissue, col=cluster)) +
           geom_point( data = tissue_df[ which( tissue_df$specific.BCR == TRUE), ], aes( x=UMAP_1.tissue, y=UMAP_2.tissue), col="black", size = 2) +
           ggtitle( paste( "UMAP enbedding of", current_tissue, "with cells using\nBCR", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                           "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or "))) +
           theme_classic()
  )
  
  ggsave( file=file.path( PATH_ANALYSIS_OUTPUT, paste0( "UMAP_", current_tissue, "_BCRspecific_", paste( c( SELECTED_V.SEGMENT.EXACT.HEAVY, SELECTED_V.SEGMENT.EXACT.LIGHT), collapse="_"), "_WithClusters.svg" )))
  
  # Look at the expression of marker genes between cells using specific BCR and others
  
  # -- Get the expression data from tissue
  current_expression_tissue_df = data.frame( t( variable_genes_expression_list[[ current_tissue]]))
  
  # -- Limit the data to the marker genes
  current_expression_tissue_df = current_expression_tissue_df[ , MARKER_GENES_LIST]
  
  # -- Indicates which cell use Specific BCR
  current_expression_tissue_df$bcr.specific = "NotBCRSpe"
  current_expression_tissue_df[ intersect( row.names( current_expression_tissue_df), cells_with_vgene_specific), "bcr.specific"] = "BCRSpe"
  

  for( current_mouse in levels( clonotype_and_seurat_and_tissue_data$mouse)){
    
    # -- Get the cell of the current mouse
    current_mouse_cells = clonotype_and_seurat_and_tissue_data[ which( clonotype_and_seurat_and_tissue_data$mouse == current_mouse), "cell.id"]
    
    # -- Limitate the expression data by mouse
    current_mouse_expression_tissue_df = current_expression_tissue_df[ intersect( row.names( current_expression_tissue_df), current_mouse_cells), ]
    
    # -- Reshape the dataframe to long format for ggplot
    current_mouse_expression_tissue_long_df = reshape2::melt( current_mouse_expression_tissue_df, varying = MARKER_GENES_LIST)
    
    # -- Look at the number of cells having a positive expression for each genes in each BCR group
    current_mouse_expression_tissue_long_df$positive.expression = "ZeroExp"
    current_mouse_expression_tissue_long_df[ which( current_mouse_expression_tissue_long_df$value >0), "positive.expression"] = "PositiveExp"
    mouse_table_df = data.frame( table( current_mouse_expression_tissue_long_df$variable, current_mouse_expression_tissue_long_df$bcr.specific, current_mouse_expression_tissue_long_df$positive.expression,
                        dnn = c( "Gene", "BCR specific", "Positive Expression")))

    # -- Accumulate the info in a global table
    mouse_table_df$mouse = current_mouse
    mouse_table_df$tissue = current_tissue
    all_tissue_mice_table_df = rbind( all_tissue_mice_table_df, mouse_table_df)

  }

}

all_tissue_mice_table_df = all_tissue_mice_table_df[ , c( "Gene", "tissue", "mouse", "BCR.specific", "Positive.Expression", "Freq")]

# Plot the numner/percentage of cells in groups defined by Cxcr3 expression and BCR specificity
all_tissue_mice_table_df_Cxcr3 = all_tissue_mice_table_df[ which( all_tissue_mice_table_df$Gene == "Cxcr3"),]
all_tissue_mice_table_df_Cxcr3$tissue = factor( all_tissue_mice_table_df_Cxcr3$tissue)
all_tissue_mice_table_df_Cxcr3$mouse = factor( all_tissue_mice_table_df_Cxcr3$mouse)

all_tissue_mice_table_df_Cxcr3_pct = all_tissue_mice_table_df_Cxcr3 %>%
  group_by(mouse, tissue) %>%
  mutate(pct = 100*Freq / sum( Freq))

ggplot( all_tissue_mice_table_df_Cxcr3) + 
  geom_bar( aes( x=interaction( mouse, tissue), fill= interaction( Positive.Expression, BCR.specific), y=Freq), stat="identity") +
  theme_classic() + theme( axis.text.x = element_text( angle = 45, hjust = 1)) +
  labs( fill = "Group", y="Count", x="Mouse/Tissue") +
  ggtitle( "Cell count of Cxcr3 expression/BCR specificity groups")

ggplot( all_tissue_mice_table_df_Cxcr3_pct) + 
  geom_bar( aes( x=interaction( tissue, mouse), fill= interaction( Positive.Expression, BCR.specific), y=pct), stat="identity") +
  theme_classic() + theme( axis.text.x = element_text( angle = 45, hjust = 1))+
  labs( fill = "Group", y="Percentage", x="Mouse/Tissue") +
  ggtitle( "Cell percentage of Cxcr3 expression/BCR specificity groups")

print( htmltools::tagList( datatable( all_tissue_mice_table_df_Cxcr3_pct, caption = "Cell count and percentage of Cxcr3 expression/BCR specificity groups")))

for( current_tissue in names( PATH_CLUSTER_MAPPING_FILE)){
  
  cat( "<H5>Testing if Cxcr3+ cells are enriched in BCR Specific group in tissue", current_tissue, "</H5>")
  cat( "<H6>with GLM/Anova test</H6>")
  current_df = all_tissue_mice_table_df_Cxcr3[ which( all_tissue_mice_table_df_Cxcr3$tissue == current_tissue),]
  model1 = glm( formula = Freq ~ BCR.specific + Positive.Expression, family = "poisson", 
                data = current_df)
  model2 = glm( formula = Freq ~ BCR.specific + Positive.Expression + BCR.specific:Positive.Expression, family = "poisson", 
                data = current_df)
  
  cat("<BR><b>Model 1 : Freq ~ BCR.specific + Positive.Expression</b><BR>")
  print( xtable::xtable( summary( model1), display = c( "d", "f", "f", "f", "e")), type = "html")
  
  cat("<BR><b>Model 2 : Freq ~ BCR.specific + Positive.Expression + BCR.specific:Positive.Expression</b><BR>")
  print( xtable::xtable( summary( model2), display = c( "d", "f", "f", "f", "e")), type = "html")
  
  cat("<BR><b>Anova between the two models with LRT test</b><BR>")
  print( xtable::xtable( anova( model1, model2, test="LRT"), display = c( "d", "d", "f", "d", "f", "e")), type = "html")
  
  cat( "<H6>with  Cochran-Mantel-Haenszel (CMH) test</H6><BR>")
  
  cmh_test = vcdExtra::CMHtest( Freq~ BCR.specific + Positive.Expression + mouse, current_df, overall = TRUE)
  
  for (name in names(cmh_test)) {
    cat( "<BR><b>Mouse:", name, "</b><BR>")
    print( xtable::xtable( cmh_test[[name]]$table, display = c( "d", "f", "d", "e")), type = "html")
    cat("<BR>")
  }
}



# ................................................................

cat("<H4>Estimation of dispersion of NP+ cells over tissue.cluster</H4>")

# Compute the confusion matrix of putative NP+ cells in tissue.cluster
# Considering cells using BCR combination IGHV9-3/IGHV9-4 and IGKV5-48 represent 20% of cells in the custom dataset, we can estimate the number of NP+ cells in the current dataset by multiplying the number ofcells with CR combination IGHV9-3/IGHV9-4 and IGKV5-48 by 5 and make a confusion matrix with the other cells against the tissue.clusters.
# This operation considers that the dispersion of NP+ cells of other BCR combinations alon gthe tissue.cluster is similar to the dispersion of the population of BCR IGHV9-3/IGHV9-4 and IGKV5-48

cat("<BR>Considering cells using BCR combination IGHV9-3/IGHV9-4 and IGKV5-48 represent 20% of cells in the custom dataset, we can estimate the number of NP+ cells in the current dataset by multiplying by 5 the number ofcells with BCR combination IGHV9-3/IGHV9-4 and IGKV5-48 and make a confusion matrix with the other cells against the tissue.clusters.")
cat ("This operation considers that the dispersion of NP+ cells of other BCR combinations along the tissue.cluster is similar to the dispersion of the population of BCR IGHV9-3/IGHV9-4 and IGKV5-48")

# Compute the putative number of NP+ cells and of other cells
tissue.cluster.spe_confusion_df = 5* tissue.cluster.spe_df
tissue.cluster.spe_confusion_df = rbind( tissue.cluster.spe_confusion_df, nb_cell_tissue.cluster.spe)
tissue.cluster.spe_confusion_df[ 2, ] = tissue.cluster.spe_confusion_df[ 2, ] - tissue.cluster.spe_confusion_df[ 1, ]
row.names( tissue.cluster.spe_confusion_df) = c( "Putative NP+", "Putative Non NP+")

# Show the result in table
datatable( t(tissue.cluster.spe_confusion_df), 
           caption= paste( "Dispersion of cells using BCR specific V-genes versus non-NP+ cells in clusters "))

# Compute a chi2 test on the confusion table
chisq_tissue.cluster.spe_confusion_df = chisq.test(  tissue.cluster.spe_confusion_df)

# Show the Pearson residuals of the chi2 test as a heatmap
pheatmap::pheatmap( t( chisq_tissue.cluster.spe_confusion_df$residuals), 
                    cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
                    main = paste( "Pearson residuals of chi2 test on \ncells count in cluster using BCR\n", paste( SELECTED_V.SEGMENT.EXACT.HEAVY, collapse=" or "),
                                   "and", paste( SELECTED_V.SEGMENT.EXACT.LIGHT, collapse=" or ")))


