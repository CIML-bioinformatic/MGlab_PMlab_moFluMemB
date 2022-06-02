# ##############################################################################
# This script aims at analyzing specific clonotypes idenitifies by their
# V gene usage
# ##############################################################################

## @knitr analyze_specific_clonotypes

# Order the tables in decreasing V-gene usage
V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF = V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF[ order( V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF$total, decreasing = TRUE), ]
V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF = V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF[ order( V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF$total, decreasing = TRUE), ]

print( htmltools::tagList( datatable( V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF[ 1:20,], caption="Top 20 Heavy-chain gene usage")))
print( htmltools::tagList( datatable( V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF[ 1:20,], caption="Top 20 Light-chain gene usage")))


# ..............................................................................
# ..............................................................................

cat("<H4>Analysis of AG specificity of cells using some specific V-genes</H4>")

# Get the v-gene names of the top usage
best_v.segment.exact.heavy = c( "IGHV9-3" )
best_v.segment.exact.light = c( "IGKV5-48")

# Get the cells using the top v-genes
cells_with_vgene_specific = phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$v.segment.exact.heavy %in% best_v.segment.exact.heavy
                                                                             & phenotype_clonotype_seurat_df$v.segment.exact.light %in% best_v.segment.exact.light)]
# Dispersion of cells in the Ag specifity groups
ag.spe_df = as.data.frame.matrix( t( table( phenotype_clonotype_seurat_df[ which( phenotype_clonotype_seurat_df$plate.bcid %in% cells_with_vgene_specific), "HAxNP"])))

# Compute the raio of cell in each Ag-specificity group
nb_cell_ag.spe = vector()
for( current_ag.spe in names( ag.spe_df)){
  nb_cell_ag.spe = append( nb_cell_ag.spe, length( phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$HAxNP == current_ag.spe)]))
}
ag.spe_ratio_df = signif( 100* ag.spe_df / nb_cell_ag.spe,3)

print( htmltools::tagList( datatable( ag.spe_df, caption= paste( "Ag specificity of cells using BCR\n", paste( best_v.segment.exact.heavy, collapse="or"),
                                                                 "and", paste( best_v.segment.exact.light, collapse=" or")))))

# ..............................................................................
# ..............................................................................

# Get the v-gene names of the top usage
best_v.segment.exact.heavy = c( "IGHV9-4" )
best_v.segment.exact.light = c( "IGKV5-48")

# Get the cells using the top v-genes
cells_with_vgene_specific = phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$v.segment.exact.heavy %in% best_v.segment.exact.heavy
                                                                             & phenotype_clonotype_seurat_df$v.segment.exact.light %in% best_v.segment.exact.light)]
# Dispersion of cells in the Ag specifity groups
ag.spe_df = as.data.frame.matrix( t( table( phenotype_clonotype_seurat_df[ which( phenotype_clonotype_seurat_df$plate.bcid %in% cells_with_vgene_specific), "HAxNP"])))

# Compute the raio of cell in each Ag-specificity group
nb_cell_ag.spe = vector()
for( current_ag.spe in names( ag.spe_df)){
  nb_cell_ag.spe = append( nb_cell_ag.spe, length( phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$HAxNP == current_ag.spe)]))
}
ag.spe_ratio_df = signif( 100* ag.spe_df / nb_cell_ag.spe,3)

print( htmltools::tagList( datatable( ag.spe_df, caption= paste( "Ag specificity of cells using BCR\n", paste( best_v.segment.exact.heavy, collapse="or"),
                                                                 "and", paste( best_v.segment.exact.light, collapse=" or")))))

# ..............................................................................
# ..............................................................................

# Get the v-gene names of the top usage
best_v.segment.exact.heavy = c( "IGHV9-3", "IGHV9-4" )
best_v.segment.exact.light = c( "IGKV5-48", "IGKV5-43", "IGKV5-37")

# Get the cells using the top v-genes
cells_with_vgene_specific = phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$v.segment.exact.heavy %in% best_v.segment.exact.heavy
                                                                             & phenotype_clonotype_seurat_df$v.segment.exact.light %in% best_v.segment.exact.light)]
# Dispersion of cells in the Ag specifity groups
ag.spe_df = as.data.frame.matrix( t( table( phenotype_clonotype_seurat_df[ which( phenotype_clonotype_seurat_df$plate.bcid %in% cells_with_vgene_specific), "HAxNP"])))

# Compute the raio of cell in each Ag-specificity group
nb_cell_ag.spe = vector()
for( current_ag.spe in names( ag.spe_df)){
  nb_cell_ag.spe = append( nb_cell_ag.spe, length( phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$HAxNP == current_ag.spe)]))
}
ag.spe_ratio_df = signif( 100* ag.spe_df / nb_cell_ag.spe,3)

print( htmltools::tagList( datatable( ag.spe_df, caption= paste( "Ag specificity of cells using BCR", paste( best_v.segment.exact.heavy, collapse=" or "),
                                                                 "and", paste( best_v.segment.exact.light, collapse=" or")))))

# ..............................................................................
# ..............................................................................

cat("<H4>Analysis of AG specificity of cells useng the top used V-genes</H4>")

# Loop over an increasing number of top v-gene usage
specificity_evolution_vgenetop_df = data.frame()
specificity_evolution_vgenetop_ratio_df = data.frame()
for( top_vgene_usage in 1:20){

  # Get the v-gene names of the top usage
  best_v.segment.exact.heavy = c( row.names( V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF)[c(1,3)])
  best_v.segment.exact.light = c( row.names( V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF)[1:1])
  best_v.segment.exact.light = c( "IGKV5-48", "IGKV5-43", "IGKV5-37")

  # Get the cells using the top v-genes
  cells_with_vgene_specific = phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$v.segment.exact.heavy %in% best_v.segment.exact.heavy
                                                                            & phenotype_clonotype_seurat_df$v.segment.exact.light %in% best_v.segment.exact.light)]
  # Dispersion of cells in the Ag specifity groups
  ag.spe_df = as.data.frame.matrix( t( table( phenotype_clonotype_seurat_df[ which( phenotype_clonotype_seurat_df$plate.bcid %in% cells_with_vgene_specific), "HAxNP"])))
  
  # Compute the raio of cell in each Ag-specificity group
  nb_cell_ag.spe = vector()
  for( current_ag.spe in names( ag.spe_df)){
    nb_cell_ag.spe = append( nb_cell_ag.spe, length( phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$HAxNP == current_ag.spe)]))
  }
  ag.spe_ratio_df = signif( 100* ag.spe_df / nb_cell_ag.spe,3)
  
  
  # look at the Ag specificity of the cells ysing the top v-genes and accumulate it in a single dataframe
 specificity_evolution_vgenetop_df = rbind( specificity_evolution_vgenetop_df,   
                                            ag.spe_df)
 specificity_evolution_vgenetop_ratio_df = rbind( specificity_evolution_vgenetop_ratio_df,   
                                            ag.spe_ratio_df)
 
}
# Rename the row of the global dataframe 
row.names( specificity_evolution_vgenetop_df) = seq( 1, nrow( specificity_evolution_vgenetop_df))
row.names( specificity_evolution_vgenetop_ratio_df) = seq( 1, nrow( specificity_evolution_vgenetop_ratio_df))

# Normalize the row in percentages
specificity_evolution_vgenetop_df_pct = t(apply(specificity_evolution_vgenetop_df, 1, function(x)(signif( 100*x/sum(x), 3))))

# Print the count table and the percentage table
print( htmltools::tagList( datatable( specificity_evolution_vgenetop_df,
                                      caption = "Number of cells in fonction of the number of top V-gene usage\nin AG affinity group")))

print( htmltools::tagList( datatable( specificity_evolution_vgenetop_df_pct,
                                      caption = "Percentage of cells in fonction of the number of top V-gene usage\nin AG affinity group")))

print( htmltools::tagList( datatable( specificity_evolution_vgenetop_ratio_df,
                                      caption = "Percentage of cells in AG affinity group\nin fonction of the number of top V-gene usage")))


pheatmap( t( specificity_evolution_vgenetop_df), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells in fonction of the number of top V-gene usage\nin AG affinity group",
          display_numbers = TRUE, number_format = "%d")

pheatmap( t( specificity_evolution_vgenetop_df_pct), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Percentage of cells in fonction of the number of top V-gene usage\nin AG affinity group")

pheatmap( t( specificity_evolution_vgenetop_ratio_df), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Percentage of cells in AG affinity group\nin fonction of the number of top V-gene usage")

# Look at the chi2 statistic on the whole table
chisq_specificity_evolution_vgenetop_df = chisq.test(  specificity_evolution_vgenetop_df)
cat("<BR>")

pheatmap( t( chisq_specificity_evolution_vgenetop_df$residuals), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nNumber of cells in fonction of the number of top V-gene usage\nin AG affinity group")


# .......................................................................................
# .......................................................................................

cat("<H4>Analysis of Mouse specificity of cells using the top used V-genes</H4>")

# Loop over an increasing number of top v-gene usage
mouse_evolution_vgenetop_df = data.frame()
mouse_evolution_vgenetop_ratio_df = data.frame()
for( top_vgene_usage in 1:20){
  
  # Get the v-gene names of the top usage
  best_v.segment.exact.heavy = c( row.names( V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF)[1:top_vgene_usage])
  best_v.segment.exact.light = c( row.names( V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF)[1:top_vgene_usage])
  
  # Get the cells using the top v-genes
  cells_with_vgene_specific = phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$v.segment.exact.heavy %in% best_v.segment.exact.heavy
                                                                               & phenotype_clonotype_seurat_df$v.segment.exact.light %in% best_v.segment.exact.light)]
  # Dispersion of cells in the mice
  mouse.spe_df = as.data.frame.matrix( t( table( phenotype_clonotype_seurat_df[ which( phenotype_clonotype_seurat_df$plate.bcid %in% cells_with_vgene_specific), "mouse"])))
  
  # Compute the raio of cell in each Ag-specificity group
  nb_cell_mouse.spe = vector()
  for( current_mouse.spe in names( mouse.spe_df)){
    nb_cell_mouse.spe = append( nb_cell_mouse.spe, length( phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$mouse == current_mouse.spe)]))
  }
  mouse.spe_ratio_df = signif( 100* mouse.spe_df / nb_cell_mouse.spe,3)
  
  
  # look at the Ag specificity of the cells ysing the top v-genes and accumulate it in a single dataframe
  mouse_evolution_vgenetop_df = rbind( mouse_evolution_vgenetop_df,   
                                       mouse.spe_df)
  mouse_evolution_vgenetop_ratio_df = rbind( mouse_evolution_vgenetop_ratio_df,   
                                             mouse.spe_ratio_df)
  
}
# Rename the row of the global dataframe 
row.names( mouse_evolution_vgenetop_df) = seq( 1, nrow( mouse_evolution_vgenetop_df))
row.names( mouse_evolution_vgenetop_ratio_df) = seq( 1, nrow( mouse_evolution_vgenetop_ratio_df))

# Normalize the row in percentages
mouse_evolution_vgenetop_df_pct = t(apply(mouse_evolution_vgenetop_df, 1, function(x)(signif( 100*x/sum(x), 3))))

# Print the count table and the percentage table
print( htmltools::tagList( datatable( mouse_evolution_vgenetop_df,
                                      caption = "Number of cells in fonction of the number of top V-gene usage\nin mice")))

print( htmltools::tagList( datatable( mouse_evolution_vgenetop_df_pct,
                                      caption = "Percentage of cells in fonction of the number of top V-gene usage\nin mice")))

print( htmltools::tagList( datatable( mouse_evolution_vgenetop_ratio_df,
                                      caption = "Percentage of cells in mice\nin fonction of the number of top V-gene usage")))


pheatmap( t( mouse_evolution_vgenetop_df), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells in fonction of the number of top V-gene usage\nin mice",
          display_numbers = TRUE, number_format = "%d")

pheatmap( t( mouse_evolution_vgenetop_df_pct), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Percentage of cells in fonction of the number of top V-gene usage\nin mice")

pheatmap( t( mouse_evolution_vgenetop_ratio_df), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Percentage of cells in mice\nin fonction of the number of top V-gene usage")

# Look at the chi2 statistic on the whole table
chisq_mouse_evolution_vgenetop_df = chisq.test(  mouse_evolution_vgenetop_df)
cat("<BR>")

pheatmap( t( chisq_mouse_evolution_vgenetop_df$residuals), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nNumber of cells in fonction of the number of top V-gene usage\nin mice")


# .......................................................................................
# .......................................................................................

cat("<H4>Analysis of Clonotype size category of cells using the top used V-genes</H4>")

# Loop over an increasing number of top v-gene usage
clonotype.size.category_evolution_vgenetop_df = data.frame()
clonotype.size.category_evolution_vgenetop_ratio_df = data.frame()
for( top_vgene_usage in 1:20){
  
  # Get the v-gene names of the top usage
  best_v.segment.exact.heavy = c( row.names( V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF)[1:top_vgene_usage])
  best_v.segment.exact.light = c( row.names( V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF)[1:top_vgene_usage])
  
  # Get the cells using the top v-genes
  cells_with_vgene_specific = phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$v.segment.exact.heavy %in% best_v.segment.exact.heavy
                                                                               & phenotype_clonotype_seurat_df$v.segment.exact.light %in% best_v.segment.exact.light)]
  # Dispersion of cells in the mice
  clonotype.size.category.spe_df = as.data.frame.matrix( t( table( phenotype_clonotype_seurat_df[ which( phenotype_clonotype_seurat_df$plate.bcid %in% cells_with_vgene_specific), "clonotype.size.category"])))
  
  # Compute the raio of cell in each Ag-specificity group
  nb_cell_clonotype.size.category.spe = vector()
  for( current_clonotype.size.category.spe in names( clonotype.size.category.spe_df)){
    nb_cell_clonotype.size.category.spe = append( nb_cell_clonotype.size.category.spe, length( phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$clonotype.size.category == current_clonotype.size.category.spe)]))
  }
  clonotype.size.category.spe_ratio_df = signif( 100* clonotype.size.category.spe_df / nb_cell_clonotype.size.category.spe,3)
  
  
  # look at the Ag specificity of the cells ysing the top v-genes and accumulate it in a single dataframe
  clonotype.size.category_evolution_vgenetop_df = rbind( clonotype.size.category_evolution_vgenetop_df,   
                                       clonotype.size.category.spe_df)
  clonotype.size.category_evolution_vgenetop_ratio_df = rbind( clonotype.size.category_evolution_vgenetop_ratio_df,   
                                             clonotype.size.category.spe_ratio_df)
  
}
# Rename the row of the global dataframe 
row.names( clonotype.size.category_evolution_vgenetop_df) = seq( 1, nrow( clonotype.size.category_evolution_vgenetop_df))
row.names( clonotype.size.category_evolution_vgenetop_ratio_df) = seq( 1, nrow( clonotype.size.category_evolution_vgenetop_ratio_df))

# Normalize the row in percentages
clonotype.size.category_evolution_vgenetop_df_pct = t(apply(clonotype.size.category_evolution_vgenetop_df, 1, function(x)(signif( 100*x/sum(x), 3))))

# Print the count table and the percentage table
print( htmltools::tagList( datatable( clonotype.size.category_evolution_vgenetop_df,
                                      caption = "Number of cells in fonction of the number of top V-gene usage\nin Clonotype size")))

print( htmltools::tagList( datatable( clonotype.size.category_evolution_vgenetop_df_pct,
                                      caption = "Percentage of cells in fonction of the number of top V-gene usage\nin Clonotype size")))

print( htmltools::tagList( datatable( clonotype.size.category_evolution_vgenetop_ratio_df,
                                      caption = "Percentage of cells in Clonotype size\nin fonction of the number of top V-gene usage")))


pheatmap( t( clonotype.size.category_evolution_vgenetop_df), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells in fonction of the number of top V-gene usage\nin Clonotype size",
          display_numbers = TRUE, number_format = "%d")

pheatmap( t( clonotype.size.category_evolution_vgenetop_df_pct), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Percentage of cells in fonction of the number of top V-gene usage\nin Clonotype size")

pheatmap( t( clonotype.size.category_evolution_vgenetop_ratio_df), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Percentage of cells in Clonotype size\nin fonction of the number of top V-gene usage")

# Look at the chi2 statistic on the whole table
chisq_clonotype.size.category_evolution_vgenetop_df = chisq.test(  clonotype.size.category_evolution_vgenetop_df)
cat("<BR>")

pheatmap( t( chisq_clonotype.size.category_evolution_vgenetop_df$residuals), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nNumber of cells in fonction of the number of top V-gene usage\nin Clonotype size")


# .......................................................................................
# .......................................................................................

cat("<H4>Analysis of Cluster of cells using the top used V-genes</H4>")

# Loop over an increasing number of top v-gene usage
seurat_clusters_evolution_vgenetop_df = data.frame()
seurat_clusters_evolution_vgenetop_ratio_df = data.frame()
for( top_vgene_usage in 1:20){
  
  # Get the v-gene names of the top usage
  best_v.segment.exact.heavy = c( row.names( V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF)[1:top_vgene_usage])
  best_v.segment.exact.light = c( row.names( V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF)[1:top_vgene_usage])
  
  # Get the cells using the top v-genes
  cells_with_vgene_specific = phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$v.segment.exact.heavy %in% best_v.segment.exact.heavy
                                                                               & phenotype_clonotype_seurat_df$v.segment.exact.light %in% best_v.segment.exact.light)]
  # Dispersion of cells in the mice
  seurat_clusters.spe_df = as.data.frame.matrix( t( table( phenotype_clonotype_seurat_df[ which( phenotype_clonotype_seurat_df$plate.bcid %in% cells_with_vgene_specific), "seurat_clusters"])))
  
  # Compute the raio of cell in each Ag-specificity group
  nb_cell_seurat_clusters.spe = vector()
  for( current_seurat_clusters.spe in names( seurat_clusters.spe_df)){
    nb_cell_seurat_clusters.spe = append( nb_cell_seurat_clusters.spe, length( phenotype_clonotype_seurat_df$plate.bcid[ which( phenotype_clonotype_seurat_df$seurat_clusters == current_seurat_clusters.spe)]))
  }
  seurat_clusters.spe_ratio_df = signif( 100* seurat_clusters.spe_df / nb_cell_seurat_clusters.spe,3)
  
  
  # look at the Ag specificity of the cells ysing the top v-genes and accumulate it in a single dataframe
  seurat_clusters_evolution_vgenetop_df = rbind( seurat_clusters_evolution_vgenetop_df,   
                                                         seurat_clusters.spe_df)
  seurat_clusters_evolution_vgenetop_ratio_df = rbind( seurat_clusters_evolution_vgenetop_ratio_df,   
                                                               seurat_clusters.spe_ratio_df)
  
}
# Rename the row of the global dataframe 
row.names( seurat_clusters_evolution_vgenetop_df) = seq( 1, nrow( seurat_clusters_evolution_vgenetop_df))
row.names( seurat_clusters_evolution_vgenetop_ratio_df) = seq( 1, nrow( seurat_clusters_evolution_vgenetop_ratio_df))

# Normalize the row in percentages
seurat_clusters_evolution_vgenetop_df_pct = t(apply(seurat_clusters_evolution_vgenetop_df, 1, function(x)(signif( 100*x/sum(x), 3))))

# Print the count table and the percentage table
print( htmltools::tagList( datatable( seurat_clusters_evolution_vgenetop_df,
                                      caption = "Number of cells in fonction of the number of top V-gene usage\nin Cluster")))

print( htmltools::tagList( datatable( seurat_clusters_evolution_vgenetop_df_pct,
                                      caption = "Percentage of cells in fonction of the number of top V-gene usage\nin Cluster")))

print( htmltools::tagList( datatable( seurat_clusters_evolution_vgenetop_ratio_df,
                                      caption = "Percentage of cells in Cluster\nin fonction of the number of top V-gene usage")))


pheatmap( t( seurat_clusters_evolution_vgenetop_df), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Number of cells in fonction of the number of top V-gene usage\nin Cluster",
          display_numbers = TRUE, number_format = "%d")

pheatmap( t( seurat_clusters_evolution_vgenetop_df_pct), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Percentage of cells in fonction of the number of top V-gene usage\nin Cluster")

pheatmap( t( seurat_clusters_evolution_vgenetop_ratio_df), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Percentage of cells in Cluster\nin fonction of the number of top V-gene usage")

# Look at the chi2 statistic on the whole table
chisq_seurat_clusters_evolution_vgenetop_df = chisq.test(  seurat_clusters_evolution_vgenetop_df)
cat("<BR>")

pheatmap( t( chisq_seurat_clusters_evolution_vgenetop_df$residuals), cellwidth = 15, cellheight = 15, cluster_rows = FALSE, cluster_cols = FALSE,
          main = "Pearson residuals of chi2 test on \nNumber of cells in fonction of the number of top V-gene usage\nin Cluster")




