# ##############################################################################
# This script aims at analyzing the pehnotype and genotype in 
#  multi-dimension aspects
# ##############################################################################

## @knitr analyze_AGspecificity

# #########################################################################################################
# Group the cells by phenotype : non-Ag-specific (HA-NP-), HA-specific (HA+NP-), NP-specific (HA-NP+). 
# Further segregate non-Ag-specific into the three phenotypic subsets (1-CCR6+CXCR3-, 2-CCR6+CXCR3+ and 3-CCR6-CXCR3-).
# Compare BCR features in those 5 subsets : isotype, IGHV / IGKorL mutations, IGHV / IGKorL gene usage
# #########################################################################################################

# Define the Groups of cells
NON_AG_SPE_GROUPS = c( "nonAgSpe_CXCR3p.CCR6p", "nonAgSpe_CXCR3n.CCR6p", "nonAgSpe_CXCR3n.CCR6n")
SPE_GROUPS = c( "HASpe", "NPSpe")
GROUP_ORDER = c( NON_AG_SPE_GROUPS, SPE_GROUPS)
TOP_VGENE_USAGE = 10

# Store data to a new dataframe
multi_dimension_1_df = phenotype_clonotype_seurat_df
multi_dimension_1_df$Group = NA

# -- cells non-Ag-specific (HA-NP-) and CXCR3+CCR6+
multi_dimension_1_df[ which( multi_dimension_1_df$HAxNP == PHENOTYPE_HAxNP_NAMES[[ "NegNeg"]] &
                             multi_dimension_1_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "PosPos"]), "Group"] = "nonAgSpe_CXCR3p.CCR6p"
# -- cells non-Ag-specific (HA-NP-) and CXCR3+CCR6-
multi_dimension_1_df[ which( multi_dimension_1_df$HAxNP == PHENOTYPE_HAxNP_NAMES[[ "NegNeg"]] &
                               multi_dimension_1_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "NegPos"]), "Group"] = "nonAgSpe_CXCR3n.CCR6p"
# -- cells non-Ag-specific (HA-NP-) and CXCR3-CCR6-
multi_dimension_1_df[ which( multi_dimension_1_df$HAxNP == PHENOTYPE_HAxNP_NAMES[[ "NegNeg"]] &
                               multi_dimension_1_df$CXCR3xCCR6 == PHENOTYPE_CXCR3xCCR6_NAMES[ "NegNeg"]), "Group"] = "nonAgSpe_CXCR3n.CCR6n"
# -- cells HA-specific (HA+NP-)
multi_dimension_1_df[ which( multi_dimension_1_df$HAxNP == PHENOTYPE_HAxNP_NAMES[[ "PosNeg"]]), "Group"] = "HASpe"
# -- cells NP-specific (HA-NP+)
multi_dimension_1_df[ which( multi_dimension_1_df$HAxNP == PHENOTYPE_HAxNP_NAMES[[ "NegPos"]]), "Group"] = "NPSpe"

# Remove cells without a Group
multi_dimension_1_df = multi_dimension_1_df[ which( !is.na( multi_dimension_1_df$Group)),]

# Reorder groups suitably
multi_dimension_1_df$Group = factor( multi_dimension_1_df$Group, levels = GROUP_ORDER)

# Associate Groups and mouse of origin
multi_dimension_1_df$Mouse.Group = paste( multi_dimension_1_df$mouse, multi_dimension_1_df$Group, sep=".")

# ................................................................................................
## CELLS in GROUPS
# ................................................................................................

cat("<HR><H4>Dispersion of cells in the GROUPS</H4>")

# Plot the counts of cells in each group
ggplot( multi_dimension_1_df) + 
  geom_bar( aes( x=Group, fill=Group)) + theme_classic() +
  theme( axis.text.x = element_text( angle=45, hjust = 1)) +
  ggtitle( "Count of cells in each group")


# ................................................................................................
## CLONOTYPE SIZE CATEGORY versus GROUPS
# ................................................................................................

cat("<HR><H4>CLONOTYPE SIZE CATEGORY versus GROUPS</H4>")

# Show the dispersion of clonotype size category of heavy chain over Groups in a datatable
clonotype.size.category_group_table = as.data.frame.matrix( table( multi_dimension_1_df[, c( "Group", "clonotype.size.category")]))
clonotype.size.category_group_table_line_total = apply( clonotype.size.category_group_table, 1, sum)
clonotype.size.category_group_table_column_total = apply( clonotype.size.category_group_table, 2, sum)
clonotype.size.category_group_table = clonotype.size.category_group_table[ clonotype.size.category_group_table_line_total > 0, clonotype.size.category_group_table_column_total > 0]
print( htmltools::tagList( datatable( clonotype.size.category_group_table, caption = "Dispersion of clonotype size category among groups")))

# Show the dispersion of groups over isotype of heavy chain in percentage in a cumulative barplot
clonotype.size.category_group_wide = clonotype.size.category_group_table
clonotype.size.category_group_wide = 100*clonotype.size.category_group_wide / clonotype.size.category_group_table_line_total
clonotype.size.category_group_wide$Group = row.names( clonotype.size.category_group_wide)
clonotype.size.category_group_wide = reshape2::melt( clonotype.size.category_group_wide, id.vars = c( "Group"))
clonotype.size.category_group_wide$Group = factor( clonotype.size.category_group_wide$Group, levels = GROUP_ORDER)
print( ggplot( clonotype.size.category_group_wide) +
         geom_bar( aes( x = Group, fill = variable, y=value), stat="identity") + 
         labs( fill = "Clonotype size category", x= "Group", y = "Percentage") + 
         theme_minimal() + theme( axis.text.x = element_text( angle=45, hjust = 1)) +
         ggtitle( "Percentage of cell in group for clonotype size category")
)

# ................................................................................................
## CLONOTYPE SIZE CATEGORY versus GROUPS by MOUSE
# ................................................................................................

cat("<HR><H4>CLONOTYPE SIZE CATEGORY versus GROUPS by MOUSE</H4>")

# Show the dispersion of clonotype size category of heavy chain over Groups in a datatable
clonotype.size.category_mouse.group_table = as.data.frame.matrix( table( multi_dimension_1_df[, c( "Mouse.Group", "clonotype.size.category")]))
clonotype.size.category_mouse.group_table_line_total = apply( clonotype.size.category_mouse.group_table, 1, sum)
clonotype.size.category_mouse.group_table_column_total = apply( clonotype.size.category_mouse.group_table, 2, sum)
clonotype.size.category_mouse.group_table = clonotype.size.category_mouse.group_table[ clonotype.size.category_mouse.group_table_line_total > 0, clonotype.size.category_mouse.group_table_column_total > 0]
print( htmltools::tagList( datatable( clonotype.size.category_mouse.group_table, caption = "Dispersion of clonotype size category among groups by mouse")))

# Show the dispersion of groups over isotype of heavy chain in percentage in a cumulative barplot
clonotype.size.category_mouse.group_wide = clonotype.size.category_mouse.group_table
clonotype.size.category_mouse.group_wide = 100*clonotype.size.category_mouse.group_wide / clonotype.size.category_mouse.group_table_line_total
clonotype.size.category_mouse.group_wide$Mouse.Group = row.names( clonotype.size.category_mouse.group_wide)
clonotype.size.category_mouse.group_wide = reshape2::melt( clonotype.size.category_mouse.group_wide, id.vars = c( "Mouse.Group"))
print( ggplot( clonotype.size.category_mouse.group_wide) +
         geom_bar( aes( x = Mouse.Group, fill = variable, y=value), stat="identity") + 
         labs( fill = "Clonotype size category", x= "Mouse.Group", y = "Percentage") + 
         theme_minimal() + theme( axis.text.x = element_text( angle=45, hjust = 1)) +
         ggtitle( "Percentage of cell in group by mouse for clonotype size category")
)

# ................................................................................................
## HEAVY-CHAIN ISOTYPE versus GROUPS
# ................................................................................................

cat("<HR><H4>HEAVY-CHAIN ISOTYPE versus GROUPS</H4>")

# Show the dispersion of isotypes of heavy chain over Groups in a datatable
isotype_heavy_group_table = as.data.frame.matrix( table( multi_dimension_1_df[, c( "Group", "isotype.heavy")]))
isotype_heavy_group_table_line_total = apply( isotype_heavy_group_table, 1, sum)
isotype_heavy_group_table_column_total = apply( isotype_heavy_group_table, 2, sum)
isotype_heavy_group_table = isotype_heavy_group_table[ isotype_heavy_group_table_line_total > 0, isotype_heavy_group_table_column_total > 0]
print( htmltools::tagList( datatable( t(isotype_heavy_group_table), caption = "Dispersion of heavy chain isotypes over groups (all cells)")))

# Show the dispersion of groups over isotype of heavy chain in percentage in a cumulative barplot
isotype_heavy_group_wide = isotype_heavy_group_table
isotype_heavy_group_wide = 100*isotype_heavy_group_wide / isotype_heavy_group_table_line_total
isotype_heavy_group_wide$Group = row.names( isotype_heavy_group_wide)
isotype_heavy_group_wide = reshape2::melt( isotype_heavy_group_wide, id.vars = c( "Group"))
isotype_heavy_group_wide$Group = factor( isotype_heavy_group_wide$Group, levels = GROUP_ORDER)
print( ggplot( isotype_heavy_group_wide) +
         geom_bar( aes( x = Group, fill = variable, y=value), stat="identity") + 
         scale_fill_manual( values = isotype_heavy_palette) +
         labs( fill = "Isotype", x= "Group", y = "Percentage") + 
         theme_minimal() + theme( axis.text.x = element_text( angle=45, hjust = 1)) +
         ggtitle( "Percentage of cell in group for each heavy chain isotype")
)


# ................................................................................................
## LIGHT-CHAIN ISOTYPE versus GROUPS
# ................................................................................................

cat("<HR><H4>LIGHT-CHAIN ISOTYPE versus GROUPS</H4>")

# Show the dispersion of isotypes of light chain over Groups in a datatable
isotype_light_group_table = as.data.frame.matrix( table( multi_dimension_1_df[, c( "Group", "isotype.light")]))
isotype_light_group_table_line_total = apply( isotype_light_group_table, 1, sum)
isotype_light_group_table_column_total = apply( isotype_light_group_table, 2, sum)
isotype_light_group_table = isotype_light_group_table[ isotype_light_group_table_line_total > 0, isotype_light_group_table_column_total > 0]
print( htmltools::tagList( datatable( t(isotype_light_group_table), caption = "Dispersion of light chain isotypes over groups (all cells)")))

# Show the dispersion of groups over isotype of light chain in percentage in a cumulative barplot
isotype_light_group_wide = isotype_light_group_table
isotype_light_group_wide = 100*isotype_light_group_wide / isotype_light_group_table_line_total
isotype_light_group_wide$Group = row.names( isotype_light_group_wide)
isotype_light_group_wide = reshape2::melt( isotype_light_group_wide, id.vars = c( "Group"))
isotype_light_group_wide$Group = factor( isotype_light_group_wide$Group, levels = GROUP_ORDER)
print( ggplot( isotype_light_group_wide) +
         geom_bar( aes( x = Group, fill = variable, y=value), stat="identity") + 
         scale_fill_manual( values = isotype_light_palette) +
         labs( fill = "Isotype", x= "Group", y = "Percentage") + 
         theme_minimal() + theme( axis.text.x = element_text( angle=45, hjust = 1)) +
         ggtitle( "Percentage of cell in group for each light chain isotype")
)


# ................................................................................................
## HEAVY AND LIGHT-CHAIN NUMBER OF MUTATION versus GROUPS
# ................................................................................................


# Compose groups of mices to analyze : all mice together and then one mouse at a time
mice_groups = list( as.character( levels( multi_dimension_1_df$mouse)))
for( mice in levels( clonotype_seurat_df$mouse)){
  mice_groups = c( mice_groups, mice)
}

cat("<HR><H4>HEAVY-CHAIN NUMBER OF MUTATION versus GROUPS</H4>")

# Compute the percentage table at display the plots
for( current_mices in mice_groups){
  
  # Show the dispersion of number of mutation of heavy chain over groups in a datatable
  mutation.qual_heavy_group_df = multi_dimension_1_df[ which( multi_dimension_1_df$mouse %in% current_mices), c( "Group", "mutation.qual.heavy")]
  mutation.qual_heavy_group_df$mutation.qual.heavy = as.numeric( as.character( mutation.qual_heavy_group_df$mutation.qual.heavy))
  mutation.qual_heavy_group_table = as.data.frame.matrix( table( mutation.qual_heavy_group_df))
  mutation.qual_heavy_group_table_line_total = apply( mutation.qual_heavy_group_table, 1, sum)
  mutation.qual_heavy_group_table_column_total = apply( mutation.qual_heavy_group_table, 2, sum)
  mutation.qual_heavy_group_table = mutation.qual_heavy_group_table[ mutation.qual_heavy_group_table_line_total > 0, mutation.qual_heavy_group_table_column_total > 0]
  print( htmltools::tagList( datatable( t( mutation.qual_heavy_group_table), caption = "Dispersion of mutation number over groups (all cells)")))
  
  # Show the dispersion of groups CDR3 length of heavy chain in percentage in a cumulative barplot
  mutation.qual_heavy_group_wide = mutation.qual_heavy_group_table
  mutation.qual_heavy_group_wide = 100*mutation.qual_heavy_group_wide / mutation.qual_heavy_group_table_line_total
  mutation.qual_heavy_group_wide$Group = row.names( mutation.qual_heavy_group_wide)
  mutation.qual_heavy_group_wide = reshape2::melt( mutation.qual_heavy_group_wide, id.vars = c( "Group"))
  mutation.qual_heavy_group_wide$variable = as.numeric( as.character( mutation.qual_heavy_group_wide$variable))
  
  print( ggplot( mutation.qual_heavy_group_wide) +
           geom_bar( aes( x = variable, fill = Group, y=value), stat="identity", position="dodge") + 
           labs( fill = "Group", x= "# mutation", y = "Percentage") + 
           theme_minimal()  +
           theme(axis.text.x = element_text(angle = 45)) +
           ggtitle( paste( "Percentage of cell in group for each heavy chain number of mutation", 
                           "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_heavy_group_df),")"))
  )
  
  print( ggplot( mutation.qual_heavy_group_wide) +
           geom_bar( aes( x = variable, fill = Group, y=value), stat="identity", position="dodge") + 
           labs( fill = "Group", x= "# mutation", y = "Percentage") + 
           facet_wrap( . ~ Group) +
           theme_minimal()  +
           theme(axis.text.x = element_text(angle = 45)) +
           ggtitle( paste( "Percentage of cell in group for each heavy chain number of mutation", 
                           "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_heavy_group_df),")"))
  )
  
}


cat("<HR><H4>LIGHT-CHAIN NUMBER OF MUTATION versus GROUPS</H4>")

# Compute the percentage table at display the plots
for( current_mices in mice_groups){
  
  # Show the dispersion of number of mutation of light chain over groups in a datatable
  mutation.qual_light_group_df = multi_dimension_1_df[ which( multi_dimension_1_df$mouse %in% current_mices), c( "Group", "mutation.qual.light")]
  mutation.qual_light_group_df$mutation.qual.light = as.numeric( as.character( mutation.qual_light_group_df$mutation.qual.light))
  mutation.qual_light_group_table = as.data.frame.matrix( table( mutation.qual_light_group_df))
  mutation.qual_light_group_table_line_total = apply( mutation.qual_light_group_table, 1, sum)
  mutation.qual_light_group_table_column_total = apply( mutation.qual_light_group_table, 2, sum)
  mutation.qual_light_group_table = mutation.qual_light_group_table[ mutation.qual_light_group_table_line_total > 0, mutation.qual_light_group_table_column_total > 0]
  print( htmltools::tagList( datatable( t( mutation.qual_light_group_table), caption = "Dispersion of mutation number over groups (all cells)")))
  
  # Show the dispersion of groups CDR3 length of light chain in percentage in a cumulative barplot
  mutation.qual_light_group_wide = mutation.qual_light_group_table
  mutation.qual_light_group_wide = 100*mutation.qual_light_group_wide / mutation.qual_light_group_table_line_total
  mutation.qual_light_group_wide$Group = row.names( mutation.qual_light_group_wide)
  mutation.qual_light_group_wide = reshape2::melt( mutation.qual_light_group_wide, id.vars = c( "Group"))
  mutation.qual_light_group_wide$variable = as.numeric( as.character( mutation.qual_light_group_wide$variable))
  
  print( ggplot( mutation.qual_light_group_wide) +
           geom_bar( aes( x = variable, fill = Group, y=value), stat="identity", position="dodge") + 
           labs( fill = "Group", x= "# mutation", y = "Percentage") + 
           theme_minimal()  +
           theme(axis.text.x = element_text(angle = 45)) +
           ggtitle( paste( "Percentage of cell in group for each light chain number of mutation", 
                           "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_light_group_df),")"))
  )
  
  print( ggplot( mutation.qual_light_group_wide) +
           geom_bar( aes( x = variable, fill = Group, y=value), stat="identity", position="dodge") + 
           labs( fill = "Group", x= "# mutation", y = "Percentage") + 
           facet_wrap( . ~ Group) +
           theme_minimal()  +
           theme(axis.text.x = element_text(angle = 45)) +
           ggtitle( paste( "Percentage of cell in group for each light chain number of mutation", 
                           "\nin mice(s)", paste( current_mices, collapse = ","), "(number of cells=", nrow( mutation.qual_light_group_df),")"))
  )
  
}


# ................................................................................................
## HEAVY CHAIN V GENE USAGE versus GROUPS
# ................................................................................................

cat("<HR><H4>HEAVY-CHAIN V GENE USAGE versus GROUPS</H4>")

# Get the data table
v.segment.heavy_group_df = as.data.frame.matrix( t( table( multi_dimension_1_df[, c( "Group", "v.segment.heavy")])))

# Print the data in datatable
print( htmltools::tagList( datatable( v.segment.heavy_group_df)))

# Define the correlation function for the pairs plot
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method = "spearman"), digits=2)
  txt <- paste0("R = ", r)
  text(0.5, 0.5, txt,cex= 2)
}

# Plot the pairs scatterplot
pairs( v.segment.heavy_group_df, upper.panel = panel.cor)

# Compute the correlation and prepare the data for heatmap plot
v.segment.heavy_group_cor_df <- round( cor( v.segment.heavy_group_df, method = "spearman"), 2)
v.segment.heavy_group_cor_df[ upper.tri(v.segment.heavy_group_cor_df)] = NA
v.segment.heavy_group_cor_melt_df <- reshape2::melt( v.segment.heavy_group_cor_df)

# Plot the correlation heatmap
ggplot(data = v.segment.heavy_group_cor_melt_df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile( color = "white") +
  scale_fill_gradientn( colours = c( "blue", "white", "red"),
                        space = "Lab",
                        name = "Spearman correlation",
                        na.value = "white") +
  theme_classic() + theme( axis.text.x = element_text( angle=45, hjust = 1), axis.title = element_blank()) +
  ggtitle( "Heavy chain V gene usage correlation among groups")

# ................................................................................................
## LIGHT CHAIN V GENE USAGE versus GROUPS
# ................................................................................................

cat("<HR><H4>LIGHT-CHAIN V GENE USAGE versus GROUPS</H4>")

# Get the data table
v.segment.light_group_df = as.data.frame.matrix( t( table( multi_dimension_1_df[, c( "Group", "v.segment.light")])))

# Print the data in datatable
print( htmltools::tagList( datatable( v.segment.light_group_df)))

# Define the correlation function for the pairs plot
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method = "spearman"), digits=2)
  txt <- paste0("R = ", r)
  text(0.5, 0.5, txt,cex= 2)
}

# Plot the pairs scatterplot
pairs( v.segment.light_group_df, upper.panel = panel.cor)

# Compute the correlation and prepare the data for heatmap plot
v.segment.light_group_cor_df <- round( cor( v.segment.light_group_df, method = "spearman"), 2)
v.segment.light_group_cor_df[ upper.tri(v.segment.light_group_cor_df)] = NA
v.segment.light_group_cor_melt_df <- reshape2::melt( v.segment.light_group_cor_df)

# Plot the correlation heatmap
ggplot(data = v.segment.light_group_cor_melt_df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile( color = "white") +
  scale_fill_gradientn( colours = c( "blue", "white", "red"),
                        space = "Lab",
                        name = "Spearman correlation",
                        na.value = "white") +
  theme_classic() + theme( axis.text.x = element_text( angle=45, hjust = 1), axis.title = element_blank()) +
  ggtitle( "Light chain V gene usage correlation among groups")


# ................................................................................................
## HEAVY CHAIN EXACT V GENE USAGE versus GROUPS
# ................................................................................................

cat("<HR><H4>HEAVY-CHAIN EXACT V GENE USAGE versus GROUPS (Ag specific versus non-specific)</H4>")

# Look at the NPSpe group
# .........................

cat("<BR><H5>Heavy-Chain Exact V gene usage in NP Specific group</H5>")

# Get the data table in wide format
v.segment.exact.heavy_groupNPSPe_wide_df = as.data.frame.matrix( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == "NPSpe"), c( "Mouse.Group", "v.segment.exact.heavy")])))
# Sum the usage on mice
v.segment.exact.heavy_groupNPSPe_wide_df$total = apply( v.segment.exact.heavy_groupNPSPe_wide_df, 1, sum)
# Remove the V segment with a null total
v.segment.exact.heavy_groupNPSPe_wide_df = v.segment.exact.heavy_groupNPSPe_wide_df[ v.segment.exact.heavy_groupNPSPe_wide_df$total > 0, ]
# Order the usage in decreasing total
v.segment.exact.heavy_groupNPSPe_wide_df = v.segment.exact.heavy_groupNPSPe_wide_df[ order( v.segment.exact.heavy_groupNPSPe_wide_df$total, decreasing = TRUE),]
# Print the data in datatable
print( htmltools::tagList( datatable( v.segment.exact.heavy_groupNPSPe_wide_df)))

# Get the data table in long format 
v.segment.exact.heavy_groupNPSPe_long_df = as.data.frame( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == "NPSpe"), c( "Mouse.Group", "v.segment.exact.heavy")])))
# Order the V segment in decreasing usage
v.segment.exact.heavy_groupNPSPe_long_df$v.segment.exact.heavy = factor( v.segment.exact.heavy_groupNPSPe_long_df$v.segment.exact.heavy, levels = row.names( v.segment.exact.heavy_groupNPSPe_wide_df))
# Bar plot of the V gene usage per mouse in decreasing total number
ggplot( v.segment.exact.heavy_groupNPSPe_long_df) + 
  geom_bar( aes( x =v.segment.exact.heavy, y=Freq, fill=Mouse.Group), stat = "identity") +
  theme_classic() + theme( legend.position = "top") + theme( axis.text.x = element_text( angle = 45, hjust = 1, size = 4))

# Get the 5 best V gene usage
best_v.segment.exact.heavy_groupNPSPe = row.names( v.segment.exact.heavy_groupNPSPe_wide_df)[1:TOP_VGENE_USAGE]

# Look at the HASpe group
# .........................

cat("<BR><H5>Heavy-Chain Exact V gene usage in HA Specific group</H5>")

# Get the data table in wide format
v.segment.exact.heavy_groupHASpe_wide_df = as.data.frame.matrix( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == "HASpe"), c( "Mouse.Group", "v.segment.exact.heavy")])))
# Sum the usage on mice
v.segment.exact.heavy_groupHASpe_wide_df$total = apply( v.segment.exact.heavy_groupHASpe_wide_df, 1, sum)
# Remove the V segment with a null total
v.segment.exact.heavy_groupHASpe_wide_df = v.segment.exact.heavy_groupHASpe_wide_df[ v.segment.exact.heavy_groupHASpe_wide_df$total > 0, ]
# Order the usage in decreasing total
v.segment.exact.heavy_groupHASpe_wide_df = v.segment.exact.heavy_groupHASpe_wide_df[ order( v.segment.exact.heavy_groupHASpe_wide_df$total, decreasing = TRUE),]
# Print the data in datatable
print( htmltools::tagList( datatable( v.segment.exact.heavy_groupHASpe_wide_df)))

# Get the data table in long format 
v.segment.exact.heavy_groupHASpe_long_df = as.data.frame( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == "HASpe"), c( "Mouse.Group", "v.segment.exact.heavy")])))
# Order the V segment in decreasing usage
v.segment.exact.heavy_groupHASpe_long_df$v.segment.exact.heavy = factor( v.segment.exact.heavy_groupHASpe_long_df$v.segment.exact.heavy, levels = row.names( v.segment.exact.heavy_groupHASpe_wide_df))
# Bar plot of the V gene usage per mouse in decreasing total number
ggplot( v.segment.exact.heavy_groupHASpe_long_df) + 
  geom_bar( aes( x =v.segment.exact.heavy, y=Freq, fill=Mouse.Group), stat = "identity") +
  theme_classic() + theme( legend.position = "top") + theme( axis.text.x = element_text( angle = 45, hjust = 1, size = 4))

# Get the 5 best V gene usage
best_v.segment.exact.heavy_groupHASPe = row.names( v.segment.exact.heavy_groupHASpe_wide_df)[1:TOP_VGENE_USAGE]


# Compare the NPSpe and HASpe group best usage
# .............................................

cat("<BR><H5>Compare Heavy-Chain Exact V gene usage in HA and NP Specific groups</H5>")

# Accumulate the best V gene of each group in a single dataframe
best_v.segment.exact.heavy_groupsSPe_df = data.frame( v.segment.exact.heavy = unique( c( best_v.segment.exact.heavy_groupNPSPe, best_v.segment.exact.heavy_groupHASPe)),
                                                      NPSpe = v.segment.exact.heavy_groupNPSPe_wide_df[ unique( c( best_v.segment.exact.heavy_groupNPSPe, best_v.segment.exact.heavy_groupHASPe)), "total"],
                                                      HASpe = v.segment.exact.heavy_groupHASpe_wide_df[ unique( c( best_v.segment.exact.heavy_groupNPSPe, best_v.segment.exact.heavy_groupHASPe)), "total"])
best_v.segment.exact.heavy_groupsSPe_df[ is.na( best_v.segment.exact.heavy_groupsSPe_df)] = 0

# Scatterplot of the comparative usage of the best V genes
ggplot( best_v.segment.exact.heavy_groupsSPe_df, aes( x = NPSpe, y = HASpe)) + 
  geom_point( ) +
  geom_text_repel( aes( label = v.segment.exact.heavy)) + 
  theme_classic() + ggtitle( "Best V gene usage for HASpe and NPSpe groups (total over mice)")

# Look at the nonAgeSpe groups
# .............................

cat("<BR><H5>Heavy-Chain Exact V gene usage in Non specific sub-groups</H5>")

v.segment.exact.heavy_groupNonAgeSpe_wide_df_list = list()
all_best_v.segment.exact.heavy_groupNonAgeSpe = vector()

for( current_nonspe_group in NON_AG_SPE_GROUPS){
  # Get the data table in wide format
  v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]]  = as.data.frame.matrix( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == current_nonspe_group), c( "Mouse.Group", "v.segment.exact.heavy")])))
  # Sum the usage on mice
  v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] $total = apply( v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] , 1, sum)
  # Remove the V segment with a null total
  v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]]  = v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] [ v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] $total > 0, ]
  # Order the usage in decreasing total
  v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]]  = v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] [ order( v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] $total, decreasing = TRUE),]
  # Print the data in datatable
  print( htmltools::tagList( datatable( v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] )))
  
  # Get the data table in long format 
  v.segment.exact.heavy_groupNonAgeSpe_long_df = as.data.frame( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == current_nonspe_group), c( "Mouse.Group", "v.segment.exact.heavy")])))
  # Order the V segment in decreasing usage
  v.segment.exact.heavy_groupNonAgeSpe_long_df$v.segment.exact.heavy = factor( v.segment.exact.heavy_groupNonAgeSpe_long_df$v.segment.exact.heavy, levels = row.names( v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] ))
  # Bar plot of the V gene usage per mouse in decreasing total number
  print( 
    ggplot( v.segment.exact.heavy_groupNonAgeSpe_long_df) + 
      geom_bar( aes( x =v.segment.exact.heavy, y=Freq, fill=Mouse.Group), stat = "identity") +
      theme_classic() + theme( legend.position = "top") + theme( axis.text.x = element_text( angle = 45, hjust = 1, size = 4))
  )
  # Get the 5 best V gene usage
  all_best_v.segment.exact.heavy_groupNonAgeSpe = unique( append( all_best_v.segment.exact.heavy_groupNonAgeSpe, row.names( v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] )[1:TOP_VGENE_USAGE]))
}

# Compare Heavy-Chain Exact V gene usage in Non specific sub-groups
# .................................................................

cat("<BR><H5>Compare Heavy-Chain Exact V gene usage in Non specific sub-groups</H5>")

# Compose all the NonAg results in a single dataframe for their best genes
best_v.segment.exact.heavy_groupsNonAgeSpe_df = data.frame( v.segment.exact.heavy = all_best_v.segment.exact.heavy_groupNonAgeSpe,
                                                            nonAgSpe_CXCR3p.CCR6p = v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ "nonAgSpe_CXCR3p.CCR6p"]][ all_best_v.segment.exact.heavy_groupNonAgeSpe, "total"],
                                                            nonAgSpe_CXCR3n.CCR6p = v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ "nonAgSpe_CXCR3n.CCR6p"]][ all_best_v.segment.exact.heavy_groupNonAgeSpe, "total"],
                                                            nonAgSpe_CXCR3n.CCR6n = v.segment.exact.heavy_groupNonAgeSpe_wide_df_list[[ "nonAgSpe_CXCR3n.CCR6n"]][ all_best_v.segment.exact.heavy_groupNonAgeSpe, "total"])
best_v.segment.exact.heavy_groupsNonAgeSpe_df[ is.na( best_v.segment.exact.heavy_groupsNonAgeSpe_df)] = 0
best_v.segment.exact.heavy_groupsNonAgeSpe_df$total = apply( best_v.segment.exact.heavy_groupsNonAgeSpe_df[, NON_AG_SPE_GROUPS], 1, sum)

# Print the result in a datatable
print( htmltools::tagList( datatable( best_v.segment.exact.heavy_groupsNonAgeSpe_df, rownames = FALSE)))

# Look at the pairwise scatterplot of gene usage, using color to show the third group usage
for( current_nonspe_group_index1 in 1:(length( NON_AG_SPE_GROUPS)-1)){
  for( current_nonspe_group_index2 in (current_nonspe_group_index1+1):length( NON_AG_SPE_GROUPS)){
    # Get the group of the pairwise comparison
    current_nonspe_group_1 = NON_AG_SPE_GROUPS[ current_nonspe_group_index1]
    current_nonspe_group_2 = NON_AG_SPE_GROUPS[ current_nonspe_group_index2]
    # Get the third group for color scale
    missing_group = setdiff( NON_AG_SPE_GROUPS, c( current_nonspe_group_1, current_nonspe_group_2))
    
    # Scatterplot of the comparative usage of the best V genes
    print(
      ggplot( best_v.segment.exact.heavy_groupsNonAgeSpe_df, aes_string( x = current_nonspe_group_1, y = current_nonspe_group_2, color = missing_group)) + 
        geom_point( ) +
        geom_text_repel( aes( label = v.segment.exact.heavy)) + 
        theme_classic() + ggtitle( "Best V gene usage for Non Ag specific groups (total over mice)")
    )
  }
}

# Compare Heavy-Chain Exact V gene usage in Non specific groups and Specific groups
# .................................................................................

cat("<BR><H5>Compare Heavy-Chain Exact V gene usage in Non specific groups and Specific groups</H5>")

# Look at the best gene usage in common between AgSpe and nonAge Spe groups
all_best_v.segment.exact.heavy = intersect( best_v.segment.exact.heavy_groupsNonAgeSpe_df$v.segment.exact.heavy, best_v.segment.exact.heavy_groupsSPe_df$v.segment.exact.heavy)

# Compose a dataframe with the best gene usage of in common between AgSpe and nonAge Spe groups
best_v.segment.exact.heavy_allgroups_df = data.frame( v.segment.exact.heavy = all_best_v.segment.exact.heavy,
                                                      nonAgSpe = best_v.segment.exact.heavy_groupsNonAgeSpe_df[ which( best_v.segment.exact.heavy_groupsNonAgeSpe_df$v.segment.exact.heavy %in% all_best_v.segment.exact.heavy), "total"],
                                                      HASpe = best_v.segment.exact.heavy_groupsSPe_df[ which( best_v.segment.exact.heavy_groupsSPe_df$v.segment.exact.heavy %in% all_best_v.segment.exact.heavy), "HASpe"],
                                                      NPSpe = best_v.segment.exact.heavy_groupsSPe_df[ which( best_v.segment.exact.heavy_groupsSPe_df$v.segment.exact.heavy %in% all_best_v.segment.exact.heavy), "NPSpe"])

# Print the data in a datatable
print( htmltools::tagList( datatable(best_v.segment.exact.heavy_allgroups_df, rownames = FALSE, caption = "Best V-Genes in all groups")))


# Look at the V gene usage in all groups
# .......................................

cat("<BR><H5>Heavy-Chain Exact V gene usage in all groups</H5>")

# Get the data table in wide format
v.segment.exact.heavy_allgroups_wide_df = as.data.frame.matrix( t( table( multi_dimension_1_df[ , c( "Group", "v.segment.exact.heavy")])))
# Sum the usage on mice
v.segment.exact.heavy_allgroups_wide_df$total = apply( v.segment.exact.heavy_allgroups_wide_df, 1, sum)
# Remove the V segment with a null total
v.segment.exact.heavy_allgroups_wide_df = v.segment.exact.heavy_allgroups_wide_df[ v.segment.exact.heavy_allgroups_wide_df$total > 0, ]
# Order the usage in decreasing total
v.segment.exact.heavy_allgroups_wide_df = v.segment.exact.heavy_allgroups_wide_df[ order( v.segment.exact.heavy_allgroups_wide_df$total, decreasing = TRUE),]
# Print the data in datatable
print( htmltools::tagList( datatable( v.segment.exact.heavy_allgroups_wide_df)))

# Keep in constant variable for future usage
V.SEGMENT.EXACT.HEAVY_ALLGROUPS_DF = v.segment.exact.heavy_allgroups_wide_df

# Get the data table in long format 
v.segment.exact.heavy_allgroups_long_df = as.data.frame( t( table( multi_dimension_1_df[ c( "Group", "v.segment.exact.heavy")])))
# Order the V segment in decreasing usage
v.segment.exact.heavy_allgroups_long_df$v.segment.exact.heavy = factor( v.segment.exact.heavy_allgroups_long_df$v.segment.exact.heavy, levels = row.names( v.segment.exact.heavy_allgroups_wide_df))
# Bar plot of the V gene usage per mouse in decreasing total number
ggplot( v.segment.exact.heavy_allgroups_long_df) + 
  geom_bar( aes( x =v.segment.exact.heavy, y=Freq, fill=Group), stat = "identity") +
  theme_classic() + theme( legend.position = "top") + theme( axis.text.x = element_text( angle = 45, hjust = 1, size = 4))

# Compute the correlation and prepare the data for heatmap plot
v.segment.exact.heavy_allgroups_wide_df$total = NULL
v.segment.exact.heavy_allgroups_cor_df <- round( cor( v.segment.exact.heavy_allgroups_wide_df, method = "spearman"), 2)
v.segment.exact.heavy_allgroups_cor_df[ upper.tri(v.segment.exact.heavy_allgroups_cor_df)] = NA
v.segment.exact.heavy_allgroups_cor_melt_df <- reshape2::melt( v.segment.exact.heavy_allgroups_cor_df)

# Plot the correlation heatmap and save it to file
ggplot(data = v.segment.exact.heavy_allgroups_cor_melt_df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile( color = "white") +
  scale_fill_gradientn( colours = c( "blue", "white", "red"),
                        space = "Lab",
                        name = "Spearman correlation",
                        na.value = "white") +
  theme_classic() + theme( axis.text.x = element_text( angle=45, hjust = 1), axis.title = element_blank()) +
  ggtitle( "Heavy chain V gene usage correlation\namong all groups")

ggsave( file.path( PATH_ANALYSIS_OUTPUT, "GeneUsage_VHeavyChain_AllGroups.svg"))

# ................................................................................................
## LIGHT CHAIN EXACT V GENE USAGE versus GROUPS
# ................................................................................................

cat("<BR><BR><HR><H4>LIGHT-CHAIN EXACT V GENE USAGE versus GROUPS  (Ag specific versus non-specific)</H4>")

# Look at the NPSpe group
# .........................

cat("<BR><H5>Light-Chain Exact V gene usage in NP Specific group</H5>")

# Get the data table in wide format
v.segment.exact.light_groupNPSPe_wide_df = as.data.frame.matrix( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == "NPSpe"), c( "Mouse.Group", "v.segment.exact.light")])))
# Sum the usage on mice
v.segment.exact.light_groupNPSPe_wide_df$total = apply( v.segment.exact.light_groupNPSPe_wide_df, 1, sum)
# Remove the V segment with a null total
v.segment.exact.light_groupNPSPe_wide_df = v.segment.exact.light_groupNPSPe_wide_df[ v.segment.exact.light_groupNPSPe_wide_df$total > 0, ]
# Order the usage in decreasing total
v.segment.exact.light_groupNPSPe_wide_df = v.segment.exact.light_groupNPSPe_wide_df[ order( v.segment.exact.light_groupNPSPe_wide_df$total, decreasing = TRUE),]
# Print the data in datatable
print( htmltools::tagList( datatable( v.segment.exact.light_groupNPSPe_wide_df)))

# Get the data table in long format 
v.segment.exact.light_groupNPSPe_long_df = as.data.frame( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == "NPSpe"), c( "Mouse.Group", "v.segment.exact.light")])))
# Order the V segment in decreasing usage
v.segment.exact.light_groupNPSPe_long_df$v.segment.exact.light = factor( v.segment.exact.light_groupNPSPe_long_df$v.segment.exact.light, levels = row.names( v.segment.exact.light_groupNPSPe_wide_df))
# Bar plot of the V gene usage per mouse in decreasing total number
ggplot( v.segment.exact.light_groupNPSPe_long_df) + 
  geom_bar( aes( x =v.segment.exact.light, y=Freq, fill=Mouse.Group), stat = "identity") +
  theme_classic() + theme( legend.position = "top") + theme( axis.text.x = element_text( angle = 45, hjust = 1, size = 4))

# Get the 5 best V gene usage
best_v.segment.exact.light_groupNPSPe = row.names( v.segment.exact.light_groupNPSPe_wide_df)[1:TOP_VGENE_USAGE]

# Look at the HASpe group
# .........................

cat("<BR><H5>Light-Chain Exact V gene usage in HA Specific group</H5>")

# Get the data table in wide format
v.segment.exact.light_groupHASpe_wide_df = as.data.frame.matrix( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == "HASpe"), c( "Mouse.Group", "v.segment.exact.light")])))
# Sum the usage on mice
v.segment.exact.light_groupHASpe_wide_df$total = apply( v.segment.exact.light_groupHASpe_wide_df, 1, sum)
# Remove the V segment with a null total
v.segment.exact.light_groupHASpe_wide_df = v.segment.exact.light_groupHASpe_wide_df[ v.segment.exact.light_groupHASpe_wide_df$total > 0, ]
# Order the usage in decreasing total
v.segment.exact.light_groupHASpe_wide_df = v.segment.exact.light_groupHASpe_wide_df[ order( v.segment.exact.light_groupHASpe_wide_df$total, decreasing = TRUE),]
# Print the data in datatable
print( htmltools::tagList( datatable( v.segment.exact.light_groupHASpe_wide_df)))

# Get the data table in long format 
v.segment.exact.light_groupHASpe_long_df = as.data.frame( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == "HASpe"), c( "Mouse.Group", "v.segment.exact.light")])))
# Order the V segment in decreasing usage
v.segment.exact.light_groupHASpe_long_df$v.segment.exact.light = factor( v.segment.exact.light_groupHASpe_long_df$v.segment.exact.light, levels = row.names( v.segment.exact.light_groupHASpe_wide_df))
# Bar plot of the V gene usage per mouse in decreasing total number
ggplot( v.segment.exact.light_groupHASpe_long_df) + 
  geom_bar( aes( x =v.segment.exact.light, y=Freq, fill=Mouse.Group), stat = "identity") +
  theme_classic() + theme( legend.position = "top") + theme( axis.text.x = element_text( angle = 45, hjust = 1, size = 4))

# Get the 5 best V gene usage
best_v.segment.exact.light_groupHASPe = row.names( v.segment.exact.light_groupHASpe_wide_df)[1:TOP_VGENE_USAGE]


# Compare the NPSpe and HASpe group best usage
# .............................................

cat("<BR><H5>Compare Light-Chain Exact V gene usage in HA and NP Specific groups</H5>")

# Accumulate the best V gene of each group in a single dataframe
best_v.segment.exact.light_groupsSPe_df = data.frame( v.segment.exact.light = unique( c( best_v.segment.exact.light_groupNPSPe, best_v.segment.exact.light_groupHASPe)),
                                                      NPSpe = v.segment.exact.light_groupNPSPe_wide_df[ unique( c( best_v.segment.exact.light_groupNPSPe, best_v.segment.exact.light_groupHASPe)), "total"],
                                                      HASpe = v.segment.exact.light_groupHASpe_wide_df[ unique( c( best_v.segment.exact.light_groupNPSPe, best_v.segment.exact.light_groupHASPe)), "total"])
best_v.segment.exact.light_groupsSPe_df[ is.na( best_v.segment.exact.light_groupsSPe_df)] = 0

# Scatterplot of the comparative usage of the best V genes
ggplot( best_v.segment.exact.light_groupsSPe_df, aes( x = NPSpe, y = HASpe)) + 
  geom_point( ) +
  geom_text_repel( aes( label = v.segment.exact.light)) + 
  theme_classic() + ggtitle( "Best V gene usage for HASpe and NPSpe groups (total over mice)")

# Look at the nonAgeSpe groups
# .............................

cat("<BR><H5>Light-Chain Exact V gene usage in Non specific sub-groups</H5>")

v.segment.exact.light_groupNonAgeSpe_wide_df_list = list()
all_best_v.segment.exact.light_groupNonAgeSpe = vector()

for( current_nonspe_group in NON_AG_SPE_GROUPS){
  # Get the data table in wide format
  v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]]  = as.data.frame.matrix( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == current_nonspe_group), c( "Mouse.Group", "v.segment.exact.light")])))
  # Sum the usage on mice
  v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] $total = apply( v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] , 1, sum)
  # Remove the V segment with a null total
  v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]]  = v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] [ v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] $total > 0, ]
  # Order the usage in decreasing total
  v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]]  = v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] [ order( v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] $total, decreasing = TRUE),]
  # Print the data in datatable
  print( htmltools::tagList( datatable( v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] )))
  
  # Get the data table in long format 
  v.segment.exact.light_groupNonAgeSpe_long_df = as.data.frame( t( table( multi_dimension_1_df[ which( multi_dimension_1_df$Group == current_nonspe_group), c( "Mouse.Group", "v.segment.exact.light")])))
  # Order the V segment in decreasing usage
  v.segment.exact.light_groupNonAgeSpe_long_df$v.segment.exact.light = factor( v.segment.exact.light_groupNonAgeSpe_long_df$v.segment.exact.light, levels = row.names( v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] ))
  # Bar plot of the V gene usage per mouse in decreasing total number
  print( 
    ggplot( v.segment.exact.light_groupNonAgeSpe_long_df) + 
      geom_bar( aes( x =v.segment.exact.light, y=Freq, fill=Mouse.Group), stat = "identity") +
      theme_classic() + theme( legend.position = "top") + theme( axis.text.x = element_text( angle = 45, hjust = 1, size = 4))
  )
  # Get the 5 best V gene usage
  all_best_v.segment.exact.light_groupNonAgeSpe = unique( append( all_best_v.segment.exact.light_groupNonAgeSpe, row.names( v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ current_nonspe_group]] )[1:TOP_VGENE_USAGE]))
}

# Compare Light-Chain Exact V gene usage in Non specific sub-groups
# ..................................................................

cat("<BR><H5>Compare Light-Chain Exact V gene usage in Non specific sub-groups</H5>")

# Compose all the NonAg results in a single dataframe for their best genes
best_v.segment.exact.light_groupsNonAgeSpe_df = data.frame( v.segment.exact.light = all_best_v.segment.exact.light_groupNonAgeSpe,
                                                            nonAgSpe_CXCR3p.CCR6p = v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ "nonAgSpe_CXCR3p.CCR6p"]][ all_best_v.segment.exact.light_groupNonAgeSpe, "total"],
                                                            nonAgSpe_CXCR3n.CCR6p = v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ "nonAgSpe_CXCR3n.CCR6p"]][ all_best_v.segment.exact.light_groupNonAgeSpe, "total"],
                                                            nonAgSpe_CXCR3n.CCR6n = v.segment.exact.light_groupNonAgeSpe_wide_df_list[[ "nonAgSpe_CXCR3n.CCR6n"]][ all_best_v.segment.exact.light_groupNonAgeSpe, "total"])
best_v.segment.exact.light_groupsNonAgeSpe_df[ is.na( best_v.segment.exact.light_groupsNonAgeSpe_df)] = 0
best_v.segment.exact.light_groupsNonAgeSpe_df$total = apply( best_v.segment.exact.light_groupsNonAgeSpe_df[, NON_AG_SPE_GROUPS], 1, sum)

# Print the result in a datatable
print( htmltools::tagList( datatable( best_v.segment.exact.light_groupsNonAgeSpe_df, rownames = FALSE)))

# Look at the pairwise scatterplot of gene usage, using color to show the third group usage
for( current_nonspe_group_index1 in 1:(length( NON_AG_SPE_GROUPS)-1)){
  for( current_nonspe_group_index2 in (current_nonspe_group_index1+1):length( NON_AG_SPE_GROUPS)){
    # Get the group of the pairwise comparison
    current_nonspe_group_1 = NON_AG_SPE_GROUPS[ current_nonspe_group_index1]
    current_nonspe_group_2 = NON_AG_SPE_GROUPS[ current_nonspe_group_index2]
    # Get the third group for color scale
    missing_group = setdiff( NON_AG_SPE_GROUPS, c( current_nonspe_group_1, current_nonspe_group_2))
        
    # Scatterplot of the comparative usage of the best V genes
    print(
      ggplot( best_v.segment.exact.light_groupsNonAgeSpe_df, aes_string( x = current_nonspe_group_1, y = current_nonspe_group_2, color = missing_group)) + 
        geom_point( ) +
        geom_text_repel( aes( label = v.segment.exact.light)) + 
        theme_classic() + ggtitle( "Best V gene usage for Non Ag specific groups (total over mice)")
    )
  }
}

# Compare Light-Chain Exact V gene usage in Non specific groups and Specific groups
# ..................................................................................

cat("<BR><H5>Compare Light-Chain Exact V gene usage in Non specific groups and Specific groups</H5>")

# Look at the best gene usage in common between AgSpe and nonAge Spe groups
all_best_v.segment.exact.light = intersect( best_v.segment.exact.light_groupsNonAgeSpe_df$v.segment.exact.light, best_v.segment.exact.light_groupsSPe_df$v.segment.exact.light)

# Compose a dataframe with the best gene usage of in common between AgSpe and nonAge Spe groups
best_v.segment.exact.light_allgroups_df = data.frame( v.segment.exact.light = all_best_v.segment.exact.light,
                                                      nonAgSpe = best_v.segment.exact.light_groupsNonAgeSpe_df[ which( best_v.segment.exact.light_groupsNonAgeSpe_df$v.segment.exact.light %in% all_best_v.segment.exact.light), "total"],
                                                      HASpe = best_v.segment.exact.light_groupsSPe_df[ which( best_v.segment.exact.light_groupsSPe_df$v.segment.exact.light %in% all_best_v.segment.exact.light), "HASpe"],
                                                      NPSpe = best_v.segment.exact.light_groupsSPe_df[ which( best_v.segment.exact.light_groupsSPe_df$v.segment.exact.light %in% all_best_v.segment.exact.light), "NPSpe"])

# Print the data in a datatable
print( htmltools::tagList( datatable(best_v.segment.exact.light_allgroups_df, rownames = FALSE, caption = "Best V-Genes in all groups")))

 
# Look at the V gene usage in all groups
# .......................................

cat("<BR><H5>Light-Chain Exact V gene usage in all groups</H5>")

# Get the data table in wide format
v.segment.exact.light_allgroups_wide_df = as.data.frame.matrix( t( table( multi_dimension_1_df[ , c( "Group", "v.segment.exact.light")])))
# Sum the usage on mice
v.segment.exact.light_allgroups_wide_df$total = apply( v.segment.exact.light_allgroups_wide_df, 1, sum)
# Remove the V segment with a null total
v.segment.exact.light_allgroups_wide_df = v.segment.exact.light_allgroups_wide_df[ v.segment.exact.light_allgroups_wide_df$total > 0, ]
# Order the usage in decreasing total
v.segment.exact.light_allgroups_wide_df = v.segment.exact.light_allgroups_wide_df[ order( v.segment.exact.light_allgroups_wide_df$total, decreasing = TRUE),]
# Print the data in datatable
print( htmltools::tagList( datatable( v.segment.exact.light_allgroups_wide_df)))

# Keep in constant variable for future usage
V.SEGMENT.EXACT.LIGHT_ALLGROUPS_DF = v.segment.exact.light_allgroups_wide_df

# Get the data table in long format 
v.segment.exact.light_allgroups_long_df = as.data.frame( t( table( multi_dimension_1_df[ c( "Group", "v.segment.exact.light")])))
# Order the V segment in decreasing usage
v.segment.exact.light_allgroups_long_df$v.segment.exact.light = factor( v.segment.exact.light_allgroups_long_df$v.segment.exact.light, levels = row.names( v.segment.exact.light_allgroups_wide_df))
# Bar plot of the V gene usage per mouse in decreasing total number
ggplot( v.segment.exact.light_allgroups_long_df) + 
  geom_bar( aes( x =v.segment.exact.light, y=Freq, fill=Group), stat = "identity") +
  theme_classic() + theme( legend.position = "top") + theme( axis.text.x = element_text( angle = 45, hjust = 1, size = 4))

# Compute the correlation and prepare the data for heatmap plot
v.segment.exact.light_allgroups_wide_df$total = NULL
v.segment.exact.light_allgroups_cor_df <- round( cor( v.segment.exact.light_allgroups_wide_df, method = "spearman"), 2)
v.segment.exact.light_allgroups_cor_df[ upper.tri(v.segment.exact.light_allgroups_cor_df)] = NA
v.segment.exact.light_allgroups_cor_melt_df <- reshape2::melt( v.segment.exact.light_allgroups_cor_df)

# Plot the correlation heatmap and save it to file
ggplot(data = v.segment.exact.light_allgroups_cor_melt_df, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile( color = "white") +
  scale_fill_gradientn( colours = c( "blue", "white", "red"),
                        space = "Lab",
                        name = "Spearman correlation",
                        na.value = "white") +
  theme_classic() + theme( axis.text.x = element_text( angle=45, hjust = 1), axis.title = element_blank()) +
  ggtitle( "Light chain V gene usage correlation\namong all groups")

ggsave( file.path( PATH_ANALYSIS_OUTPUT, "GeneUsage_VLightChain_AllGroups.svg"))




# ................................................................................................
## CLONAL OVERLAP BETWEEN AG GROUPS
# ................................................................................................

cat("<BR><H4>Clonotype overlap between Ag groups</H4>")

# Get the list of clonotype without the mouse names
clonotype.no.mouse = gsub( "moLung\\d.", "", multi_dimension_1_df$full.clonotype.similarity.symbol, fixed = FALSE)

# Compute the intersection of the clonotype among groups
clonotype.no.mouse_allgroups_table = as.data.frame.matrix( table( clonotype.no.mouse, multi_dimension_1_df$Group))

# List the clonotype for each group
area = list()
for( group in names( clonotype.no.mouse_allgroups_table)){
  area[[ group]] = row.names( clonotype.no.mouse_allgroups_table)[ which( clonotype.no.mouse_allgroups_table[ , group] > 0)]
}

# Compute the table of pairwise intersection between groups
intersection_allgroups_df = data.frame()
for( group1_index in 1:length( names( clonotype.no.mouse_allgroups_table))){
  for( group2_index in 1 : length( names( clonotype.no.mouse_allgroups_table))){
    group1 = names( clonotype.no.mouse_allgroups_table)[ group1_index]
    group2 = names( clonotype.no.mouse_allgroups_table)[ group2_index]
    if( group1_index > group2_index){
      intersection = length( intersect( area[[ group1]], area[[ group2]]))
    }else{
      intersection = NA
    }
    intersection_allgroups_df = rbind( intersection_allgroups_df, data.frame( group1 = group1,
                                                                              group2 = group2,
                                                                              intersect = intersection))
  }
}

# Plot the table of pairwise intersection between groups and save it to file
ggplot(data = intersection_allgroups_df, aes(x=group1, y=group2, fill=intersect)) + 
  geom_tile( color = "white") +
  scale_fill_gradientn( colours = c( "blue", "white", "red"),
                        space = "Lab",
                        name = "Count",
                        na.value = "white") +
  theme_classic() + theme( axis.text.x = element_text( angle=45, hjust = 1), axis.title = element_blank()) +
  ggtitle( "Shared clontype counts among all groups")

ggsave( file.path( PATH_ANALYSIS_OUTPUT, "SharedClonotypeAcrossAgGroups_HeatMap_Counts.svg"))


# Draw a 5-sets Venn Diagram of shared clonotype between the Ag groups
venn.plot.all = VennDiagram::draw.quintuple.venn(area1 = length( area[[ 1]]), 
           area2 = length( area[[ 2]]), 
           area3 = length( area[[ 3]]), 
           area4 = length( area[[ 4]]), 
           area5 = length( area[[ 5]]), 
           n12 = length( intersect( area[[ 1]], area[[ 2]])), 
           n13 = length( intersect( area[[ 1]], area[[ 3]])), 
           n14 = length( intersect( area[[ 1]], area[[ 4]])), 
           n15 = length( intersect( area[[ 1]], area[[ 5]])),
           n23 = length( intersect( area[[ 2]], area[[ 3]])), 
           n24 = length( intersect( area[[ 2]], area[[ 4]])), 
           n25 = length( intersect( area[[ 2]], area[[ 5]])), 
           n34 = length( intersect( area[[ 3]], area[[ 4]])), 
           n35 = length( intersect( area[[ 3]], area[[ 5]])), 
           n45 = length( intersect( area[[ 4]], area[[ 5]])), 
           n123 = length( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 3]])), 
           n124 = length( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 4]])), 
           n125 = length( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 5]])), 
           n134 = length( intersect( intersect(area[[ 1]], area[[ 3]]), area[[ 4]])),
           n135 = length( intersect( intersect(area[[ 1]], area[[ 3]]), area[[ 5]])), 
           n145 = length( intersect( intersect(area[[ 1]], area[[ 4]]), area[[ 5]])), 
           n234 = length( intersect( intersect(area[[ 2]], area[[ 3]]), area[[ 4]])), 
           n235 = length( intersect( intersect(area[[ 2]], area[[ 3]]), area[[ 5]])), 
           n245 = length( intersect( intersect(area[[ 2]], area[[ 4]]), area[[ 5]])), 
           n345 = length( intersect( intersect(area[[ 3]], area[[ 4]]), area[[ 5]])), 
           n1234 = length( intersect( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 3]]), area[[ 4]])), 
           n1235 = length( intersect( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 3]]), area[[ 5]])),
           n1245 = length( intersect( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 4]]), area[[ 5]])), 
           n1345 = length( intersect( intersect( intersect(area[[ 1]], area[[ 3]]), area[[ 4]]), area[[ 5]])), 
           n2345 = length( intersect( intersect( intersect(area[[ 2]], area[[ 3]]), area[[ 4]]), area[[ 5]])), 
           n12345 = length( intersect( intersect( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 3]]), area[[ 4]]), area[[ 5]])),
           category = names( area),
           fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
           cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
           margin = 0.3, ind = FALSE
)
svg( file.path( PATH_ANALYSIS_OUTPUT, "SharedClonotypeAcrossAgGroups_AllGroups_VennDiagram.svg"), width = 7, height = 7)
  grid.draw( venn.plot.all)
dev.off()


venn.plot.HASpe = VennDiagram::draw.quad.venn( area1 = length( area[[ 1]]), 
                             area2 = length( area[[ 2]]), 
                             area3 = length( area[[ 3]]), 
                             area4 = length( area[[ 4]]), 
                             n12 = length( intersect( area[[ 1]], area[[ 2]])), 
                             n13 = length( intersect( area[[ 1]], area[[ 3]])), 
                             n14 = length( intersect( area[[ 1]], area[[ 4]])), 
                             n23 = length( intersect( area[[ 2]], area[[ 3]])), 
                             n24 = length( intersect( area[[ 2]], area[[ 4]])), 
                             n34 = length( intersect( area[[ 3]], area[[ 4]])), 
                             n123 = length( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 3]])), 
                             n124 = length( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 4]])), 
                             n134 = length( intersect( intersect(area[[ 1]], area[[ 3]]), area[[ 4]])),
                             n234 = length( intersect( intersect(area[[ 2]], area[[ 3]]), area[[ 4]])), 
                             n1234 = length( intersect( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 3]]), area[[ 4]])),
                             category = names( area)[ 1:4],
                             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                             cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
                             margin = 0.3, ind = FALSE
                             )

svg( file.path( PATH_ANALYSIS_OUTPUT, "SharedClonotypeAcrossAgGroups_HASpevsGroups_VennDiagram.svg"), width = 7, height = 7)
grid.draw( venn.plot.HASpe)
dev.off()

venn.plot.NPSpe = VennDiagram::draw.quad.venn( area1 = length( area[[ 1]]), 
                             area2 = length( area[[ 2]]), 
                             area3 = length( area[[ 3]]), 
                             area4 = length( area[[ 5]]), 
                             n12 = length( intersect( area[[ 1]], area[[ 2]])), 
                             n13 = length( intersect( area[[ 1]], area[[ 3]])), 
                             n14 = length( intersect( area[[ 1]], area[[ 5]])), 
                             n23 = length( intersect( area[[ 2]], area[[ 3]])), 
                             n24 = length( intersect( area[[ 2]], area[[ 5]])), 
                             n34 = length( intersect( area[[ 3]], area[[ 5]])), 
                             n123 = length( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 3]])), 
                             n124 = length( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 5]])), 
                             n134 = length( intersect( intersect(area[[ 1]], area[[ 3]]), area[[ 5]])),
                             n234 = length( intersect( intersect(area[[ 2]], area[[ 3]]), area[[ 5]])), 
                             n1234 = length( intersect( intersect( intersect(area[[ 1]], area[[ 2]]), area[[ 3]]), area[[ 5]])),
                             category = names( area)[ c(1:3,5)],
                             fill = c("dodgerblue", "goldenrod1", "darkorange1", "orchid3"),
                             cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "orchid3"),
                             margin = 0.3, ind = FALSE
)

svg( file.path( PATH_ANALYSIS_OUTPUT, "SharedClonotypeAcrossAgGroups_NPSpevsGroups_VennDiagram.svg"), width = 7, height = 7)
grid.draw( venn.plot.NPSpe)
dev.off()

# ................................................................................................
## EXPORT MULTIDIMENSIONNAL DATA TO FILE INCLUDING HEAVY AND LIGHT CHAIN SEQUENCES
# ................................................................................................

# Declare the file where the export is done
export_file = file.path( PATH_ANALYSIS_OUTPUT, "AGspecificity_with_sequences.tsv")
cat("<BR><H4>Exporting AG specificity analysis</H4>")
cat("<bExporting AG specificity analysis and all data including chain sequences to<BR>", export_file, "</b>")

# Add the chain sequences to dataframe
export_multi_dimention =multi_dimension_1_df
export_multi_dimention$sequence.heavy = heavy_and_light_chains_df[ export_multi_dimention$plate.bcid, "sequence.heavy"]
export_multi_dimention$sequence.light = heavy_and_light_chains_df[ export_multi_dimention$plate.bcid, "sequence.light"]

# Write the file
write.table( export_multi_dimention,
             file = file.path( PATH_ANALYSIS_OUTPUT, "AGspecificity_with_sequences.tsv"),
             sep ="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


