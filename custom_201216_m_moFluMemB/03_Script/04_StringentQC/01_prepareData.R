# #########################################
# This script reads and filters sc10x  data
# #########################################


# -----------------------------------------------
# -----------------------------------------------
# READ DATA
# -----------------------------------------------
# -----------------------------------------------

## @knitr loadData

cat("<BR>Loading data fom", PATH_RNA_COUNT_TABLE)

# Load the data from file
rna_raw_count_df = read.table( PATH_RNA_COUNT_TABLE, header = TRUE, sep=",", quote = '"')

# Get a unique cellID keeping only the cell plate ID and barcode ID
rna_raw_count_df$cell.plate.bcid = unlist( sapply( rna_raw_count_df$UniqueCellID, function( code){
  
  # Convert factor to char
  code = as.character( code)
  
  # get the first "_plate" index since the plate+barcodeID is after
  index_underscore = regexpr( "_plate", code)
  if( index_underscore > 0){
    code = substr( code , start = index_underscore+1, stop = nchar( code))
  }
  
  return( code)
}), use.names = FALSE)

# Put the plateID+barcodeID as row name
row.names( rna_raw_count_df) = rna_raw_count_df$cell.plate.bcid

# Create the Seurat object from the expresion matrix
custom.rna.seurat = CreateSeuratObject( counts = t( rna_raw_count_df[ , !names( rna_raw_count_df) %in% c( "cell.plate.bcid", "UniqueCellID")]), 
                                        min.cells = LOAD_MIN_CELLS, 
                                        min.features = LOAD_MIN_FEATURES, 
                                        project = SAMPLE_NAME);

# Attribute a numeric ID to each cell (easier than barcode)
custom.rna.seurat[[ "numID"]] = 1:length( Cells( custom.rna.seurat));




cat("<BR><BR>Loading meta-data fom", PATH_CELL_PREPROCESSING_METADATA)

preprocessing_metadata_df = read.table( file = PATH_CELL_PREPROCESSING_METADATA, 
                                        header = TRUE, 
                                        sep=",", 
                                        stringsAsFactors = FALSE, 
                                        quote = '"',
                                        row.names = 1)

# Find the plate and barcodeID of each cell to create a simplified unique identifier of the cell
# Get it from the UniqueCellID which is composed of 201216_m_moFluMemB_plate1_BC80
preprocessing_metadata_df$cell.plate.bcid = unlist( sapply( preprocessing_metadata_df$UniqueCellID, function( code){
  
  # get the first "_plate" index since the plate+barcodeID is after
  index_underscore = regexpr( "_plate", code)
  if( index_underscore > 0){
    code = substr( code , start = index_underscore+1, stop = nchar( code))
  }
  
  return( code)
}), use.names = FALSE)

# Put cell.plateid;bcid as row.names
row.names( preprocessing_metadata_df) = preprocessing_metadata_df$cell.plate.bcid

# Compare the cells in the meta-data with the cells in the Seurat object
if( length( setdiff( Cells( custom.rna.seurat), preprocessing_metadata_df$cell.plate.bcid)) > 0 || 
            length( setdiff( preprocessing_metadata_df$cell.plate.bcid, Cells( custom.rna.seurat)))  > 0 ){
  cat("<BR><B>ERROR:</B>The list of cells are different between expression data and meta-data")
  stop( "ERROR:The list of cells are different between expression data and meta-data")
}

# Order the dataframe like in the Seurat object
preprocessing_metadata_df = preprocessing_metadata_df[ Cells( custom.rna.seurat), ]

# Add Mouse of origin as meta-data in Seurat object
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, metadata = preprocessing_metadata_df$Mouse_ID, col.name = "mouse")

# Add cell plate as meta-data in Seurat object
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, metadata = preprocessing_metadata_df$Plate, col.name = "plate")

# Add tissue as meta-data in Seurat object
custom.rna.seurat = AddMetaData( object = custom.rna.seurat, metadata = rep( "LG", nrow( preprocessing_metadata_df)), col.name = "tissue")



# -----------------------------------------------
# -----------------------------------------------
# FILTER DATA
# -----------------------------------------------
# -----------------------------------------------

## @knitr filterData

# FILTER THE CELLS USING THE GENE EXPRESSION MATRIX
# ====================================================

cat("<H4>filter the cells using gene expression</H4>")

### Identify mitocondrial genes in matrix
mito.genes = grep( pattern = "^mt\\.", x = rownames(GetAssayData(object = custom.rna.seurat, slot = "counts")), value = TRUE, ignore.case = TRUE)
if(length(mito.genes)==0){
  warning( "No mitochondrial genes could be identified in this dataset.");
} else{
  # Compute percentage of mitochondrial transcripts in each cell and 
  # Add the mitocondrial gene percentage as meta information in the Seurat object
  # Note : Specify features instead of pattern so we can use ignore.case in grep call
  custom.rna.seurat[["percent.mito"]] <- PercentageFeatureSet(custom.rna.seurat, features=mito.genes)
}

## Identify Ribisomal genes
# Compute percentage of mitochondrial transcripts in each cell and 
# Add the ribosomal gene percentage as meta information in the Seurat object
# Note : Specify features instead of pattern so we can use ignore.case in grep call
ribo.genes = grep( pattern = "^Rp[sl][[:digit:]]", x = rownames(GetAssayData(object = custom.rna.seurat, slot = "counts")), value = TRUE, ignore.case = FALSE)
if(length(ribo.genes)==0){
  warning( "No ribosomal genes could be identified in this dataset.");
} else{
  custom.rna.seurat[["percent.ribo"]] <- PercentageFeatureSet(custom.rna.seurat, features=ribo.genes)
}

### Identify cells that will be rejected based on specified thresholds
# ....................................................................

# Identify cells to reject based on UMI numbers
nUMI.drop = logical( length( Cells(custom.rna.seurat)));
if( ! is.null( FILTER_UMI_MIN)){
  nUMI.drop = nUMI.drop | (custom.rna.seurat[["nCount_RNA", drop=TRUE]] < FILTER_UMI_MIN);
}
if( ! is.null( FILTER_UMI_MAX)){
  nUMI.drop = nUMI.drop | (custom.rna.seurat[["nCount_RNA", drop=TRUE]] > FILTER_UMI_MAX);
}

# Identify cells to reject based on number of expressed genes
nGene.drop = logical( length(Cells(custom.rna.seurat)));
if( ! is.null( FILTER_FEATURE_MIN)){
  nGene.drop = nGene.drop | (custom.rna.seurat[["nFeature_RNA", drop=TRUE]] < FILTER_FEATURE_MIN);
}
if( ! is.null( FILTER_FEATURE_MAX)){
  nGene.drop = nGene.drop | (custom.rna.seurat[["nFeature_RNA", drop=TRUE]] > FILTER_FEATURE_MAX);
}

# Identify cells with high percentage of mitocondrial genes
mito.drop = logical( length(Cells(custom.rna.seurat)));
if( length(mito.genes) && (! is.null(FILTER_MITOPCT_MAX))){
  mito.drop = (custom.rna.seurat[["percent.mito", drop=TRUE]] > FILTER_MITOPCT_MAX);
}

# Identify cells with low percentage of ribosomal genes
ribo.drop = logical( length(Cells(custom.rna.seurat)));
if( length( ribo.genes) && (! is.null( FILTER_RIBOPCT_MIN))){
  ribo.drop = (custom.rna.seurat[["percent.ribo", drop=TRUE]] < FILTER_RIBOPCT_MIN);
}

### Plot distributions of #UMIs, #Genes and %Mito among cells
# ....................................................................

# Gather data to be visualized together (cell name + numID + metrics)
cellsData = cbind( "Cell" = colnames( custom.rna.seurat), # cell names from rownames conflict with search field in datatable
                   custom.rna.seurat[[c( "numID", "nCount_RNA", "nFeature_RNA")]],
                   "percent.mito" = if(length( mito.genes)) as.numeric(format(custom.rna.seurat[["percent.mito", drop = TRUE]], digits = 5)) else NULL,
                   "percent.ribo" = if(length( ribo.genes)) as.numeric(format(custom.rna.seurat[["percent.ribo", drop = TRUE]], digits = 5)) else NULL);

# Create text to show under cursor for each cell
hoverText = do.call(paste, c(Map( paste, c( "",
                                            "Cell ID: ", "# UMIs: ", "# Genes: ",
                                            if(length( mito.genes)) "% Mito: ",
                                            if(length( ribo.genes)) "% Ribo: "), cellsData, sep=""), sep="\n"))

# Generate plotly violin/jitter panels for #umis, #genes, and %mitochondrial stats
panelWidth = 180;
panelHeight = 400;
lypanel_umis  = plotViolinJitter(cbind( cellsData, outliers = nUMI.drop), xAxisFormula = ~as.numeric(1), yAxisFormula = ~nCount_RNA, pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, xTicklabels = "# UMIs", thresholdHigh = FILTER_UMI_MAX, thresholdLow = FILTER_UMI_MIN, panelWidth = panelWidth, panelHeight = panelHeight);
lypanel_genes = plotViolinJitter(cbind( cellsData, outliers = nGene.drop), xAxisFormula = ~as.numeric(1), yAxisFormula = ~nFeature_RNA, pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, xTicklabels = "# Genes", thresholdHigh = FILTER_FEATURE_MAX, thresholdLow = FILTER_FEATURE_MIN, panelWidth = panelWidth, panelHeight = panelHeight);
lypanel_mitos = if(length(mito.genes)) plotViolinJitter(cbind( cellsData, outliers = mito.drop), xAxisFormula = ~as.numeric(1), yAxisFormula = ~percent.mito, pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, xTicklabels = "% Mito", thresholdHigh = FILTER_MITOPCT_MAX, panelWidth = panelWidth, panelHeight = panelHeight) else NULL;
lypanel_ribos = if(length(ribo.genes)) plotViolinJitter(cbind( cellsData, outliers = ribo.drop), xAxisFormula = ~as.numeric(1), yAxisFormula = ~percent.ribo, pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, xTicklabels = "% Ribo", thresholdHigh = FILTER_RIBOPCT_MIN, panelWidth = panelWidth, panelHeight = panelHeight) else NULL;


# Set panels as a list and define plotly config
distribPanelsList = list( lypanel_umis, lypanel_genes, lypanel_mitos, lypanel_ribos);
distribPanelsList = lapply( distribPanelsList, config, displaylogo = FALSE, 
                            toImageButtonOptions = list(format='svg'), 
                            modeBarButtons = list(list('toImage'), list('zoom2d', 'pan2d', 'resetScale2d')));

# Control layout using flex because subplot is limited regarding plot sizing and alignment
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
     div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
          lapply( distribPanelsList, div, style = paste("flex : 0 0 auto; margin: 5px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
#)

### Plot scatterplots of pairs of #UMIs, #Genes and %Mito among cells
# ....................................................................

lypanel_feature_vs_numi= plotScatterPlotly(cbind( cellsData, outliers = nUMI.drop | nGene.drop), xAxisFormula = ~nCount_RNA, yAxisFormula = ~nFeature_RNA, xAxisTitle = "# UMI", yAxisTitle = "# Features", pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, panelWidth = 380, panelHeight = 380);
lypanel_ribo_vs_mito = if(length(ribo.genes) && length( mito.genes)) plotScatterPlotly(cbind( cellsData, outliers = ribo.drop | mito.drop), xAxisFormula = ~percent.mito, yAxisFormula = ~percent.ribo, xAxisTitle = "% Mito genes", yAxisTitle = "% Ribo genes", pointsColorFormula = ~ifelse(outliers, "#FF0000", "#444444"), hoverText = hoverText, panelWidth = 380, panelHeight = 380) else NULL;

# Set panels as a list and define plotly config
scatterPanelsList = list( lypanel_ribo_vs_mito, lypanel_feature_vs_numi );
scatterPanelsList = lapply( scatterPanelsList, config, displaylogo = FALSE, 
                            toImageButtonOptions = list(format='svg'), 
                            modeBarButtons = list(list('toImage'), list('zoom2d', 'pan2d', 'resetScale2d')));

# Control layout using flex because subplot is limited regarding plot sizing and alignment
# 'browsable' required in console, not in script/document
#browsable(
div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
     div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
          lapply( scatterPanelsList, div, style = paste("flex : 0 0 auto; margin: 5px;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
#)

### Display information on the removed cells
# ....................................................................

cat( "<br>Number of cells removed based on number of UMIs:", sum( nUMI.drop));
cat( "<br>Number of cells removed based on number of genes:", sum( nGene.drop));
if(exists( "mito.drop")) cat( "<br>Number of cells removed based on high percentage of mitochondrial transcripts:", sum( mito.drop));
if(exists( "ribo.drop")) cat( "<br>Number of cells removed based on low percentage of ribosomal transcripts:", sum( ribo.drop));

expression_cell.plate.bcid_toexclude = (nUMI.drop | nGene.drop | (if(exists( "mito.drop")) mito.drop else FALSE) | (if(exists( "ribo.drop")) ribo.drop else FALSE));

cat("<BR>Number of cells to exclude :", sum( expression_cell.plate.bcid_toexclude))
cat("<BR>Number of cells to keep:", length( expression_cell.plate.bcid_toexclude) - sum( expression_cell.plate.bcid_toexclude))


# FILTER THE CELLS USING THE META-DATA FROM PRE-PROCESSING ANALYSIS
# =================================================================

cat("<H4>filter the cells using metadata from preprocessing analysis</H4>")

# Look at the QC on the meta-data

## Look at the number of Feature_RNA
qc_plot1 = ggplot( preprocessing_metadata_df) + 
  geom_violin( aes( x = Donor, y=nFeature_RNA, fill=Donor)) +
  geom_jitter( aes( x = Donor, y=nFeature_RNA, fill=Donor), width = 0.3, size = 1) +
  geom_hline( yintercept = FILTER_FEATURE_MIN, color = "red") +
  theme_classic() + theme( legend.position = "None") +
  ggtitle( "Distribution of nFeature_RNA \nand minimal QC limit")

cell.plate.bcid_toofewFeature = preprocessing_metadata_df$nFeature_RNA < FILTER_FEATURE_MIN
cell.plate.bcid_toofewFeature[ is.na( cell.plate.bcid_toofewFeature)] = TRUE
names( cell.plate.bcid_toofewFeature) = preprocessing_metadata_df$cell.plate.bcid

## Look at the percentage of mitochondrial genes
qc_plot2 = ggplot( preprocessing_metadata_df) + 
  geom_violin( aes( x = Donor, y=percent.mito, fill=Donor)) +
  geom_jitter( aes( x = Donor, y=percent.mito, fill=Donor), width = 0.3, size = 1) +
  geom_hline( yintercept = FILTER_MITOPCT_MAX/100, color = "red") +
  theme_classic() + theme( legend.position = "None") +
  ggtitle( "Distribution of mitochondrial percentage \nand maximal QC limit")

cell.plate.bcid_toomanyMito = preprocessing_metadata_df$percent.mito > (FILTER_MITOPCT_MAX/100)
cell.plate.bcid_toomanyMito[ is.na( cell.plate.bcid_toomanyMito)] = TRUE
names( cell.plate.bcid_toomanyMito) = preprocessing_metadata_df$cell.plate.bcid

## Look at the percentage of ERCC
qc_plot3 = ggplot( preprocessing_metadata_df) + 
  geom_violin( aes( x = Donor, y=percent.ERCC, fill=Donor)) +
  geom_jitter( aes( x = Donor, y=percent.ERCC, fill=Donor), width = 0.3, size = 1) +
  geom_hline( yintercept = FILTER_ERCC_MAX/100, color = "red") +
  theme_classic() + theme( legend.position = "None") +
  ggtitle( "Distribution of ERCC percentage \nand maximal QC limit")

cell.plate.bcid_toomanyERCC = preprocessing_metadata_df$percent.ERCC > (FILTER_ERCC_MAX/100)
cell.plate.bcid_toomanyERCC[ is.na( cell.plate.bcid_toomanyERCC)] = TRUE
names( cell.plate.bcid_toomanyERCC) = preprocessing_metadata_df$cell.plate.bcid

## Look at the Accuracy pearcon correlation value
qc_plot4 = ggplot( preprocessing_metadata_df) + 
  geom_violin( aes( x = Donor, y=AccuracyPearsonERCC, fill=Donor)) +
  geom_jitter( aes( x = Donor, y=AccuracyPearsonERCC, fill=Donor), width = 0.3, size = 1) +
  geom_hline( yintercept = FILTER_ACCURACYPEARSONERCC_MIN, color = "red") +
  theme_classic() + theme( legend.position = "None") +
  ggtitle( "Distribution of Accuracy Pearson correlation value \nand minimal QC limit")

cell.plate.bcid_toolowAccuracyPearsonERCC = (preprocessing_metadata_df$AccuracyPearsonERCC < FILTER_ACCURACYPEARSONERCC_MIN)
cell.plate.bcid_toolowAccuracyPearsonERCC[ is.na( cell.plate.bcid_toolowAccuracyPearsonERCC)] = TRUE
names( cell.plate.bcid_toolowAccuracyPearsonERCC) = preprocessing_metadata_df$cell.plate.bcid

#Plot the graph of the QC analysis
qc_plot1 + qc_plot2 + qc_plot3 + qc_plot4

# Show the number of cell filtered by the QC
cat("<BR>Number of cells with too few Feature_RNA :", sum( cell.plate.bcid_toofewFeature))
cat("<BR>Number of cells with too high mitochondrial gene percentage :", sum( cell.plate.bcid_toomanyMito))
cat("<BR>Number of cells with too high ERCC percentage :", sum( cell.plate.bcid_toomanyERCC))
cat("<BR>Number of cells with too low Accuracy Pearson correlation value :", sum( cell.plate.bcid_toolowAccuracyPearsonERCC))

metadata_cell.plate.bcid_toexclude = (cell.plate.bcid_toofewFeature | cell.plate.bcid_toomanyMito | cell.plate.bcid_toomanyERCC | cell.plate.bcid_toolowAccuracyPearsonERCC)

cat("<BR>Number of cells to exclude :", sum( metadata_cell.plate.bcid_toexclude))
cat("<BR>Number of cells to keep:", length( metadata_cell.plate.bcid_toexclude) - sum( metadata_cell.plate.bcid_toexclude))



# MERGE THE TWO QC ANALYSIS
# =================================================================

cat( "<H4>Merge the QC methods to filter out the cells with al the filters</H4>")

# Compare the list of cells to exclude between the two methods
# ....................................................................

diff_nfeature = length( nGene.drop) - sum(nGene.drop == cell.plate.bcid_toofewFeature)
if( diff_nfeature > 0){
  cat("<BR><BR><B>WARNING:</B>List of cells excluded due to a low number of RNA_feature are different. Number of different cells:", diff_nfeature)
}

diff_pctmito = length( mito.drop) - sum(mito.drop == cell.plate.bcid_toomanyMito)
if( diff_pctmito > 0){
  cat("<BR><BR><B>WARNING:</B>List of cells excluded due to a high percentage of mitochondrial genes are different. Number of different cells:", diff_pctmito)
}

### Identify the cells to exclude as union of cells with low nb UMI, low nb of expressed genes, 
### high percentage of mitocondrial genes or low percentage of ribosomal genes
# ....................................................................

custom.rna.seurat[["outlier"]] = (  nUMI.drop | 
                                    nGene.drop | 
                                    (if(exists( "mito.drop")) mito.drop else FALSE) | 
                                    (if(exists( "ribo.drop")) ribo.drop else FALSE) |
                                    cell.plate.bcid_toomanyERCC |
                                    cell.plate.bcid_toolowAccuracyPearsonERCC);

cat("<br><br>Removed cells after filters:", sum( unlist(custom.rna.seurat[["outlier"]] )));
cat("<br>Remaining cells after filters:", sum( ! unlist(custom.rna.seurat[["outlier"]] )));

### Exclude the outliers cells from the Seurat RNA object
# ....................................................................
cat("<BR><br>Number of cells in Seurat object before filtering:" , length( Cells( custom.rna.seurat)))
custom.rna.seurat = subset( custom.rna.seurat, subset = (outlier == FALSE))
cat("<BR><br>Number of cells in Seurat object after filtering:" , length( Cells( custom.rna.seurat)))

### Record which cells got rejected
# ....................................................................

# Export the excluded cells to file
write.table( data.frame( cellid = names( which( nUMI.drop))), file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_nbUMI.txt")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table( data.frame( cellid = names( which( nGene.drop))), file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_nbGene.txt")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table( data.frame( cellid = names( which( mito.drop))), file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctMito.txt")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table( data.frame( cellid = names( which( ribo.drop))), file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctRibo.txt")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table( data.frame( cellid = names( which( cell.plate.bcid_toomanyERCC))), file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_pctERCC.txt")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
write.table( data.frame( cellid = names( which( cell.plate.bcid_toolowAccuracyPearsonERCC))), file= file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "excluded_cells_AccuracyPearsonERCC.txt")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")


# NORMALIZE DATA
# --------------

## @knitr normalizeData

# Normalize RNA data
custom.rna.seurat = NormalizeData( object = custom.rna.seurat,
                                  normalization.method = DATA_NORM_METHOD,
                                  scale.factor = DATA_NORM_SCALEFACTOR,
                                  verbose = .VERBOSE);

custom.rna.seurat = ScaleData( object    = custom.rna.seurat,
                              do.center = DATA_CENTER,
                              do.scale  = DATA_SCALE,
                              vars.to.regress = DATA_VARS_REGRESS,
                              verbose = .VERBOSE);


