# ############################################################################################
# This script aims to compare the gene expression accross tissues
# ############################################################################################

## @knitr marker_genes_between_tissue

# Set the Sample identity as the active Idents
Idents( sc10x.rna.seurat) = "MULTI_Sample"

# Get the list of tissue names from the list of sample name
tissue_only = sapply( Idents( sc10x.rna.seurat), function( sample){
  
  # Get the current sample name as string
  sample_str = as.character( sample)
  
  # Get the dataframe position of the current sample name in the dataframe containing both tissue and mouse information
  indexes = which( MOUSE_TO_SAMPLE_DF == sample_str, arr.ind = TRUE)
  
  # If the sample name is found in the dataframe, extract the name of the tissue corresponding to the sample as the corresponding row name of the dataframe
  # If the sample name is not found in the dataframe, return NA
  if( length(indexes) == 2){
    return( row.names(MOUSE_TO_SAMPLE_DF)[  indexes[ 1]])
  }else{
    return( NA)
  }
})

# Add the Tissue information as metadata in the Seurat object
sc10x.rna.seurat = AddMetaData( sc10x.rna.seurat, metadata = tissue_only, col.name = "Sample_Tissue")

# Find markers between tissues
Idents( sc10x.rna.seurat) = "Sample_Tissue"
tissues_markers_df = FindAllMarkers( sc10x.rna.seurat, assay = "RNA")
ordered_tissues_markers_df = tissues_markers_df[ , c( "gene", "cluster", "avg_logFC", "pct.1", "pct.2", "p_val", "p_val_adj")]

# Show the datatable of marker genes
datatable( ordered_tissues_markers_df,
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Avg logFC", "% tissue1", "% other tissues", "p-value", "Adjusted p-value")
           # extensions = c('Buttons', 'Select'),
           # options = list( dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
           #                autoWidth = FALSE,
           #                buttons = exportButtonsListDT,
           #                columnDefs = list( 
           #                  list( # Center all columns except first one
           #                    targets = 1:(ncol( variableAnnotationsDT)-1),
           #                    className = 'dt-center'),
           #                  list( # Set renderer function for 'float' type columns
           #                    targets = 1:(ncol( variableAnnotationsDT)-1),
           #                    render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}"))), 
           #                #fixedColumns = TRUE, # Does not work well with filter on this column
           #                #fixedHeader = TRUE, # Does not work well with 'scrollX'
           #                lengthMenu = list(c( 10, 50, 100, -1),
           #                                  c( 10, 50, 100, "All")),
           #                orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
           #                processing = TRUE, 
           #                #search.regex= TRUE, # Does not work well with 'search.smart'
           #                search.smart = TRUE,
           #                select = TRUE, # Enable ability to select rows
           #                scrollCollapse = TRUE,
           #                scroller = TRUE,  # Only load visible data
           #                scrollX = TRUE,
           #                scrollY = "525px",
           #                stateSave = TRUE)
           )
