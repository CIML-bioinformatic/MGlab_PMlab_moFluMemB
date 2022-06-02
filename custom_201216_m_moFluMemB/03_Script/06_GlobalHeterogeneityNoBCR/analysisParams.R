###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

LITERAL_TITLE = "Heterogeneity analysis through all samples with removal of BCR genes"
ANALYSIS_STEP_NAME = "06_GlobalHeterogeneityNoBCR"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

PATH_SEURAT_FILTERED_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT,  "04_QC")

# Maximum number of variable features to keep for PCA analysis
VARIABLE_FEATURES_MAXNB   = 2000;

# Maximum number of variable features to keep for table presentation
VARIABLE_FEATURES_SHOWTOP = 200;

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION = 0.5;

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 20;  # Number of dimensions to use from PCA results

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 10;       # Number of marker genes to show in report and tables (NULL for all)


#### List of genes to be monitored

MONITORED_GENES = list( 
                        SELECTION      = c( "Ccr6", "Cnn3", "Lmo2", "Fcer2a", "Serpinb1a", "Cxcr3", "Fscn1", "Plac8", "Cst3", "Itga4", "Ahr",
                                            "Rgs1", "Lmnb1", "Lta", "Ffar2", "S100a6", "Cd55", "Klf2", "Vim", "Tuba1a", "Marcksl1", "Apoe", 
                                            "Cd24a", "Fcmr", "Ms4a6c", "Lck", "Nfatc1", "Zbtb32", "Cr2", "Il9r", "Dtx1", "S1pr1", "S1pr3", 
                                            "Sell", "Cd22", "Dusp2", "Slpi", "Adgre1", "Tbx21", "Itgax"),
                        GC             = c( "Aicda", "Bach2", "Aurka", "Bcl2a1b", "Bcl6", "Ccnb1", "Ccnb2", "Ccne1", "Cenpa", "Efnb1", "Ezh2", "Fas", "Gna13", "Il4i1", "Mcm6", "Mki67", "Polh", "Rad51", "S1pr2", "Slamf1"),
                        GC_DZ          = c( "Cd86", "Cxcr4"),
                        GC_LZ          = c( "Cd83", "Cxcr5"),
                        NAIVE_B        = c( "Ccr6", "Ccr7", "Cd2", "Ighd", "Cd38", "S1pr1", "Sell"),
                        NAIVE_GC_B     = c( "Ptprc", "Cd19", "Ms4a1"),
                        PLASMA         = c( "Prdm1", "Xbp1", "Igha", "Ighm", "Jchain", "Bsg", "Irf4", "Sdc1", "Xbp1", "Spn"),
                        FACS           = c( "eYFP"));


# Genes expressed by putative contaminating cells
#CONTAMINATION_GENES = c();
#MONITORED_GENES = c(MONITORED_GENES, list(CONTAMINATION  = CONTAMINATION_GENES));



