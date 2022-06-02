###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

LITERAL_TITLE = "Analysis of BCR repertoire"
ANALYSIS_STEP_NAME = "08_BCRRepertoire"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

# Path to the output of migmap to get the details on sequences of BCR
PATH_BCR_MIGMAP_LIGHTCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/4_migmap/IGKLV_filtered_for_all.csv")
PATH_BCR_MIGMAP_HEAVYCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/4_migmap/IGHV_filtered_for_all.csv")

PATH_BCR_BLAST_LIGHTCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/4_blast/10938436BCR_S5_blastn_csr_LightCstRegion.out")
PATH_BCR_BLAST_HEAVYCHAIN_ANALYSIS_FILE = file.path( PATH_EXPERIMENT_OUTPUT, "03_BCRAnalysis/4_blast/10938436BCR_S5_blastn_csr_HeavyCstRegion.out")

# The length of the FR4 kept for the similiaty sequence analysis
FR4_AA_LENGTH_CUT = 36

# Path to the Seurat object to use
PATH_SEURAT_FILTERED_OBJECT = file.path( PATH_EXPERIMENT_OUTPUT, "06_GlobalHeterogeneityNoBCR")

PATH_SEURAT_FILTERED_OBJECT_BY_TISSUE = list( 
                                              LG = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue", "sc10x.rna.seurat_LG.RDS"),
                                              LN = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue", "sc10x.rna.seurat_LN.RDS"),
                                              SP = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue", "sc10x.rna.seurat_SP.RDS")
                                          )
# Boolean indicating if the ClustalOmega must be done
# If FALSE, the analysis will try to get file from previous ClustalOmega run
# If TRUE, the ClustalOmega analysis will be launched and the possible previous resutls will be removed
FORCE_CLUSTALO = FALSE
