###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used byt the current
# analysis
#

ANALYSIS_STEP_NAME="09_CompareTissues"

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME);
LITERAL_TITLE = "Compare cells in different tissues"



# The list of all markers genes per tissue and per cluster in tissue
PATH_ALL_MARKER_GENES_TABLE_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/allMarkers_LG.tab"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/allMarkers_LN.tab"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/allMarkers_SP.tab")
)

# The list of top markers genes per tissue and per cluster in tissue
PATH_TOP_MARKER_GENES_TABLE_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/topMarkersDT_LG.tab"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/topMarkersDT_LN.tab"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/topMarkersDT_SP.tab")
)

# The list of UMAP embedding per tissue
PATH_UMAP_EMBEDDING_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/UMAP_Embedding_LG.tab"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/UMAP_Embedding_LN.tab"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/UMAP_Embedding_SP.tab")
)

# The list of cluster mapping per tissue
PATH_CLUSTER_MAPPING_FILE = list(
  LG = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/clusters_LG.tab"),
  LN = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/clusters_LN.tab"),
  SP = file.path( PATH_EXPERIMENT_OUTPUT, "07_GlobalHeterogeneityNoBCRByTissue/clusters_SP.tab")
)
