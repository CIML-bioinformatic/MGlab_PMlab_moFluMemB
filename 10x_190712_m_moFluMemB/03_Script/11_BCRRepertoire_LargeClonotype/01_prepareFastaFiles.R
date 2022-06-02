# ###########################################################################################################
# This script aims to obtain the fasta files for IGH and IGKL root cells and for all clonootypes cells
# ###########################################################################################################

library( "seqinr")

# Load the data on clonotype cells and their IG sequences
clonotype_cell_data_df = read.csv(  PATH_LARGE_CLONOTYPE_CELL_FILE, stringsAsFactors = FALSE, header = TRUE, sep=",",check.names=FALSE)

# Concatenate IGH and IGKL
clonotype_cell_data_df$IGH_IGKL_seq = paste0(clonotype_cell_data_df$IGHseq,clonotype_cell_data_df$IGKLseq)

# Define a color for the Seurat clusters
clonotype_cell_data_df$Color = sapply( clonotype_cell_data_df$seurat_clusters, function(x) if(x ==0){"red"}else{"blue"})

# Select the cell IGH information
IGH_df = clonotype_cell_data_df[ , c( "UniqueCellID", "full.clonotype.similarity.symbol", "IGHseq", "Root")]

# Select the root cell IGH information
IGH_df_root = IGH_df[ IGH_df$Root == "Y", ]

# Select the cell IGKL information
IGKL_df = clonotype_cell_data_df[,c("UniqueCellID","full.clonotype.similarity.symbol","IGKLseq","Root")]

# Select the root cell IGKL information
IGKL_df_root = IGKL_df[ IGKL_df$Root == "Y", ]

# Split the cell data by clonotypes
clonotype_cell_data_df_byClonotype_list <- split( clonotype_cell_data_df , f = clonotype_cell_data_df$full.clonotype.similarity.symbol )

# Define the output folder for this script
PATH_SCRIPT_OUTPUT = file.path( PATH_ANALYSIS_OUTPUT, "01_FastaFiles")
if(!dir.exists(PATH_SCRIPT_OUTPUT)){
  dir.create( PATH_SCRIPT_OUTPUT,showWarnings=TRUE)
}

# Write IGH sequences fasta files of root cells
write.fasta( as.list( toupper( IGH_df_root$IGHseq)), 
             IGH_df_root$UniqueCellID, 
             file.path( PATH_SCRIPT_OUTPUT,"Root_IGH.fasta"),
             open = "w", nbchar = 60, as.string = TRUE)

# Write IGKL sequences fasta files of root cells
write.fasta( as.list( toupper( IGKL_df_root$IGKLseq)), 
             IGKL_df_root$UniqueCellID, 
             file.path( PATH_SCRIPT_OUTPUT,"Root_IGKL.fasta"),
             open = "w", nbchar = 60, as.string = TRUE)

# Write fasta files of concatenation of IGH/IGKL sequences by clonotypes
for( clonotype_name in names(clonotype_cell_data_df_byClonotype_list)){
  PATH_SCRIPT_OUTPUT_BY_CLONOTYPE = file.path( PATH_SCRIPT_OUTPUT ,clonotype_name)
  dir.create( PATH_SCRIPT_OUTPUT_BY_CLONOTYPE, showWarnings=FALSE)
    
  write.fasta (as.list( toupper( clonotype_cell_data_df_byClonotype_list[[ clonotype_name]]$IGH_IGKL_seq)),
               clonotype_cell_data_df_byClonotype_list[[ clonotype_name]]$UniqueCellID,
               file.path( PATH_SCRIPT_OUTPUT_BY_CLONOTYPE, "merged_IGH_IGKL.fasta"),
               open = "w", nbchar = 60, as.string = TRUE)
    
  cell_color_df = clonotype_cell_data_df_byClonotype_list[[ clonotype_name]][ , c("UniqueCellID","Color")]
  # add uca in 1 first row
  cell_color_df <- rbind(c("uca","black"), df_tmp)
    
  write.table( cell_color_df ,file.path( PATH_SCRIPT_OUTPUT_BY_CLONOTYPE, "color.csv"),
               row.names=FALSE, col.names=FALSE, sep =',', quote = FALSE)
}

#
#print(IGKL_df_root$UniqueCellID)
#print(IGKL_df_root$UniqueCellID)




