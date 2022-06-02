###############################################################################
# This file defines SAMPLE parameters as global variables that will be loaded
# before analysis starts. It should only define values for variables that were
# not set in the PROJECT parameters file, or values that must differ from the
# globally defined parameters (commented by default here).
#




#### General

SAMPLE_NAME  = "10935372";

SUBSAMPLE_BY_TISSUE = list( LG = c( "M1LG", "M2LG", "M3LG"),
                            LN = c( "M1LN", "M2LN", "M3LN"),
                            SP = c( "M1SP", "M2SP", "M3SP")
                      )
 
MOUSE_TO_SAMPLE_DF = data.frame( M1 = c( "M1LG", "M1LN", "M1SP"),
                                 M2 = c( "M2LG", "M2LN", "M2SP"),
                                 M3 = c( "M3LG", "M3LN", "M3SP"))
row.names( MOUSE_TO_SAMPLE_DF) = c( "LG", "LN", "SP")

#### Input / Output


#### Add manually other parameters that need to be customized for this sample 
#### Use 'globalParams.R' as template

