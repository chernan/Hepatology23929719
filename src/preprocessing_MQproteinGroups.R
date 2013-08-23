
#############################################################################################
#
# SliceSILAC pipeline : Step 0
#
# Clean and filter data from a MaxQuant proteinGroups.txt 
# 
# Includes: 
# - removing Contaminants (MaxQuant label: Contaminant)
# - removing Reverse (MaxQuant label: Reverse)
# - removing Only Identified by Site = TRUE (MaxQuant label: Only.identified.by.site)
# - removing Ratio Count < 3 (MaxQuant label: Ratio.H.L.Count)
#
#
# Usage:
#  Rscript preprocessing_MQproteinGroups.R <input_file_name> <output_file_name>
# 
# 
# Script version : 1.0
# Date: 22-08-2013
# Author: Celine Hernandez
# Contact: wwwpaf@unil.ch
#
#############################################################################################

## Read command line arguments

options <- commandArgs(trailingOnly = TRUE)
if(length(options) != 2) 
    stop('\n\nUsage:\n  Rscript preprocessing_MQproteinGroups.R <input_file_name> <output_file_name>\n\nThis script needs two arguments.\nSee internal documentation for more details and examples of command lines.\n\n')
## location of the MaxQuant ProteinGroups.txt
inputFileName <- options[1]
## output file
outputFileName <- options[2]

## Open input file and load all data
dataMQ <- read.table(inputFileName, stringsAsFactors=FALSE, quote="", na='',
                     check.names=FALSE,
                     row.names=NULL, header=TRUE, sep="\t", fill=TRUE, comment.char="")


################################################################################ 
## Preliminary notes
# 
# There was a lot of manual treatment on the original MaxQuant proteinGroup.txt files 
# and between step 1 and 2, before obtaining the final images.
# Other contaminants can be found in data/potential_skin_contaminants.xls, were they used?
# 
# From the information I could gather, I was not able to reproduce this manual processing, unfortunately.
# Informations are stored below, as well as my attempts in comments (coherent with the output 
# I could find), but as soon as I try to go further I get discrepancies 
# with my reference file...
# 
# Note to self: NO manual processing! This is not good for reproducibility!
# 
# 
# 4045
# Data were analysed with MaxQuant 1.0.13.13 against the IPI_human database 
# with standard parameters ; the protein_groups table was then processed further through the following steps
# 1) remove reverse database hits
# 2) manually annotate presumable contaminants from skin proteins based on 
#  following criteria : 
#   I) preexisting list from other SILAC experiments 
#   ii) SILAC ratio very low since they are not isotope labelled 
#   iii) found in several fractions all along the gel
# 3) eliminate contaminants identified in step 2) and protein groups with ratio count =0
# 4) copy worksheet and delete proteins with H/L ratio greater than 1.2 or 
#   lower than 0.8=> “less than 0.8” and “more than 1.2” sheets ; 
#   for all these delete proteins with ratio count <3 ; 
#   delete band info to make data simpler
# 5)  main table : run through “gel_mobility” script (pasted below) ; 
#   only protein groups kept that had both one band with ratio<0.9 AND one band 
#   with ratio >1.1. All these with a minimum of 2 peptides for confidence 
#   of quantification
# 
# 4047
# Data were analysed with MaxQuant 1.0.13.13 against the IPI_human database 
# with standard parameters ; the protein_groups table was then processed further through the following steps
# 1) remove reverse database hits
# 2) manually annotate presumable contaminants from skin proteins based on 
#  following criteria : 
#   I) preexisting list from other SILAC experiments 
#   ii) SILAC ratio very low since they are not isotope labelled 
#   iii) found in several fractions all along the gel
# 3)  eliminate contaminants identified in step 2) and protein groups with ratio count =0
# 4) copy worksheet and delete proteins with H/L ratio greater than 0.5 or 
#   lower than 2.0=> “ more2.0_less0.5 “ sheet ; 
#   for all these delete proteins with ratio count <3 ; 
#   delete band info to make data simpler
# 5)  main table : run through “gel_mobility” script (pasted below) ; 
#   only protein groups kept that had both one band with ratio<1 AND one band 
#   with ratio >2.0. All these with a minimum of 2 peptides for confidence 
#   of quantification
################################################################################ 
# 
# ## 1 Reverse
# reverseName <- "Reverse"
# ## 2 Contaminants
# silacRatioName <- "Ratio.H.L.Normalized"
# contaminantName <- "Contaminant"
# 
# bandCountNames <- grep(x=names(dataMQ), pattern="Ratio.H.L.Count.band[0-9]{1,2}", perl=TRUE, value=TRUE)
# 
# filter <- apply(dataMQ, 1, 
#                 FUN=function(xVals){
#                     reverseOK <- is.na(xVals[reverseName]) | xVals[reverseName] == ''
#                     
#                     silacRatioNotLowOK <- as.numeric(xVals[silacRatioName]) > 0.03 | is.na(xVals[silacRatioName])
#                     silacRatioNotHighOK <- as.numeric(xVals[silacRatioName]) < 12.0 | is.na(xVals[silacRatioName])
#                     contaminantOK <- is.na(xVals[contaminantName]) | xVals[contaminantName] == ''
#                     
#                     return(reverseOK & silacRatioNotLowOK & silacRatioNotHighOK & contaminantOK)
#                 })





################################################################################ 
## A sound pre-processing ?
# 
# Consists in removing Reverse and Contaminants, as well as protein groups 
# with less than 3 identified peptides, and appearing in one band only.
# 

## Contaminants, Reverse
contaminantsNames <- c("Contaminant", "Reverse")
## Only identified by site
onlyIdBySiteName <- "Only identified by site"
## At least 3 peptides identified for all bands
allCountName <- "Ratio H/L Count"
minPeptidesCount <- 3
## At least 2 bands
bandCountNames <- grep(x=names(dataMQ), pattern="Ratio H/L Count band[0-9]{1,2}", perl=TRUE, value=TRUE)


## Create filter

filter <- apply(dataMQ, 1, 
                FUN=function(xVals){
                    counts <- as.numeric(xVals[allCountName])
                    countOK <- (counts >= minPeptidesCount & !is.na(counts))
                    
                    onlyIdBySiteOK <- xVals[onlyIdBySiteName] == "False"
                    contaminantOK <- all(is.na(xVals[contaminantsNames]) | xVals[contaminantsNames] == '')
                    
                    bandNbOK <- length(which(xVals[bandCountNames] != 0)) >= 2
                    return(contaminantOK & onlyIdBySiteOK & countOK & bandNbOK)
                })




## Write filtered table
write.table(x=dataMQ[filter, ], file=outputFileName, sep="\t", na='', quote=FALSE, row.names=FALSE)

