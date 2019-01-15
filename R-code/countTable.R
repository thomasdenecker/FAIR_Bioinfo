################################################################################
# Boxplot
# Thomas Denecker
# Janvier 2019
################################################################################

#===============================================================================
# Library
#===============================================================================

library(dplyr)
library(tidyverse)

#===============================================================================
# Function
#===============================================================================

#===============================================================================
# Main
#===============================================================================

conditions <- read.csv2("/home/rstudio/conditions.txt", sep ="\t", header = T)

listDFCondA <- list()
listDFCondB <- list()

for (i in 1:nrow(conditions)) {
  if(conditions[i,"Condition_Name"] == "CondA"){
    listDFCondA[[as.character(conditions[i,1])]] = read.table(paste0("/home/rstudio/Project/htseq/count_",as.character(conditions[i,1]) , ".txt"), 
                                                              header = F, sep ="\t", col.names = c('Features', paste0("CondA_",as.character(conditions[i,1]))))
  } else {
    listDFCondB[[as.character(conditions[i,1])]] = read.table(paste0("/home/rstudio/Project/htseq/count_",as.character(conditions[i,1]) , ".txt"), 
                                                              header = F, sep ="\t", col.names = c('Features', paste0("CondB_",as.character(conditions[i,1]))))
  }
}

countTableConA = listDFCondA %>% reduce(full_join, by = "Features")
countTableConB = listDFCondB %>% reduce(full_join, by = "Features")

countTable = full_join(countTableConA, countTableConB, by= "Features")
rownames(countTable) = countTable[,1]
countTable = countTable[,-1]
countTable = data.matrix(countTable)
countTable = countTable[ -which(rownames(countTable) %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")), ]

write.table(countTable, "countTable.txt", quote = F, sep = "\t")