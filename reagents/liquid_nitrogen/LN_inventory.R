#set working directory to source file location
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# Load libraries
library(openxlsx)
library(dplyr)
library(readxl)

#Gets names of excel files for LN tanks and puts them in a character vector
LN_file_names <- "../../../000_IntFys_secure/Liquid_nitrogen/LN inventory.xlsx"

#Stores excel files as workbooks into a list
LN_tank_workbooks <- list()

for (i in 1 : length(LN_file_names)){
  LN_tank_workbooks <- c(LN_tank_workbooks, loadWorkbook(LN_file_names[i]))
}

#Creates R dataframes with names based on LNtank_Rack
LN_dataframes <- list()

k <- 1

for (i in 1 : length(LN_file_names)){
  for (j in 1 : length(sheets(LN_tank_workbooks[[i]]))){
    LN_dataframes[[k]] <-
      read.xlsx(xlsxFile = LN_file_names[i],
                sheet = j,
                startRow = 1)
    names(LN_dataframes)[k] <-
      paste0(sheets(LN_tank_workbooks[[i]])[j])
    k <- k + 1
  }
}

rm(i, j, k)


#Pulls out dataframes from list of dataframes for direct manipulation with sql
for (i in 1 : length(LN_dataframes)){
  assign(x = names(LN_dataframes)[i],
         value = LN_dataframes[[i]])
}

rm(i)

columns_names <- data.frame(
  colnames(Egg13_ABC),
  colnames(Egg13_DEF),
  colnames(Egg13_GHI),
  colnames(Egg13_JKL),
  colnames(Egg13_MNOP),
  colnames(Egg14_A),
  colnames(Egg14_BC)
)

#put all tabs in one file
Egg_all <- LN_dataframes[[1]]
for (i in 2:length(LN_dataframes)){
  Egg_all <- rbind(Egg_all, LN_dataframes[[i]])
}
colnames(Egg_all) <- gsub("\\.", " ", colnames(Egg_all))

saveRDS(Egg_all, "full_inventory.Rds")

#save date of last update
last_update <- format(Sys.Date(), "%B %d, %Y")
saveRDS(last_update, "last_update.Rds")
