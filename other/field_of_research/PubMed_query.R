library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
#--------------------------------------------------
AllResearch <- "all[sb]"
MainField <- "(obesity OR diabetes)"
SubFieldA <- "inflammation"
SubFieldB <- "exercise"
SubFieldC <- '"skeletal muscle"'

#install.packages("RISmed")
library(RISmed)

#all publications on pubmed
count_AllResearch <- 31821468 #EUtilsSummary(AllResearch, type="esearch", db="pubmed")@count
Sys.sleep(1) 
count_MainField <- EUtilsSummary(MainField, type="esearch", db="pubmed")@count
Sys.sleep(1) 
count_SubFieldA <- EUtilsSummary(paste(MainField, SubFieldA, sep=" AND "), type="esearch", db="pubmed")@count
Sys.sleep(1) 
count_SubFieldB <- EUtilsSummary(paste(MainField, SubFieldB, sep=" AND "), type="esearch", db="pubmed")@count
Sys.sleep(1) 
count_SubFieldC <- EUtilsSummary(paste(MainField, SubFieldC, sep=" AND "), type="esearch", db="pubmed")@count
Sys.sleep(1) 
count_SubFieldAB <- EUtilsSummary(paste(MainField, SubFieldA, SubFieldB, sep=" AND "), type="esearch", db="pubmed")@count
Sys.sleep(1) 
count_SubFieldAC <- EUtilsSummary(paste(MainField, SubFieldA, SubFieldC, sep=" AND "), type="esearch", db="pubmed")@count
Sys.sleep(1) 
count_SubFieldBC <- EUtilsSummary(paste(MainField, SubFieldB, SubFieldC, sep=" AND "), type="esearch", db="pubmed")@count
Sys.sleep(1) 
count_SubFieldABC <- EUtilsSummary(paste(MainField, SubFieldA, SubFieldB, SubFieldC, sep=" AND "), type="esearch", db="pubmed")@count


library(eulerr)
VennDiag <- euler(c("Total" = (count_AllResearch)^(1/1.2),
                    "Total&Main" = (count_MainField)^(1/1.2),
                    "Total&Main&A" = (count_SubFieldA)^(1/1.2),
                    "Total&Main&B" = (count_SubFieldB)^(1/1.2), 
                    "Total&Main&C" = (count_SubFieldC)^(1/1.2),
                    "Total&Main&A&B" = (count_SubFieldAB)^(1/1.2), 
                    "Total&Main&A&C" = (count_SubFieldAC)^(1/1.2),
                    "Total&Main&B&C" = (count_SubFieldBC)^(1/1.2),
                    "Total&Main&A&B&C" = (count_SubFieldABC)^(1/1.2)
                    ),
                  shape = "ellipse")
Fields <- c("All publications in biology and medicine", MainField, SubFieldA, SubFieldB, SubFieldC)
Fields <- gsub("\\(|\\)", "", Fields)
Fields <- gsub('\\"', "", Fields)
Fields <- gsub('OR', "&", Fields)
Fields <- gsub(' ', "\n", Fields)
Fields

png(file="venn.png", res=600,
    height=20, width = 20, unit="cm")
plot(VennDiag, 
     alpha=0.5, 
     edges = list(col=c("black", "#0072B2", "#E69F00", "#009E73", "black", "grey", "grey", "gray", "yellow"),
                  lex=1),
     fill=c("white", "#56B4E9", "#E69F00", "#009E73", "black", "#739f3a", "#d98c54", "#668c8d", "yellow"),
     labels = list(labels = Fields,
                   font = 1,
                   pos=4,
                   cex = c(1,1,0.5,0.5,0.5)))
dev.off()



svg(file="venn.svg",
    height=8, width = 8)
plot(VennDiag, 
     alpha=0.5, 
     edges = list(col=c("black", "#0072B2", "#E69F00", "#009E73", "black", "grey", "grey", "gray", "yellow"),
                  lex=1),
     fill=c("white", "#56B4E9", "#E69F00", "#009E73", "black", "#739f3a", "#d98c54", "#668c8d", "yellow"),
     labels = list(labels = Fields,
                   font = 1,
                   pos=4,
                   cex = c(1,1,0.5,0.5,0.5)))
dev.off()
