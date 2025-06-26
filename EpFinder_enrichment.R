
#Setp1: merged in gwas results in BMD results-----

library(data.table)
library(readxl)
library(writexl)
library(dplyr)
library(cmapR)

#read in data
bmd <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/BMD_EPFinder_prediction.tsv")
ebmd <- read.csv("/mnt/Storage2/lishalin/EPFinder_Enrichment/cojo_bf3_EPFinder_final.csv")

ebmd_need <- ebmd[,c(5,4,2,8, 10,11,12,13)]
bmd_need <- bmd[,c(5,7,8)]

#For BMD, the column "SNPID hg38" looks like this: 1:8387662:T:C
#seperate it into four columns: chr, position, ref, alt.
bmd_need[, c("chr", "position", "ref", "alt") := tstrsplit(`SNPID hg38`, ":", fixed = TRUE)]

#There are some rows that are duplicate except for column 3
#that are duplicate in ebmd_need. 
#For these rows, only keep the row with the 
#highest value in column 3
# Remove duplicates except for column 3, keeping the row with the highest value in column 3
setDT(ebmd_need)
# ...existing code...
library(dplyr)

ebmd_need1 <- ebmd_need %>%
  group_by(across(1:2)) %>%
  slice_max(order_by = .[[3]], n = 1, with_ties = FALSE) %>%
  ungroup()

  # ...existing code...
ebmd_need2 <- ebmd_need %>%
  group_by(SNPID, Predicted_Gene) %>%
  slice_max(order_by = Prediction, n = 1, with_ties = FALSE) %>%
  ungroup()


#Define locus flanking
ebmd_need2$pos_min <- ebmd_need2$Position - 50000
ebmd_need2$pos_max <- ebmd_need2$Position + 50000

#add SNP resource
ebmd_need2$Existence  <- "eBMD"

#Now look at BMD_need data. Look at the Column chr and position. Compare it with the ebmd_need1 data. 
#If same chr number, and the position falls within the pos_min and pos_max range, then: 
#Assign the L.Bin from ebmd_need2 to bmd_need1 (create a new column called L.Bin in bmd_need1).
#If the positions fall into multiple flanking regions of the L.Bin in ebmd_need2
#Then append the LBin numbers by;, and add that character to the L.Bin in bmd_need1.

# Assume bmd_need1 is a copy of bmd_need with unique rows
bmd_need1 <- copy(bmd_need)

# Make sure columns are numeric for comparison
bmd_need1$position <- as.numeric(bmd_need1$position)
ebmd_need2$Position <- as.numeric(ebmd_need2$Position)
ebmd_need2$pos_min <- as.numeric(ebmd_need2$pos_min)
ebmd_need2$pos_max <- as.numeric(ebmd_need2$pos_max)

# Assign a unique L.Bin to each row in ebmd_need2 if not already present
if (!"L.Bin" %in% names(ebmd_need2)) {
  ebmd_need2$L.Bin <- seq_len(nrow(ebmd_need2))
}

# For each row in bmd_need1, find all L.Bin(s) in ebmd_need2 where chr matches and position is within pos_min/pos_max
library(dplyr)

bmd_need1 <- bmd_need1 %>%
  rowwise() %>%
  mutate(
    L.Bin = paste(
      ebmd_need2$L.Bin[
        ebmd_need2$chr == chr &
        position >= ebmd_need2$pos_min &
        position <= ebmd_need2$pos_max
      ],
      collapse = ";"
    )
  ) %>%
  ungroup()

#Rearrange bmd_need_unique_new columns so that the format matches the ebmd_need2
ebmd_for_merge <- ebmd_need2[,-c(9:10)]

bmd_for_merge <- bmd_need_unique_new
bmd_for_merge$Existence <- "BMD"
colnames(bmd_for_merge)<- c("SNPID","Predicted_Gene","Prediction","chr","Position","Ref","Alt","L.BIN","Existence")

#Rearrange column sequence of bmd_for_merge to the sequence in ebmd_for_merge
colnames(ebmd_for_merge)
bmd_for_merge <- bmd_for_merge[,c("SNPID","Predicted_Gene","Prediction","L.BIN","chr","Position","Ref","Alt","Existence")]


#Now merge the two datasets together.
combined_data <- rbind(ebmd_for_merge, bmd_for_merge)

#Now start the filtering process.
#For rows that contains the same number for L.BIN, Keep the row with the highest value in column 3 (Prediction).
#Add a column called Source, and assign the value "eBMD" or "BMD" to the rows accordingly.
#Also if this process happens, then change the value in Existence to "eBMD; BMD".
# ...existing code...

# Now start the filtering process.
# For rows that contain the same number for L.BIN, keep the row with the highest value in column 3 (Prediction).
# Add a column called Source, and assign the value "eBMD" or "BMD" to the rows accordingly.
# If a row is present in both, change Existence to "eBMD; BMD".

# ...existing code...
combined_data1 <- combined_data %>%
  group_by(L.BIN) %>%
  mutate(
    is_empty = is.na(L.BIN) | L.BIN == ""
  ) %>%
  filter(
    is_empty | Prediction == max(Prediction, na.rm = TRUE)
  ) %>%
  mutate(
    Source = ifelse(is_empty, NA, Existence),
    Existence = ifelse(is_empty, Existence, "eBMD; BMD")
  ) %>%
  ungroup() %>%
  select(-is_empty)


combined_data1$Source<-ifelse(is.na(combined_data1$Source), combined_data1$Existence, combined_data1$Source)

write.csv(combined_data1, "/mnt/Storage2/lishalin/EPFinder_Enrichment/Combine_EPFinder_Results_06262025.csv",row.names = F)
#Now we have the combined data, we can start the enrichment analysis.
subset_combine <- combined_data1[,c(2,3)]

subset_combine1 <- subset_combine %>%
  arrange(desc(Prediction))

combined_data1_max <- subset_combine1 %>%
  group_by(Prediction_Gene) %>%
  slice_max(order_by = Prediction, n = 1, with_ties = FALSE) %>%
  ungroup()

write.table(combined_data1_max, file = "/mnt/Storage2/lishalin/EPFinder_Enrichment/Gene_for_enrichment_Ranked_06262025.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#Use GSEA to do the enrichment analysis
#transform all tsv file to excel file.
library(data.table)
#H
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/New_GSEAPreRank_Results/gsea_report_for_na_pos_H.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/New_GSEAPreRank_Results/gsea_report_for_na_neg_H.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_H_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_H_neg.xlsx")


#C2
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_C2.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_C2.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C2_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C2_neg.xlsx")


#C3
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_C3.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_C3.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C3_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C3_neg.xlsx")

#C5
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_C5.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_C5.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C5_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C5_neg.xlsx")


#C8
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_C8.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_C8.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C8_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C8_neg.xlsx")






















#Now get the list of genes
gene_set <- unique(combined_data$Prediction_Gene)
write_grp(gene_set, "/mnt/Storage2/lishalin/EPFinder_Enrichment/Gene_for_enrichment.grp")

library(GSEABase)
library(GSEA)
library(clusterProfiler)

#Use MSigDb
h_set <- read.gmt("/mnt/Storage2/lishalin/EPFinder_Enrichment/h.all.v2024.1.Hs.symbols.gmt")
c2_set <- read.gmt("/mnt/Storage2/lishalin/EPFinder_Enrichment/c2.all.v2024.1.Hs.symbols.gmt")
c3_set <- read.gmt("/mnt/Storage2/lishalin/EPFinder_Enrichment/c3.all.v2024.1.Hs.symbols.gmt")
c5_set <- read.gmt("/mnt/Storage2/lishalin/EPFinder_Enrichment/c5.all.v2024.1.Hs.symbols.gmt")
c8_set <- read.gmt("/mnt/Storage2/lishalin/EPFinder_Enrichment/c8.all.v2024.1.Hs.symbols.gmt")

h_result <- as.data.frame(enricher(gene = gene_set, TERM2GENE = h_set, pvalueCutoff = 0.05, pAdjustMethod = "fdr"))
c2_result <- as.data.frame(enricher(gene = gene_set, TERM2GENE = c2_set, pvalueCutoff = 0.05, pAdjustMethod = "fdr"))
c3_result <- as.data.frame(enricher(gene = gene_set, TERM2GENE = c3_set, pvalueCutoff = 0.05, pAdjustMethod = "fdr"))
c5_result <- as.data.frame(enricher(gene = gene_set, TERM2GENE = c5_set, pvalueCutoff = 0.05, pAdjustMethod = "fdr"))
c8_result <- as.data.frame(enricher(gene = gene_set, TERM2GENE = c8_set, pvalueCutoff = 0.05, pAdjustMethod = "fdr"))

write.csv(h_result, "/mnt/Storage2/lishalin/EPFinder_Enrichment/NEW_H_result.csv")
write.csv(c2_result, "/mnt/Storage2/lishalin/EPFinder_Enrichment/NEW_C2_result.csv")
write.csv(c3_result, "/mnt/Storage2/lishalin/EPFinder_Enrichment/NEW_C3_result.csv")
write.csv(c5_result, "/mnt/Storage2/lishalin/EPFinder_Enrichment/NEW_C5_result.csv")
write.csv(c8_result, "/mnt/Storage2/lishalin/EPFinder_Enrichment/NEW_C8_result.csv")


#Now also write out combined-data as .rnk file to be put into GSEA software for analysis
subset_combine <- combined_data[,c(2,3)]
subset_combine <- subset_combine %>%
  arrange(desc(Prediction_Score))


write.table(subset_combine, file = "/mnt/Storage2/lishalin/EPFinder_Enrichment/Gene_for_enrichment_Ranked.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)





bmd_collapsed  <-unique(bmd_result2[,c(1:2)])


targe_gene_list <- bmd_result %>%
  group_by(RSID) %>%
  reframe(Target_Gene = paste(unique(Target_Gene), collapse = ";"))
  
#Merge target gene back to the data
bmd_result_final <- merge(targe_gene_list, bmd_collapsed, by="RSID")
colnames(bmd_result_final)[6]<-"SNPID"

gwas_result <- read_excel("/mnt/Storage2/lishalin/EPFinder_Enrichment/GEFOSII_GWAS_genes_new with Zebrafish Orthologs.xlsx",sheet=2)

gwas_result1 <- gwas_result[,c(6,7,17:22)]
unique_gwas <- unique(gwas_result1)

mergetable <- merge(bmd_result_final, unique_gwas, by=c("SNPID","RSID"))
unique_mergetable <- unique(mergetable)


write_xlsx(mergetable, "/mnt/Storage2/lishalin/EPFinder_Enrichment/BMD_EPFinderResults_Merged.xlsx")

#Step2: combine the gene from BMD and eBMD-----
bmd <- read_excel("/mnt/Storage2/lishalin/EPFinder_Enrichment/BMD_EPFinderResults_Merged.xlsx")
ebmd <- read.csv("/mnt/Storage2/lishalin/EPFinder_Enrichment/eBMD_EPFinder_positive_prediction.csv")

#BMD prediction score
gene1 <- bmd[,c(7:8)]
colnames(gene1)<-c("Predicted_Gene","Prediction")

#eBMD score
gene2 <- ebmd[,c(5,3)]
colnames(gene2)<-c("Predicted_Gene","Prediction")

library(tidyverse)

gene1_max <- gene1 %>%
  group_by(Predicted_Gene) %>%
  slice_max(order_by = Prediction, n = 1, with_ties = FALSE)

gene2_max <- gene2 %>%
  group_by(Predicted_Gene) %>%
  slice_max(order_by = Prediction, n = 1, with_ties = FALSE)

combine_data <- rbind(gene1_max, gene2_max)

combine_data_max <- combine_data %>%
  group_by(Predicted_Gene) %>%
  slice_max(order_by = Prediction, n = 1, with_ties = FALSE) %>% 
  arrange(desc(Prediction))

#Step3: Write as .rnk file (tab-delimited, no row names or quotes)-----
write.table(combine_data_max, file = "/mnt/Storage2/lishalin/EPFinder_Enrichment/Gene_for_enrichment_Ranked.rnk", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


gene_set <- union(gene1, gene2)
write_grp(gene_set, "/mnt/Storage2/lishalin/EPFinder_Enrichment/Gene_for_enrichment.grp")

#Step3: GSEA analysis
gene_set <- as.list(gene_set)



#Step4
library(data.table)
h_data <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/H_over_representation.txt")
c2_data <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/C2_over_representation.txt")
c3_data <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/C3_over_representation.txt")
c5_data <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/C5_over_representation.txt")
c8_data <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/C8_over_representation.txt")

write_xlsx(h_data,"/mnt/Storage2/lishalin/EPFinder_Enrichment/H_over_representation.xlsx")
write_xlsx(c2_data,"/mnt/Storage2/lishalin/EPFinder_Enrichment/C2_over_representation.xlsx")
write_xlsx(c3_data,"/mnt/Storage2/lishalin/EPFinder_Enrichment/C3_over_representation.xlsx")
write_xlsx(c5_data,"/mnt/Storage2/lishalin/EPFinder_Enrichment/C5_over_representation.xlsx")
write_xlsx(c8_data,"/mnt/Storage2/lishalin/EPFinder_Enrichment/C8_over_representation.xlsx")


#Same thing for the GSEA preranked results
library(data.table)
#H
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_H.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_H.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_H_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_H_neg.xlsx")


#C2
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_C2.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_C2.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C2_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C2_neg.xlsx")


#C3
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_C3.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_C3.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C3_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C3_neg.xlsx")

#C5
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_C5.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_C5.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C5_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C5_neg.xlsx")


#C8
data_pos <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_pos_C8.tsv")
data_neg <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/gsea_report_for_na_neg_C8.tsv")

write_xlsx(data_pos,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C8_pos.xlsx")
write_xlsx(data_neg,"/mnt/Storage2/lishalin/EPFinder_Enrichment/GSEA_preranked_C8_neg.xlsx")




