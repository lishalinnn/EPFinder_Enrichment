
#Setp1: merged in gwas results in BMD results-----

library(data.table)
library(readxl)
library(writexl)
library(dplyr)
library(cmapR)

bmd_result <- fread("/mnt/Storage2/lishalin/EPFinder_Enrichment/BMD_EPFinder_positive_prediction.tsv")
colnames(bmd_result)[4]<-"Target_Gene"
bmd_result1 <- bmd_result[,c(1,2,3,5,6,8,9)]
bmd_result2 <- bmd_result1[,c(4,6,7)]
bmd_result2$source <- "BMD"



#Do the same thing for eBMD.
ebmd <- read.csv("/mnt/Storage2/lishalin/EPFinder_Enrichment/eBMD_EPFinder_positive_prediction.csv")
ebmd_1 <- ebmd[,c(2,3,5)]
ebmd_1$source <- "eBMD"

#Now rename bmd and ebmd columns to make sure they match for the sake of merging.
colnames(bmd_result2)<- c("SNP_id","Prediction_Gene","Prediction_Score","Source")
colnames(ebmd_1)<- c("SNP_id","Prediction_Score","Prediction_Gene","Source")

#
filtered_bmd <- bmd_result2 %>%
  group_by(SNP_id, Prediction_Gene) %>%
  slice_max(order_by = Prediction_Score, n = 1, with_ties = FALSE) %>%
  ungroup()

#
filtered_ebmd <- ebmd_1 %>%
  group_by(SNP_id, Prediction_Gene) %>%
  slice_max(order_by = Prediction_Score, n = 1, with_ties = FALSE) %>%
  ungroup()

#Rowbind them together into one dataset. 
combined_data <- rbind(filtered_bmd,filtered_ebmd)

#Do one more filter to see whether there are same pair shows up in both BMD and eBMD results. 
combined_clean <- combined_data %>%
  group_by(SNP_id, Prediction_Gene) %>%
  summarise(
    n_sources = n_distinct(Source),
    Combined_Source = ifelse(n_sources > 1, paste(sort(unique(Source)), collapse = "; "), unique(Source)),
    Prediction_Score = max(Prediction_Score),
    Prediction_Score_Source = Source[which.max(Prediction_Score)],
    .groups = 'drop'
  )

#No overlap. Now can write out the combined_data file. 
write.csv(combined_data, "/mnt/Storage2/lishalin/EPFinder_Enrichment/Prediction_Gene_Combined_05192025.csv")


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




 