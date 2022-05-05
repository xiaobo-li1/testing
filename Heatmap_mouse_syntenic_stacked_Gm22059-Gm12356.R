# Counts

library(edgeR)
library("pheatmap")
library("RColorBrewer")
library("ComplexHeatmap")
library(circlize)

set.seed(123)

rawcounts <- read.table("Mouse_read_counts_heatmap_seekr.txt",sep="\t",header=TRUE)
rownames(rawcounts) <- rawcounts[,1]


rnaseqMatrix <- round(rawcounts[,c(2:19)]) + 1

samples <- c("HBS00_1","HBS00_2","HBS00_3","HBS00_4",
             "IAV00_1","IAV00_2","IAV00_3","IAV00_4",
             "HBS09_1","HBS09_2","HBS09_3","HBS09_4",
             "IAV09_1","IAV09_2","IAV09_3","IAV09_4",
             "HBS18_1","HBS18_2","HBS18_3","HBS18_4",
             "IAV18_1","IAV18_2","IAV18_3","IAV18_4",
             "HBS21_1","HBS21_2","HBS21_3","HBS21_4",
             "IAV21_1","IAV21_2","IAV21_3","IAV21_4",                               
             "HBS24_1","HBS24_2","HBS24_3","HBS24_4",
             "IAV24_1","IAV24_2","IAV24_3","IAV24_4",                                     
             "HBS27_1","HBS27_2","HBS27_3","HBS27_4",
             "IAV27_1","IAV27_2","IAV27_3","IAV27_4",                                
             "HBS30_1","HBS30_2","HBS30_3","HBS30_4",
             "IAV30_1","IAV30_2","IAV30_3","IAV30_4",                                  
             "HBS36_1","HBS36_2","HBS36_3","HBS36_4",
             "IAV36_1","IAV36_2","IAV36_3","IAV36_4")

genes <- rownames(rnaseqMatrix)

# Calculate fold changes for candidate gene list
#aggregate( CPM ~ Group, boxplot.tnfa, mean)
Group <- c("Mouse_resection_12hr","Mouse_resection_12hr","Mouse_resection_12hr",
           "Mouse_resection_24hr","Mouse_resection_24hr","Mouse_resection_24hr",
           "Mouse_resection_48hr","Mouse_resection_48hr","Mouse_resection_48hr",
           "Mouse_sham_12hr","Mouse_sham_12hr","Mouse_sham_12hr",
           "Mouse_sham_24hr","Mouse_sham_24hr","Mouse_sham_24hr",
           "Mouse_sham_48hr","Mouse_sham_48hr","Mouse_sham_48hr")

design <- model.matrix(~0+Group)

colnames(design) <- levels(Group)

dge <- DGEList(counts=rnaseqMatrix,group=Group)

dge <- calcNormFactors(dge)

dge <- estimateGLMTrendedDisp(dge,verbose=TRUE)


my.contrasts <- makeContrasts(
  Twelve = Mouse_resection_12hr-Mouse_sham_12hr,
  #Twelve_02 = Mouse_resection_12hr_2-Mouse_sham_12hr_2,
  #Twelve_03 = Mouse_resection_12hr_3-Mouse_sham_12hr_3,
  Twentyfour = Mouse_resection_24hr-Mouse_sham_24hr,
  #Twentyfour_02 = Mouse_resection_24hr_2-Mouse_sham_24hr_2,
  #Twentyfour_03 = Mouse_resection_24hr_3-Mouse_sham_24hr_3,
  Fortyeight = Mouse_resection_48hr-Mouse_sham_48hr,
  #Fortyeight_02 = Mouse_resection_48hr_2-Mouse_sham_48hr_2,
  #Fortyeight_03 = Mouse_resection_48hr_3-Mouse_sham_48hr_3,
  levels=levels(factor(Group)))

fit <- glmFit(dge, design)

lrt_Twelve <- glmLRT(fit, contrast=my.contrasts[,"Twelve"])
lrt_Twentyfour <- glmLRT(fit, contrast=my.contrasts[,"Twentyfour"])
lrt_Fortyeight <- glmLRT(fit, contrast=my.contrasts[,"Fortyeight"])
# lrt_INF_21 <- glmLRT(fit, contrast=my.contrasts[,"INF_21"])
# lrt_INF_24 <- glmLRT(fit, contrast=my.contrasts[,"INF_24"])
# lrt_INF_27 <- glmLRT(fit, contrast=my.contrasts[,"INF_27"])
# lrt_INF_30 <- glmLRT(fit, contrast=my.contrasts[,"INF_30"])
# lrt_INF_36 <- glmLRT(fit, contrast=my.contrasts[,"INF_36"])


# generation of the first heatmap for chromosome 3
# Candidate genes
#candidates <- read.table("nyu_unique_genes_sorted.txt")
candidates <- read.table("Heatmap_mouse_synetic_rank_Gm22059-Gm12356.txt", header = TRUE)
colnames(candidates)[1] <- "Gene.ID"
candidates[,3] = candidates[,3]/10000
candidates_chr3 <- candidates[1:34,]

# names(candidates) <- c("ID","Zfish_symbol","Human_symbol","Screen1","Screen2","Sum","RankDiff","Rank1","Rank2","RankProd","Rank")
#order
candidates_chr3$order <- c(1:34)
table_Twelve <- topTags(lrt_Twelve, n=Inf)
table_Twentyfour <- topTags(lrt_Twentyfour, n=Inf)
table_Fortyeight <- topTags(lrt_Fortyeight, n=Inf)
# table_INF_21 <- topTags(lrt_INF_21, n=Inf)
# table_INF_24 <- topTags(lrt_INF_24, n=Inf)
# table_INF_27 <- topTags(lrt_INF_27, n=Inf)
# table_INF_30 <- topTags(lrt_INF_30, n=Inf)
# table_INF_36 <- topTags(lrt_INF_36, n=Inf)
candidates_chr3_Twelve <- merge(candidates_chr3,table_Twelve,by.x=1,by.y=0,all.x=TRUE)
candidates_chr3_Twentyfour <- merge(candidates_chr3,table_Twentyfour,by.x=1,by.y=0,all.x=TRUE)
candidates_chr3_Fortyeight <- merge(candidates_chr3,table_Fortyeight,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_21 <- merge(candidates,table_INF_21,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_24 <- merge(candidates,table_INF_24,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_27 <- merge(candidates,table_INF_27,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_30 <- merge(candidates,table_INF_30,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_36 <- merge(candidates,table_INF_36,by.x=1,by.y=0,all.x=TRUE)




FC_for_heatmap_chr3 <- as.data.frame(cbind(candidates_chr3_Twelve$logFC,candidates_chr3_Twentyfour$logFC,candidates_chr3_Fortyeight$logFC)) 
names(FC_for_heatmap_chr3) <- c("12hr","24hr","48hr")
row.names(FC_for_heatmap_chr3) <- candidates_chr3_Twelve$Gene.ID
FC_for_heatmap_chr3_labeled <- merge(candidates_chr3,FC_for_heatmap_chr3,by.x=1,by.y=0)
FC_for_heatmap_chr3_labeled <- FC_for_heatmap_chr3_labeled[order(FC_for_heatmap_chr3_labeled$order),]
FC_for_heatmap_chr3_labeled <- FC_for_heatmap_chr3_labeled[,-4]
distance <- as.matrix(FC_for_heatmap_chr3_labeled$Distance)
rownames(distance) <- FC_for_heatmap_chr3_labeled$Gene_name
FC_for_heatmap_chr3 <- FC_for_heatmap_chr3_labeled[,c(4:6)]
row.names(FC_for_heatmap_chr3) <- FC_for_heatmap_chr3_labeled[,2]

FC_for_heatmap_chr3[is.na(FC_for_heatmap_chr3)] <- 0



# row_ha = rowAnnotation(Rank = anno_barplot(distance))
# row_ha = rowAnnotation(Distance = anno_barplot(distance, axis = FALSE, axis_param = list(side = "top")))
row_ha = rowAnnotation(' ' = anno_barplot(distance, axis = FALSE, ylim = range(-50, 50)), width=unit(2, "cm"))

# pdf("FC_heatmap_chr3.pdf",width=6,height=6) # change height to 12
Heatmap(as.matrix(FC_for_heatmap_chr3),
        column_names_max_height = unit(10, "cm"), row_names_max_width = unit(6, "cm"), height = nrow(as.matrix(FC_for_heatmap_chr3))*unit(3, "mm"),
        cluster_rows = FALSE,
        column_names_gp = gpar(fontsize = 14),row_names_gp = gpar(fontsize = 8),
        show_heatmap_legend = TRUE,row_dend_reorder = FALSE,
        right_annotation = row_ha,
        heatmap_legend_param = list(at = c(-1.5, -1.0, -0.50, 0, 0.50, 1.0, 1.5, 2.0, 2.5, 3.00), 
                                    labels = c("-1.50", "-1.00", "-0.50", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "3.00")),
        cluster_columns = FALSE,
        col = colorRamp2(c(-1.7, 0, 2.5), c("blue","white", "red")),
        #       column_hclust_height=unit(10,"mm"),
        name = "log2(Fold Change)",
        gap=unit(0,"mm")
)
dev.off()
ht1 = Heatmap(as.matrix(FC_for_heatmap_chr3),
              column_names_max_height = unit(10, "cm"), row_names_max_width = unit(6, "cm"), cluster_rows = FALSE,
              column_names_gp = gpar(fontsize = 14),row_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = TRUE, row_dend_reorder = FALSE,
              right_annotation = row_ha,
              heatmap_legend_param = list(at = c(-1.5, -1.0, -0.50, 0, 0.50, 1.0, 1.5, 2.0, 2.5, 3.0), 
                                          labels = c("-1.50", "-1.00", "-0.50", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "3.00")),
              cluster_columns = FALSE,
              col = colorRamp2(c(-1.7, 0, 2.5), c("blue","white", "red")),
              #       column_hclust_height=unit(10,"mm"),
              name = "log2(Fold Change)",
              gap=unit(0,"mm"))


# generation of the second heatmap for chromosome 7
# Candidate genes
#candidates <- read.table("nyu_unique_genes_sorted.txt")
# candidates <- read.table("Heatmap_mouse_seekr_rank.txt", header = TRUE)
# colnames(candidates)[1] <- "Gene.ID"
candidates_chr12 <- candidates[35:68,]
# order
candidates_chr12$order <- c(1:34)
# names(candidates) <- c("ID","Zfish_symbol","Human_symbol","Screen1","Screen2","Sum","RankDiff","Rank1","Rank2","RankProd","Rank")

table_Twelve <- topTags(lrt_Twelve, n=Inf)
table_Twentyfour <- topTags(lrt_Twentyfour, n=Inf)
table_Fortyeight <- topTags(lrt_Fortyeight, n=Inf)
# table_INF_21 <- topTags(lrt_INF_21, n=Inf)
# table_INF_24 <- topTags(lrt_INF_24, n=Inf)
# table_INF_27 <- topTags(lrt_INF_27, n=Inf)
# table_INF_30 <- topTags(lrt_INF_30, n=Inf)
# table_INF_36 <- topTags(lrt_INF_36, n=Inf)
candidates_chr12_Twelve <- merge(candidates_chr12,table_Twelve,by.x=1,by.y=0,all.x=TRUE)
candidates_chr12_Twentyfour <- merge(candidates_chr12,table_Twentyfour,by.x=1,by.y=0,all.x=TRUE)
candidates_chr12_Fortyeight <- merge(candidates_chr12,table_Fortyeight,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_21 <- merge(candidates,table_INF_21,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_24 <- merge(candidates,table_INF_24,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_27 <- merge(candidates,table_INF_27,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_30 <- merge(candidates,table_INF_30,by.x=1,by.y=0,all.x=TRUE)
# candidates_INF_36 <- merge(candidates,table_INF_36,by.x=1,by.y=0,all.x=TRUE)

FC_for_heatmap_chr12 <- as.data.frame(cbind(candidates_chr12_Twelve$logFC,candidates_chr12_Twentyfour$logFC,candidates_chr12_Fortyeight$logFC)) 
names(FC_for_heatmap_chr12) <- c("12hr","24hr","48hr")
row.names(FC_for_heatmap_chr12) <- candidates_chr12_Twelve$Gene.ID
FC_for_heatmap_chr12_labeled <- merge(candidates_chr12,FC_for_heatmap_chr12,by.x=1,by.y=0)
FC_for_heatmap_chr12_labeled <- FC_for_heatmap_chr12_labeled[order(FC_for_heatmap_chr12_labeled$order),]
FC_for_heatmap_chr12_labeled <- FC_for_heatmap_chr12_labeled[,-4]
distance_chr12 <- as.matrix(FC_for_heatmap_chr12_labeled$Distance)
rownames(distance_chr12) <- FC_for_heatmap_chr12_labeled$Gene_name
FC_for_heatmap_chr12 <- FC_for_heatmap_chr12_labeled[,c(4:6)]
row.names(FC_for_heatmap_chr12) <- FC_for_heatmap_chr12_labeled[,2]

FC_for_heatmap_chr12[is.na(FC_for_heatmap_chr12)] <- 0



row_ha_chr12 = rowAnnotation("Distance 10kb" = anno_barplot(distance_chr12, ylim = range(-50, 50)), width=unit(2, "cm"))

# pdf("FC_heatmap_chr12.pdf",width=6,height=6)
Heatmap(as.matrix(FC_for_heatmap_chr12),
        column_names_max_height = unit(10, "cm"), row_names_max_width = unit(6, "cm"), cluster_rows = FALSE,
        column_names_gp = gpar(fontsize = 14),row_names_gp = gpar(fontsize = 8),
        show_heatmap_legend = TRUE,row_dend_reorder = FALSE,
        right_annotation = row_ha_chr12,
        heatmap_legend_param = list(at = c(-1.5, -1.0, -0.50, 0, 0.50, 1.0, 1.5, 2.0, 2.5, 3.00), 
                                    labels = c("-1.50", "-1.00", "-0.50", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "3.00")),
        cluster_columns = FALSE,
        col = colorRamp2(c(-1.7, 0, 2.5), c("blue","white", "red")), 
        #       column_hclust_height=unit(10,"mm"),
        name = "log2(Fold Change)",
        gap=unit(0,"mm")
)

dev.off()
ht2 = Heatmap(as.matrix(FC_for_heatmap_chr12),
              column_names_max_height = unit(10, "cm"), row_names_max_width = unit(6, "cm"), cluster_rows = FALSE,
              column_names_gp = gpar(fontsize = 14),row_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = TRUE, row_dend_reorder = FALSE,
              right_annotation = row_ha_chr12,
              heatmap_legend_param = list(at = c(-1.5, -1.0, -0.50, 0, 0.50, 1.0, 1.5, 2.0, 2.5, 3.00), 
                                          labels = c("-1.50", "-1.00", "-0.50", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "3.00")),
              cluster_columns = FALSE,
              col = colorRamp2(c(-1.7, 0, 2.5), c("blue","white", "red")),
              #       column_hclust_height=unit(10,"mm"),
              name = "log2(Fold Change)",
              gap=unit(0,"mm")
)

ht_list = ht1 %v% ht2

pdf("Mouse_syntenic_heatmap_Gm22059-Gm12356.pdf",width=6,height=12)
draw(ht_list)

dev.off()
