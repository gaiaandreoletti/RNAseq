rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire from RNAseq
###
### CITATION: 
###
### PROCESS: Diff Expression analysis uisng the RNASeq data only coding genes
###           
### DESCRIP: Analysis with the data from the QC process
###         
###
### 
### 
############################################################################################
library("VSURF")
library("randomForest")
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("DESeq2")
library("glmnet")

triage_data=read.table ('~/Documents/Phyton_Score/prova.txt',header = T)

#triage_data <- cbind(fpkm_file$Sham_discordance_score,fpkm_file$TACvehicle_discordance_score,fpkm_file$TACj_discordance_score)
#exp_data <- cbind(fpkm_file$Sham_expression_value,fpkm_file$TAC_vehicle_expression_value,fpkm_file$TACjq1_expression_value)
#rownames(exp_data) <- fpkm_file$Gene

head(triage_data)
dim(triage_data)
colnames(triage_data) <- c("SHAM", "TACveh", "TACjq1")
rownames(triage_data) <- fpkm_file$Gene
head(triage_data)
dim(triage_data)

### Limma 
# Create object for use in limma and edgeR analysis
d <- DGEList(triage_data)

# calculate normalization factors
d <- calcNormFactors(d)
d
# filter genes below X cpms in all samples 
cutoff <- 1
drop <- which(apply(cpm(d), 1, max) < cutoff)
d <- d[-drop,] 
dim(d) # number of genes left


# Group information
replicate <- c("1","2","3","1","2","3","1","2","3")
Condition <- c("TAC","TAC","TAC","TAC","TAC","TAC","SHAM","SHAM","SHAM")
treatment <- c("JQ1","JQ1","JQ1","veh","veh","veh","veh","veh","veh")

#group <- interaction(Condition, treatment)
group <- interaction(Condition)
group

# Quick MDS plot
plotMDS(d, col = as.numeric(group))

###################### Limma-voom analysis 
mm <- model.matrix(~0 + group) # specify model with no intercept for easier contrasts
y <- voom(d, mm, plot = T)
mtext(side = 3, line = 0.5, text = "1-Factor Model Without Intercept")
fit <- lmFit(y$counts, mm)
head(coef(fit))

# Comparison between cultivars C and I5 at time 6
contr <- makeContrasts(groupSHAM - groupTACjq1, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp2$adj.P.Val < 0.05)) # number of DE genes
write.table(tmp2, file = "I5_v_C_time6.txt", row.names = F, sep = "\t", quote = F)
head(tmp2, 20)
results.mod1 <- tmp2 # for illustration only

# Comparison between times 6 and 9 for cultivar I5    
contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
tmp2 <- topTable(tmp, sort.by = "P", n = Inf)
tmp2$Gene <- rownames(tmp2)
tmp2 <- tmp2[,c("Gene","logFC","AveExpr","P.Value","adj.P.Val")]
length(which(tmp2$adj.P.Val < 0.05)) # number of DE genes
write.table(tmp2, file = "time9_v_time6_I5.txt", row.names = F, sep = "\t", quote = F)
head(tmp2, 20)




# Remove all gene which has 0 value in all sample
all <- apply(triage_data, 1, function(x) all(x==0) )
newdata <- triage_data[!all,]
#write.csv(newdata, file = "count_0_filter.csv")
head (newdata)
dim(newdata)

# remove uninformative genes keep only genes that are expressed in at least 200 count in 20 samples
#triage_data <- newdata[rowSums(newdata > 3) >= 1,]
#write.csv(triage_data, file = "rpkm.csv_0_filter_atleat_200_count.csv")
#head (triage_data)
#dim(triage_data)


########################################
## Differential expression analysis ####
#######################################

###Considering only coding genes
## annotation_coding<- annotation[which(annotation$type_gene=="protein_coding"),]

###################
## 1. Dseq2 #######
###################
library(DESeq2)
#Let's make a table with the aforementioned conditions
colData <- data.frame(treatment,Condition, replicate)
#We used the header of the first column of the expression file called cts
rownames(colData) <- colnames(triage_data)
#To make sure that R see in the Replicate column a number and not a text
colData$replicate=factor(colData$replicate,levels=c("1","2","3"))
colData$all=sapply(1:nrow(colData),function(x) paste(colData$Condition[x],colData$treatment[x],sep = "."))
colData

dds=DESeqDataSetFromMatrix(countData = round(triage_data),
                           colData = colData,
                           design = ~ all)
#This is for runnung the expression analysis with a value that is normalized for different variables (e.g. lenght, this is better than fpkm)
dds <- estimateSizeFactors(dds)
dds=DESeq(dds)
summary(dds)

plotDispEsts(dds, main="Dispersion plot")

#filtering by low counts 
##keep genes that have counts per milion values above 0.5 in at least two libraries
#dds_filter <-  dds[rowSums(fpm(dds) > 0.5) >= 3] 
##normalized with rlog
rld <- rlogTransformation(dds, fitType='local')
DESeq2::plotPCA(rld, intgroup="Condition")
norm_data_rlog<-assay(normalized_rlog)
save(norm_data_rlog,file="norm_data_rlog_filter.Rdata")


mycols <- rev(brewer.pal(8, "Blues"))
sampleDists <- as.matrix(dist(t(assay(rld))))
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=mycols,
          #ColSideColors=mycols[technique], RowSideColors=mycols[technique],
          margin=c(10, 10), main="Sample Distance Matrix")


library(DEGreport)
resreport <- degResults(dds = dds, name = "test", org = NULL,
                        do_go = FALSE, group = "all", xs = "all",
                        path_results = NULL)
#load("Data/norm_data_rlog_filter.Rdata")


res1 <- results(dds, contrast=c("all","SHAM.veh", "TAC.veh")) 
length(which(res1$padj<0.05)) 
res1 <- res1[order(res1$padj), ]
head(res1)
library(biobroom)
res_tidy1 = tidy.DESeqResults(res1)
head(res_tidy1)
res_tidy1_anno = res_tidy1 %>% inner_join(grcm38, by=c("gene"="symbol")) %>% arrange(p.value) 
res_tidy1_anno_protein_coding = res_tidy1_anno %>% filter(biotype=="protein_coding")

ggplot(res_tidy1_anno_protein_coding, aes(x=estimate, y=-log(p.value),
                     color=log(baseMean))) + geom_point(alpha=0.5) +
  ggtitle("SHAM_TACveh Volcano Plot") + theme_bw()
write.csv(res_tidy1_anno_protein_coding,"SHAM_TACveh.csv",quote=F,row.names = F)



res2 <- results(dds, contrast=c("all","SHAM.veh","TAC.JQ1")) #length(which(res2$padj<0.05)) 
res2 <- res2[order(res2$padj), ]
head(res2)
library(biobroom)
res_tidy2 = tidy.DESeqResults(res2)
head(res_tidy2)
res_tidy2_anno = res_tidy2 %>% inner_join(grcm38, by=c("gene"="symbol")) %>% arrange(p.value) 
res_tidy2_anno_protein_coding = res_tidy2_anno %>% filter(biotype=="protein_coding")

ggplot(res_tidy2_anno_protein_coding, aes(x=estimate, y=-log(p.value),
                                          color=log(baseMean))) + geom_point(alpha=0.5) +
  ggtitle("SHAM_TACjq1 Volcano Plot") + theme_bw()
write.csv(res_tidy2_anno_protein_coding,"SHAM_TACjq1.csv",quote=F,row.names = F)



res3 <- results(dds, contrast=c("all", "TAC.veh", "TAC.JQ1")) #length(which(res3$padj<0.05)) 
res3 <- res3[order(res3$padj), ]
head(res3)
res_tidy3 = tidy.DESeqResults(res3)
head(res_tidy3)
res_tidy3_anno = res_tidy3 %>% inner_join(grcm38, by=c("gene"="symbol")) %>% arrange(p.value) 
res_tidy3_anno_protein_coding = res_tidy3_anno %>% filter(biotype=="protein_coding")

ggplot(res_tidy3_anno_protein_coding, aes(x=estimate, y=-log(p.value),
                                          color=log(baseMean))) + geom_point(alpha=0.5) +
  ggtitle("TACveh_TACjq1 Volcano Plot") + theme_bw()
write.csv(res_tidy3_anno_protein_coding,"TACveh_TACjq1.csv",quote=F,row.names = F)



write.csv(res_tidy1_anno,"SHAM_TACveh_allGenes.csv",quote=F,row.names = F)
write.csv(res_tidy2_anno,"SHAM_TACjq1_allGenes.csv",quote=F,row.names = F)
write.csv(res_tidy3_anno,"TACveh_TACjq1_allGenes.csv",quote=F,row.names = F)


#####
res1 <- results(dds, contrast=c("all","SHAM.veh", "TACveh.veh")) #length(which(res1$padj<0.05)) 
res1 <- res1[order(res1$padj), ]
head(res1)

res2 <- results(dds, contrast=c("all","SHAM.veh","TACjq1.JQ1")) #length(which(res2$padj<0.05)) 
res2 <- res2[order(res2$padj), ]
head(res2)

res3 <- results(dds, contrast=c("all", "TACveh.veh", "TACjq1.JQ1")) #length(which(res3$padj<0.05)) 
res3 <- res3[order(res3$padj), ]
head(res3)





##################
### AMR vs STA ###
##################
#1. DSEQ2
allgenes <- rownames(res1)

#results
genes_sign <- allgenes[which(res1$padj < 0.05)]
id_sign<-match(genes_sign,annotation_coding$Ensemble_id)
annotation_sign<-annotation_coding[id_sign,] 
annotation_sign$log2FC<-res1$log2FoldChange[which(res1$padj < 0.05)]
results_STA_AMR_DEseq<-annotation_sign
write.csv(results_STA_AMR_DEseq, "Results/RNAseq/results_STA_AMR_DEseq_coding.csv", row.names = F)

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_AMR_STA_Dseq2_coding.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog_coding, "AMR", "STA", results_STA_AMR_DEseq, color=COLOR[c(1,3)])
dev.off()

#2. ENET with binomial distribution
id_genes<-match(annotation_coding$Ensemble_id,rownames(norm_data_rlog))
norm_data_rlog_coding<-norm_data_rlog[na.omit(id_genes),]
results_STA_AMR_ENET<-ENET_binomial(clin, "AMR", "STA", norm_data_rlog_coding) ##50 genes
write.csv(results_STA_AMR_ENET,"Results/RNAseq/results_STA_AMR_ENET_coding.csv",row.names = F)

COLOR = brewer.pal(4,"Set2")
tiff(filename = "Results/RNAseq/heatmap_AMR_STA_ENET_coding.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog_coding, "AMR", "STA", results_STA_AMR_ENET, color=COLOR[c(1,3)])
dev.off()

id_overlap<-match(results_STA_AMR_ENET$name,results_STA_AMR_DEseq$name) #44


##################
### CMR vs STA ###
##################
#1. DSEQ2
allgenes <- rownames(res2)

#results
genes_sign <- allgenes[which(res2$padj < 0.05)] ##No genes significant

#2. ENET with binomial distribution
results_STA_CMR_ENET<-ENET_binomial(clin, "CMR", "STA", norm_data_rlog_coding) ## 1 gene
write.csv(results_STA_CMR_ENET,"Results/RNAseq/results_STA_CMR_ENET_coding.csv",row.names = F)

id<-match(results_STA_CMR_ENET$Ensemble_id,rownames(norm_data_rlog))
tiff(filename = "Results/RNAseq/heatmap_CMR_STA_ENET.tiff", width = 2000, height = 2000,  res = 300)
boxplot(norm_data_rlog[id,which(clin!="AMR")]~factor(clin[which(clin!="AMR")]),col=COLOR[c(1,2)])
dev.off()

##################
### AMR vs CMR ###
##################
#1. DSEQ2
allgenes <- rownames(res3)

#results
genes_sign <- allgenes[which(res3$padj < 0.05)]
id_sign<-match(genes_sign,annotation_coding$Ensemble_id)
annotation_sign<-annotation_coding[id_sign,] 
annotation_sign$log2FC<-res3$log2FoldChange[which(res3$padj < 0.05)]
results_CMR_AMR_DEseq<-annotation_sign
write.csv(results_CMR_AMR_DEseq, "Results/RNAseq/results_CMR_AMR_DEseq_coding.csv", row.names = F)

###Plot the results
tiff(filename = "Results/RNAseq/heatmap_CMR_AMR_Dseq2_coding.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog_coding, "AMR", "CMR", results_STA_AMR_DEseq, color=COLOR[c(2,3)])
dev.off()

#2. ENET with binomial distribution
results_CMR_AMR_ENET<-ENET_binomial(clin, "AMR", "CMR", norm_data_rlog_coding) ##1074 genes
write.csv(results_CMR_AMR_ENET,"Results/RNAseq/results_CMR_AMR_ENET_coding.csv",row.names = F)

COLOR = brewer.pal(4,"Set2")
tiff(filename = "Results/RNAseq/heatmap_AMR_CMR_ENET_coding.tiff", width = 4000, height = 3000,  res = 300)
plot_heatmap(clin,norm_data_rlog, "AMR", "CMR", results_CMR_AMR_ENET, color=COLOR[c(2,3)])
dev.off()

id_overlap<-match(results_CMR_AMR_ENET$name,results_CMR_AMR_DEseq$name) #564


#########################
## AMR vs CMR vs STA ###
########################
## 1. Using Random Forest to classify the three categories (VSURF to select the gene of interest)

#fit<-VSURF(x = t(norm_data_rlog), y=clin, parallel = TRUE,ncores=4)
#save(fit,file="Results/RNAseq/Genes_VSURF.Rdata")

load("Results/RNAseq/Genes_VSURF.Rdata")

genes<-data.frame(norm_data_rlog[fit$varselect.interp,]) ##33 genes

id<-match(rownames(genes),annotation$Ensemble_id)
rownames(genes)<-annotation$name[id]

###Plot the results
id_gene<-match(annotation$name,rownames(genes))
significantResults<-genes[na.omit(id_gene),]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin)
COLOR = brewer.pal(4,"Set2")
ann_colors = list(Type = c("STA" = COLOR [1] ,"AMR" = COLOR[3], "CMR" = COLOR[2]))

tiff(filename = "Results/RNAseq/heatmap_3categ_RF.tiff", width = 3000, height = 2000, res = 300)
pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
dev.off()


## 2. Using multinomial ENET
###Applied ENET
alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(t(norm_data_rlog_coding),clin,family="multinomial",type.multinomial = "grouped"
                                                        ,standardize=TRUE,alpha=a,nfolds=5))})
xx<-rep(NA,length(alphalist))
yy<-rep(NA,length(alphalist))
for (j in 1:length(alphalist)) {
  #print(j)
  if(class(elasticnet[[j]]) != "try-error"){
    xx[j]<-elasticnet[[j]]$lambda.min
    id.cv.opt<-grep(elasticnet[[j]]$lambda.min,elasticnet[[j]]$lambda,fixed=TRUE)
    yy[j]<-elasticnet[[j]]$cvm[id.cv.opt]
  }
}
id.min<-which(yy==min(yy,na.rm=TRUE))
lambda<-xx[id.min]
alpha<-alphalist[id.min]

enet<-glmnet(t(norm_data_rlog_coding),clin,family="multinomial",type.multinomial = "grouped",standardize=TRUE,alpha=alpha,lambda=lambda)
genes<-rownames(enet$beta[[2]])[which(enet$beta[[2]]!=0)]
coef1<-enet$beta[[1]][which(enet$beta[[1]]!=0)]
coef2<-enet$beta[[2]][which(enet$beta[[2]]!=0)]
results<-annotation_coding[match(genes,annotation_coding$Ensemble_id),]

###Obtain the FC
id_sign<-match(results$Ensemble_id,rownames(norm_data_rlog_coding))
results$mean_AMR<-apply(norm_data_rlog_coding[id_sign,which(clin=="AMR")],1,mean)
results$mean_CMR<-apply(norm_data_rlog_coding[id_sign,which(clin=="CMR")],1,mean)
results$mean_STA<-apply(norm_data_rlog_coding[id_sign,which(clin=="STA")],1,mean)

write.csv(cbind(results,coef1,coef2),"Results/RNAseq//genes.enet.multinomial_coding.csv",row.names = F)

###Plot the results
id_gene<-match(genes,rownames(norm_data_rlog_coding))
significantResults<-norm_data_rlog_coding[na.omit(id_gene),]

id<-match(rownames(significantResults),annotation_coding$Ensemble_id)
rownames(significantResults)<-annotation_coding$name[id]

xt<-t(significantResults)
xts<-scale(xt)
xtst<-t(xts)
rownames <- colnames(significantResults)
annotation.col <- data.frame(row.names = rownames)
annotation.col$Type <- factor(clin)
COLOR = brewer.pal(4,"Set2")
ann_colors = list(Type = c("STA" = COLOR [1] ,"AMR" = COLOR[3], "CMR" = COLOR[2]))

tiff(filename = "Results/RNAseq/heatmap_ENET_multinomial_coding.tiff", width = 3000, height = 2000, res = 300)
out<-pheatmap(xtst, annotation = annotation.col,border_color=F, annotation_colors = ann_colors,show_rownames=T)
dev.off()

write.csv(sort(cutree(out$tree_row, k=3)),"Results/RNAseq/geneList_clusters_ENET_multinomial_coding.csv")
plot(out$tree_row)
abline(h=9, col="red", lty=2, lwd=2)


integration<-read.csv("Results/Integration/genes.enet.alphaOptimum.csv")
intersect(integration$name,results$name) ##  "ITIH4" "TRDC" 
