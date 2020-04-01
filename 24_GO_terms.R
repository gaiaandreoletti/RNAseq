########### 24 GO terms ##############
# GO_terms <- read.table("GO_terms", head = F)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(reshape2)

setwd("/Users/gaiaandreoletti/Documents/Phyton_Score/output_7Feb2019")

################################
# GO0045823_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0045823.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0045823.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0045823.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0045823.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
TACvehiclepvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( TACvehiclepvlue_triage,TACvehiclepvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0045823_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0045823 positive regulation of heart contraction",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 12)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0045823_FPKM_TACveh + theme(legend.position="bottom")  
################################
# GO0045823_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0045823.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0045823.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
TACjq1pvlue_exp = apply(data_exp, 1, rowFisher)
##pvlue_triage = append(pvlue_triage, 1.000000e+00)
##pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( TACjq1pvlue_triage,TACjq1pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0045823_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0045823 positive regulation of heart contraction",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 12)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0045823_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0045823_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0045823.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0045823.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
Sham_pvlue_exp = apply(data_exp, 1, rowFisher)
##pvlue_triage = append(pvlue_triage, 1.000000e+00)
##pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( Sham_pvlue_triage,Sham_pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0045823_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0045823 positive regulation of heart contraction",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 12)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0045823_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0045823_positive_regulation_of_heart_contraction.#pdf")
#grid.arrange(GO0045823_FPKM_TACveh,GO0045823_FPKM_TACjq1,GO0045823_FPKM_SHAM,  ncol=1)
dev.off()


merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
 merged.data = merged.data[-1,]
                            
f_melt = melt(merged.data, id = "rank")
                            #View(f_melt)
 GO0045823_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
                            #geom_point(aes(color=variable), size=1, shape = 15) +
    xlim(0, 100) +
  labs(title = "GO0045823 positive regulation of heart contraction",
                            subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
                            x = "Rank",
                            y = "- log10  pvalue") + 
                            theme_bw() + ylim(NA, 12)   + theme (panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)
                            
                            

################################
# GO0008016_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0008016.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0008016.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0008016_FPKM_Fisher__TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0008016.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0008016_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0008016 regulation of heart contraction",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0008016_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0008016_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0008016.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0008016.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0008016_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0008016 regulation of heart contraction",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0008016_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0008016_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0008016.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0008016.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0008016_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0008016 regulation of heart contraction",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0008016_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0008016_regulation_ofHeartContraction.#pdf")
#grid.arrange(GO0008016_FPKM_TACveh,GO0008016_FPKM_TACjq1,GO0008016_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0008016_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0008016 regulation of heart contraction",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)

################################
# GO0086091_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0086091.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0086091.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0086091.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086091.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086091_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086091 regulation of heart rate by cardiac conduction",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 12)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086091_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0086091_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0086091.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086091.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086091_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086091 regulation of heart rate by cardiac conduction",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 12)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086091_FPKM_TACjq1 + theme(legend.position="bottom")  


## QUI ##
################################
# GO0086091_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0086091.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086091.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086091_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086091 regulation of heart rate by cardiac conduction",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 12)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086091_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0086091_regulationOf_heartRateByCardiacConduction.#pdf")
#grid.arrange(GO0086091_FPKM_TACveh,GO0086091_FPKM_TACjq1,GO0086091_FPKM_SHAM,  ncol=1)
dev.off()


merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086091_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086091 regulation of heart rate by cardiac conduction",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0060047_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0060047.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0060047.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0060047.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0060047.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0060047_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0060047 heart contraction",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0060047_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0060047_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0060047.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0060047.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0060047_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0060047 heart contraction",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0060047_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0060047_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0060047.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0060047.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0060047_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0060047 heart contraction",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0060047_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0060047_heartContraction.#pdf")
#grid.arrange(GO0060047_FPKM_TACveh,GO0060047_FPKM_TACjq1,GO0060047_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0060047_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0060047 heart contraction",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0002448_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0002448.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0002448.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0002448.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0002448.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0002448_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0002448 mast cell mediated immunity ",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0002448_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0002448_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0002448.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0002448.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0002448_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0002448 mast cell mediated immunity ",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0002448_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0002448_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0002448.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0002448.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0002448_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0002448 mast cell mediated immunity ",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0002448_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0002448_mastCellMediatedImmunity.#pdf")
#grid.arrange(GO0002448_FPKM_TACveh,GO0002448_FPKM_TACjq1,GO0002448_FPKM_SHAM,  ncol=1)
dev.off()


merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0002448_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0002448 mast cell mediated immunity",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0002279_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0002279.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0002279.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0002279.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0002279.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0002279_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0002279 mast cell activation involved in immune response",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0002279_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0002279_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0002279.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0002279.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0002279_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0002279 mast cell activation involved in immune response",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0002279_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0002279_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0002279.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0002279.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0002279_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0002279 mast cell activation involved in immune response",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0002279_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0002279_mast_cellActivationInvolvedInImmuneResponse.#pdf")
#grid.arrange(GO0002279_FPKM_TACveh,GO0002279_FPKM_TACjq1,GO0002279_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0002279_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0002279 mast cell activation involved in immune response",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0043303_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0043303.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0043303.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0043303.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0043303.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043303_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0043303 mast cell degranulation",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0043303_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0043303_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0043303.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0043303.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043303_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0043303 mast cell degranulation",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0043303_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0043303_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0043303.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0043303.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043303_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0043303 mast cell degranulation",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0043303_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0043303_mastCellDegranulation.#pdf")
#grid.arrange(GO0043303_FPKM_TACveh,GO0043303_FPKM_TACjq1,GO0043303_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043303_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0043303 mast cell degranulation",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)

################################
# GO0060372_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0060372.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0060372.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0060372.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0060372.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0060372_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0060372 regulation of atrial cardiac muscle cell membrane repolarization",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0060372_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0060372_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0060372.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0060372.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0060372_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0060372 regulation of atrial cardiac muscle cell membrane repolarization",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0060372_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0060372_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0060372.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0060372.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0060372_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0060372 regulation of atrial cardiac muscle cell membrane repolarization",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0060372_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0060372_regulationOfAtrial_cardiac_muscleCellMembraneRepolarization.#pdf")
#grid.arrange(GO0060372_FPKM_TACveh,GO0060372_FPKM_TACjq1,GO0060372_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0060372_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0060372 regulation of atrial cardiac muscle cell membrane repolarization",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)

################################
# GO0007159_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0007159.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0007159.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0007159.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0007159.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0007159_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0007159 leukocyte cell-cell adhesion",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 13)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0007159_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0007159_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0007159.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0007159.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0007159_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0007159 leukocyte cell-cell adhesion",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 13)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0007159_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0007159_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################
triage = read.table('FPKM/GO0007159.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0007159.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0007159_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0007159 leukocyte cell-cell adhesion",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 13)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0007159_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0007159leukocytecell-celladhesion.#pdf")
#grid.arrange(GO0007159_FPKM_TACveh,GO0007159_FPKM_TACjq1,GO0007159_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0007159_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0007159 leukocyte cell-cell adhesion",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0043542_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

# DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0043542.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
# DEgenes_TACvehvsSHAM = read.table('FPKM/GO0043542.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')

#DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0043542_FPKM_Fisher__Genes_TACveh_ovExp_TACjq1.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0043542_FPKM_Fisher__DEgenes_Genes_TACveh_ovExp_SHAM.txt')
JQ1_attenuated <- read.table('FPKM/GO0043542_FPKM_Fisher__Genes_Gene_JQ1_attenuated.txt', head =F)
  

data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

data_JQ1_attenuated <- 
  data.frame(P_A=c(JQ1_attenuated$V1),
             NP_A=c(JQ1_attenuated$V2),
             P_B=c(JQ1_attenuated$V3),
             NP_B=c(JQ1_attenuated$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_data_JQ1_attenuated = apply(data_JQ1_attenuated, 1, rowFisher)

triage = read.table('FPKM/GO0043542.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0043542.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
TACvehiclepvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_data_JQ1_attenuated)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","Stess_induced","JQ1_attenuated")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043542_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0043542 endothelial cell migration",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 18)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)
GO0043542_FPKM_TACveh + theme(legend.position="bottom")  

################################
# GO0043542_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0043542.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0043542.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
TACjq1pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043542_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0043542 endothelial cell migration",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 18)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0043542_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0043542_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0043542.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0043542.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Shampvlue_triage = apply(data_triage, 1, rowFisher)
Shampvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043542_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0043542 endothelial cell migration",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 18)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0043542_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0043542_endothelial_cellMigration_22Apr2019.#pdf")
#grid.arrange(GO0043542_FPKM_TACveh,GO0043542_FPKM_TACjq1,GO0043542_FPKM_SHAM,  ncol=1)
dev.off()


###
merged.data = cbind( Shampvlue_triage, TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_Triage", "TACveh_Triage","TACjq1_Triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043542_FPKM_all <- ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(1, 100) + ylim(NA, 18) +
  labs(title = "GO0043542 endothelial cell migration",
       subtitle = "SHAM, TAC veh, and TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw()  + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


######
merged.data = cbind( Shampvlue_triage, TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_Triage", "TACveh_Triage","TACjq1_Triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0043542_FPKM_all <- ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(1, 100) + ylim(NA, 18) +
  labs(title = "GO0043542 endothelial cell migration",
       subtitle = "SHAM, TAC veh, and TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw()  + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0032757_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0032757.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0032757.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0032757.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0032757.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0032757_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0032757 positive regulation of interleukin-8 production ",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 3)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0032757_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0032757_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0032757.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0032757.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0032757_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0032757 positive regulation of interleukin-8 production ",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 3)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0032757_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0032757_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0032757.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0032757.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0032757_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0032757 positive regulation of interleukin-8 production ",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 3)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0032757_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0032757_positiveRegulation_of_interleukin-8_production.#pdf")
#grid.arrange(GO0032757_FPKM_TACveh,GO0032757_FPKM_TACjq1,GO0032757_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0032757_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0032757 positive regulation of interleukin-8 production",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0014910_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

#DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0014910.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
#DEgenes_TACvehvsSHAM = read.table('FPKM/GO0014910.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0014910_FPKM_Fisher__Genes_TACveh_ovExp_TACjq1.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0014910_FPKM_Fisher__DEgenes_Genes_TACveh_ovExp_SHAM.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0014910.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0014910.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
TACvehiclepvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0014910_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0014910 regulation of smooth muscle cell migration",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 16)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0014910_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0014910_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0014910.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0014910.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
TACjq1pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0014910_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0014910 regulation of smooth muscle cell migration",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 16)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0014910_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0014910_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0014910.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0014910.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Shampvlue_triage = apply(data_triage, 1, rowFisher)
Shampvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0014910_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0014910 regulation of smooth muscle cell migration",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 16)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0014910_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0014910_regulationOfsmoothMuscleCellMigration_23Apr2019.#pdf")
#grid.arrange(GO0014910_FPKM_TACveh,GO0014910_FPKM_TACjq1,GO0014910_FPKM_SHAM,  ncol=1)
dev.off()

####
###
merged.data = cbind( Shampvlue_triage,Shampvlue_exp, TACvehiclepvlue_triage,TACvehiclepvlue_exp,TACjq1pvlue_triage,TACjq1pvlue_exp)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_Triage","SHAM_Expression", "TACveh_Triage","TACveh_Expression","TACjq1_Triage","TACvehjq1_Expression")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(1, 100) + ylim(NA, 18) +
  labs(title = "GO0014910 regulation of smooth muscle cell migration",
       subtitle = "SHAM, TAC veh, and TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw()  + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


######
merged.data = cbind( Shampvlue_triage, TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_Triage", "TACveh_Triage","TACjq1_Triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0014910_FPKM_all <- ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(1, 100) + ylim(NA, 18) +
  labs(title = "GO0014910 regulation of smooth muscle cell migration",
       subtitle = "SHAM, TAC veh, and TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw()  + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0099623_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0099623.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0099623.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0099623.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0099623.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0099623_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0099623 regulation of cardiac muscle cell membrane repolarization ",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0099623_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0099623_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0099623.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0099623.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0099623_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0099623 regulation of cardiac muscle cell membrane repolarization ",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0099623_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0099623_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0099623.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0099623.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0099623_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0099623 regulation of cardiac muscle cell membrane repolarization ",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0099623_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0099623_regulationOfCardiacMuscleCellMembraneRepolarization.#pdf")
#grid.arrange(GO0099623_FPKM_TACveh,GO0099623_FPKM_TACjq1,GO0099623_FPKM_SHAM,  ncol=1)
dev.off()


merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0099623_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0099623 regulation of cardiac muscle cell membrane repolarization",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)



################################
# GO0032677_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0032677.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0032677.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0032677.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0032677.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0032677_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0032677 regulation of interleukin-8 production",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0032677_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0032677_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0032677.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0032677.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0032677_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0032677 regulation of interleukin-8 production",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0032677_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0032677_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0032677.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0032677.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0032677_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0032677 regulation of interleukin-8 production",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0032677_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0032677_regulationOfInterleukin-8production.#pdf")
#grid.arrange(GO0032677_FPKM_TACveh,GO0032677_FPKM_TACjq1,GO0032677_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0032677_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0032677 regulation of interleukin-8 production",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)



################################
# GO0086012_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0086012.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0086012.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0086012.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086012.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086012_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086012 membrane depolarization during cardiac muscle cell action potential ",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086012_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0086012_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0086012.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086012.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086012_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086012 membrane depolarization during cardiac muscle cell action potential ",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086012_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0086012_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0086012.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086012.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086012_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086012 membrane depolarization during cardiac muscle cell action potential ",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086012_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0086012_membraneDepolarizationDuringCardiacMuscleCellActionPotential.#pdf")
#grid.arrange(GO0086012_FPKM_TACveh,GO0086012_FPKM_TACjq1,GO0086012_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086012_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086012 membrane depolarization during cardiac muscle cell action potential",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)



################################
# GO0086003_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0086003.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0086003.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0086003.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086003.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086003_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086003 cardiac muscle cell contraction",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086003_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0086003_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0086003.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086003.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086003_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086003 cardiac muscle cell contraction",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086003_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0086003_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0086003.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0086003.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086003_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086003 cardiac muscle cell contraction",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0086003_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0086003_cardiacMuscleCellContraction.#pdf")
#grid.arrange(GO0086003_FPKM_TACveh,GO0086003_FPKM_TACjq1,GO0086003_FPKM_SHAM,  ncol=1)
dev.off()


merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0086003_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0086003 cardiac muscle cell contraction",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)

################################
# GO0010818_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0010818.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0010818.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0010818.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0010818.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0010818_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0010818 T cell chemotaxis ",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 3)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0010818_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0010818_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0010818.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0010818.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0010818_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0010818 T cell chemotaxis ",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 3)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0010818_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0010818_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0010818.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0010818.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0010818_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0010818 T cell chemotaxis ",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 3)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0010818_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0010818_Tcellchemotaxis.#pdf")
#grid.arrange(GO0010818_FPKM_TACveh,GO0010818_FPKM_TACjq1,GO0010818_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0010818_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0010818 T cell chemotaxis",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)



################################
# GO0072678_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0072678.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0072678.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0072678.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0072678.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0072678_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0072678 T cell migration ",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0072678_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0072678_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0072678.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0072678.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0072678_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0072678 T cell migration ",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0072678_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0072678_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0072678.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0072678.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0072678_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0072678 T cell migration ",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 4)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0072678_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0072678_TcellMigration.#pdf")
#grid.arrange(GO0072678_FPKM_TACveh,GO0072678_FPKM_TACjq1,GO0072678_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0072678_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0072678 T cell migration",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO2000482_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO2000482.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO2000482.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO2000482.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO2000482.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO2000482_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO2000482 regulation of interleukin-8 secretion",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO2000482_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO2000482_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO2000482.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO2000482.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO2000482_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO2000482 regulation of interleukin-8 secretion",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO2000482_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO2000482_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO2000482.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO2000482.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO2000482_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO2000482 regulation of interleukin-8 secretion",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO2000482_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO2000482_regulation_of_interleukin-8secretion.#pdf")
#grid.arrange(GO2000482_FPKM_TACveh,GO2000482_FPKM_TACjq1,GO2000482_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO2000482_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO2000482 regulation of interleukin-8 secretion",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)

################################
# GO0007015_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0007015.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0007015.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0007015.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0007015.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0007015_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0007015 actin filament organization",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 10)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0007015_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0007015_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0007015.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0007015.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0007015_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0007015 actin filament organization",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 10)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0007015_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0007015_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0007015.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0007015.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0007015_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0007015 actin filament organization",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 10)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0007015_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0007015_actin_filament_organization.#pdf")
#grid.arrange(GO0007015_FPKM_TACveh,GO0007015_FPKM_TACjq1,GO0007015_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0007015_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0007015 actin filament organization",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)



################################
# GO0050672_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0050672.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0050672.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0050672.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0050672.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0050672_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0050672 negative regulation of lymphocyte proliferation",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0050672_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0050672_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0050672.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0050672.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0050672_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0050672 negative regulation of lymphocyte proliferation",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)  +   theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0050672_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0050672_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0050672.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0050672.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0050672_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0050672 negative regulation of lymphocyte proliferation",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0050672_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0050672_negative_regulation_of_lymphocyte_proliferation.#pdf")
#grid.arrange(GO0050672_FPKM_TACveh,GO0050672_FPKM_TACjq1,GO0050672_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0050672_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0050672 negative regulation of lymphocyte proliferation",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO0033622_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO0033622.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO0033622.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO0033622.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0033622.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0033622_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0033622 integrin activation",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0033622_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO0033622_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0033622.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0033622.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0033622_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0033622 integrin activation",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0033622_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO0033622_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO0033622.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO0033622.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0033622_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0033622 integrin activation",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 5)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0033622_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO0033622_integrin_activation.#pdf")
#grid.arrange(GO0033622_FPKM_TACveh,GO0033622_FPKM_TACjq1,GO0033622_FPKM_SHAM,  ncol=1)
dev.off()


merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0033622_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO0033622 integrin activation",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


################################
# GO2000147_FPKM_Fisher__TACvehicle_discordance_score_genes #######
#############################

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

DEgenes_TACvehvsTACjq1 = read.table('FPKM/GO2000147.txt_FPKM_Fisher__DEgenes_TACvehvsTACjq1.txt.txt')
DEgenes_TACvehvsSHAM = read.table('FPKM/GO2000147.txt_FPKM_Fisher__DEgenes_TACvehvsSHAM.txt.txt')


data_DEgenes_TACvehvsSHAM <- 
  data.frame(P_A=c(DEgenes_TACvehvsSHAM$V1),
             NP_A=c(DEgenes_TACvehvsSHAM$V2),
             P_B=c(DEgenes_TACvehvsSHAM$V3),
             NP_B=c(DEgenes_TACvehvsSHAM$V4))

DEgenes_TACvehvsTACjq1 <- 
  data.frame(P_A=c(DEgenes_TACvehvsTACjq1$V1),
             NP_A=c(DEgenes_TACvehvsTACjq1$V2),
             P_B=c(DEgenes_TACvehvsTACjq1$V3),
             NP_B=c(DEgenes_TACvehvsTACjq1$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

pvlue_data_DEgenes_TACvehvsSHAM = apply(data_DEgenes_TACvehvsSHAM, 1, rowFisher)
pvlue_DEgenes_TACvehvsTACjq1 = apply(DEgenes_TACvehvsTACjq1, 1, rowFisher)

triage = read.table('FPKM/GO2000147.txt_FPKM_Fisher_TACvehicle_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO2000147.txt_FPKM_Fisher_TAC_vehicle_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehiclepvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
# #pvlue_triage = append(pvlue_triage, 1.000000e+00)
# #pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp, pvlue_data_DEgenes_TACvehvsSHAM,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","DEgenes_TACvehvsSHAM","DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO2000147_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO2000147 positive regulation of cell motility",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO2000147_FPKM_TACveh + theme(legend.position="bottom")  



################################
# GO2000147_FPKM_Fisher__TACjq1icle_discordance_score_genes
#############################

triage = read.table('FPKM/GO2000147.txt_FPKM_Fisher_TACj_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO2000147.txt_FPKM_Fisher_TACjq1_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_DEgenes_TACvehvsTACjq1)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsTACjq1")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO2000147_FPKM_TACjq1 <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO2000147 positive regulation of cell motility",
       subtitle = "TAC jq1 data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO2000147_FPKM_TACjq1 + theme(legend.position="bottom")  

################################
# GO2000147_FPKM_Fisher__SHAMicle_discordance_score_genes
#############################

triage = read.table('FPKM/GO2000147.txt_FPKM_Fisher_Sham_discordance_score_genes.txt', header = F)
exp = read.table('FPKM/GO2000147.txt_FPKM_Fisher_Sham_expression_value_genes.txt', header = F)

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

Sham_pvlue_triage = apply(data_triage, 1, rowFisher)
pvlue_exp = apply(data_exp, 1, rowFisher)
#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)
merged.data = cbind( pvlue_triage,pvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "DEgenes_TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO2000147_FPKM_SHAM <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO2000147 positive regulation of cell motility",
       subtitle = "SHAM data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() + ylim(NA, 40)   + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO2000147_FPKM_SHAM + theme(legend.position="bottom")  

#pdf("GO2000147_positiveRegulationOfcellMotility.#pdf")
#grid.arrange(GO2000147_FPKM_TACveh,GO2000147_FPKM_TACjq1,GO2000147_FPKM_SHAM,  ncol=1)
dev.off()

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO2000147_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO2000147 positive regulation of cell motility",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)




################################
# GO0005578_extr_matrix.txt #####
#############################

GO0005578_FPKM_Sham_fisher_file_triage = read.table('FPKM/GO0005578_FPKM_Fisher__Sham_discordance_score_genes.txt', header = F)
GO0005578_FPKM_Sham_fisher_file_exp = read.table('FPKM/GO0005578_FPKM_Fisher__Sham_expression_value_genes.txt', header = F)

GO0005578_FPKM_TACveh_fisher_file_triage = read.table('FPKM/GO0005578_FPKM_Fisher__TACvehicle_discordance_score_genes.txt', header = F)
GO0005578_FPKM_TACveh_fisher_file_exp = read.table('FPKM/GO0005578_FPKM_Fisher__TAC_vehicle_expression_value_genes.txt', header = F)

GO0005578_FPKM_TACjq1_fisher_file_triage = read.table('FPKM/GO0005578_FPKM_Fisher__TACj_discordance_score_genes.txt', header = F)
GO0005578_FPKM_TACjq1_fisher_file_exp = read.table('FPKM/GO0005578_FPKM_Fisher__TACjq1_expression_value_genes.txt', header = F)

## Diff express ex matrix ####
## TAC veh / SHAM & ## TAC veh / TAC jq1 ### ###

# SHAM
triage = GO0005578_FPKM_Sham_fisher_file_triage
exp = GO0005578_FPKM_Sham_fisher_file_exp

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

SHAMpvlue_triage = apply(data_triage, 1, rowFisher)
SHAMpvlue_exp = apply(data_exp, 1, rowFisher)

#pvlue_triage = append(pvlue_triage, 1.000000e+00)
#pvlue_exp = append(pvlue_exp, 1.000000e+00)

merged.data = cbind( SHAMpvlue_triage,SHAMpvlue_exp,pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression", "TACvehvsSHAM")
dim(merged.data)
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0005578_FPKM_Sham <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(1, 100) + ylim(NA, 40) +
  labs(title = "GO 0005578 - extracellular matrix",
       subtitle = "FPKM Sham data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw()  + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


GO0005578_FPKM_Sham_data = merged.data
GO0005578_FPKM_Sham + theme(legend.position="bottom")  



triage = GO0005578_FPKM_TACveh_fisher_file_triage 
exp = GO0005578_FPKM_TACveh_fisher_file_exp 

data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

rowFisher <- function(x, ...) {
  return(fisher.test(matrix(x, nrow = 2, ...), alternative = "greater")$p.value)
}

TACvehpvlue_triage = apply(data_triage, 1, rowFisher)
TACvehpvlue_exp = apply(data_exp, 1, rowFisher)

merged.data = cbind(TACvehpvlue_triage,TACvehpvlue_exp,pvlue_JQ1_attenuated, pvlue_data_DEgenes_TACvehvsSHAM)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("Triage","Expression","JQ1_attenuated", "TACvehvsSHAM")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0005578_FPKM_TACveh <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) + ylim(NA, 40) +
  labs(title = "GO 0005578 - extracellular matrix",
       subtitle = "TAC veh data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw()  + theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)


triage = GO0005578_FPKM_TACjq1_fisher_file_triage 
exp = GO0005578_FPKM_TACjq1_fisher_file_exp 
data_triage <- 
  data.frame(P_A=c(triage$V1),
             NP_A=c(triage$V2),
             P_B=c(triage$V3),
             NP_B=c(triage$V4))

data_exp <- 
  data.frame(P_A=c(exp$V1),
             NP_A=c(exp$V2),
             P_B=c(exp$V3),
             NP_B=c(exp$V4))

TACjq1pvlue_triage = apply(data_triage, 1, rowFisher)

merged.data = cbind(Sham_pvlue_triage,TACvehiclepvlue_triage,TACjq1pvlue_triage)
merged.data = as.data.frame(merged.data)
colnames(merged.data) <- c("SHAM_triage","TAC_triage", "TACjq1_triage")
merged.data$rank =  seq(1:nrow(merged.data))
dim(merged.data)
head(merged.data,10)
merged.data = merged.data[-1,]

f_melt = melt(merged.data, id = "rank")
#View(f_melt)
GO0005578_FPKM_all <-  ggplot(f_melt, aes(x = rank, y = -log10(value), colour = variable)) +
  geom_line(aes(color=variable)) +
  #geom_point(aes(color=variable), size=1, shape = 15) +
  xlim(0, 100) +
  labs(title = "GO 0005578 - extracellular matrix",
       subtitle = "SHAM, TAC veh, TAC JQ1 traige data",
       x = "Rank",
       y = "- log10  pvalue") + 
  theme_bw() +  theme (panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +theme(aspect.ratio = 1)





pdf("TRIAGE_only_GO_terms.pdf")
grid.arrange(GO0045823_FPKM_all,
             GO0008016_FPKM_all,
             GO0086091_FPKM_all, ncol=1)
grid.arrange( GO0060047_FPKM_all,
              GO0002448_FPKM_all,
              GO0002279_FPKM_all, ncol=1)
grid.arrange(GO0043303_FPKM_all,
             GO0060372_FPKM_all,
             GO0007159_FPKM_all, ncol=1)
grid.arrange(GO0043542_FPKM_all,
             GO0014910_FPKM_all, 
             GO0086012_FPKM_all, ncol=1)
grid.arrange(GO0010818_FPKM_all,
             GO0072678_FPKM_all,
             GO2000482_FPKM_all, ncol=1)
grid.arrange( GO0007015_FPKM_all,
              GO0050672_FPKM_all, 
              GO0033622_FPKM_all,ncol=1)
grid.arrange(GO2000147_FPKM_all, GO0005578_FPKM_all, e, ncol=1)
dev.off()
