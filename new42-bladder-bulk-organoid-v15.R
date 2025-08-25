####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########      ########      ########     ########             #########           ##       ############      ##
##     ##     ##     ##     ##     ##           ##     ##            ##      ##         ####           ##          ####
##     ##     ##     ##    ##      ##           ##     ##            ##       ##       ##  ##          ##         ##  ##
##     ########      #######       ########     ########             ##        ##     ##    ##         ##        ##    ##
##     ##            ##    ##      ##           ##                   ##       ##     ##########        ##       ##########
##     ##            ##     ##     ##           ##                   ##      ##     ##        ##       ##      ##        ##
##     ##            ##      ##    ########     ##                   #########     ##          ##      ##     ##          ##
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#DATE:   1-23-25
#GOAL:   Basic Parsing of GEO Dataset 
#PLAN:
setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 42 bladder cancer")

{

load("Bladder40_Organoid_paired.RData")

#####################################################
#1.1 HISTOGRAMS:  make sure data is log-normally distributed
#####################################################

hist(log10(BladderOrganoid_counts+1),breaks=50,ylim=c(0,40000),main="Gene Expression is ~ Log-Normal Distributed")

#####################################################
#1.2 BARPLOTS:  make sure each sample has good "read depth"
#####################################################


read_depth<-colSums(BladderOrganoid_counts)

barplot(sort(read_depth),las=2,cex.names = 0.6,main="Read Depth = Amount mRNA/sample",)
abline(h=10e6,lty=2)


#####################################################
#1.3 NORMALIZE COLUMNS (read depth) 
#####################################################


barplot(sort(colSums(BladderOrganoid_tpm)),las=2,cex.names = 0.6,main="Gene Expression normalized to Depth")



#####################################################
#2.1 NORMALIZE ROWS (z-score)
#####################################################

#FILTER OUT ALL ZERO ROWS:
zero_rows<-which(rowSums(BladderOrganoid_tpm)==0)
organoid_logRPM_zero<-BladderOrganoid_tpm[-zero_rows,]

BladderOrganoid_tpm<-organoid_logRPM_zero

#CONVERT TO Z-SCORE
allSamples_zscore<-t(scale(t(BladderOrganoid_tpm),center=T,scale=T))


#####################################################
#NORMALIZE TISSUE AND ORGANOID SEPARALY
#####################################################
organoid_zscore<-t(scale(t(BladderOrganoid_tpm[,grep("organoid",colnames(BladderOrganoid_tpm))]),center=T,scale=T))
tissue_zscore<-t(scale(t(BladderOrganoid_tpm[,grep("issue",colnames(BladderOrganoid_tpm))]),center=T,scale=T))

separated_organ_tissue_zscore_raw<-cbind(organoid_zscore,tissue_zscore)

separated_organ_tissue_zscore<-separated_organ_tissue_zscore_raw[-which(is.na(rowSums(separated_organ_tissue_zscore_raw))==T),]

#####################################################
#2.2 SAVE as CSV (to open in excel) and as RData file (so easy to do future R analysis on)
#####################################################


write.csv(allSamples_zscore,file="allSamples_zscore.csv")

save(list=c("allSamples_zscore","BladderOrganoid_tpm","separated_organ_tissue_zscore"),file="Bladder42_zscores.RData")

}  
  


####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ##      ###     ##       ##        ##        ##      ##   ########    ##     #######                               
##      ####     ####    ##      ####       ##         ##    ##   ##           ##    ##                                               
##     ##  ##    ## ##   ##     ##  ##      ##          ##  ##    ##           ##    ##           
##    ##    ##   ##  ##  ##    ##    ##     ##           ####      #######     ##     #######                                           
##    ########   ##   ## ##   ##########    ##            ##             ##    ##           ##                                      
##    ##    ##   ##    ####   ##      ##    ##            ##             ##    ##           ##                                          
##    ##    ##   ##     ###   ##      ##    ########      ##       #######     ##     #######                                                 
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#DATE:   1-23-25
#GOAL:   Basic Analysis of RNAseq dataset
#PLAN:

#       1.  UNSUPERVISED ANALYSIS:   
#                           -> CORRELATION MATRIX:   are replicates correlated?   what groups are there?
#                           -> PRINCIPAL COMPONENT ANALYSIS:   what are the 2 major "dimensions" of varaibility in samples
#                           -> PATHWAY ANALYSIS:  transcriptional programs
#
#################################################################################################################################################################################################################################################################### 

setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 42 bladder cancer")

{
#load("BladderOrganoid_paired.RData")
load("Bladder42_zscores.RData")

library("gplots")
library("RColorBrewer")

###############################################################################################################################################################
#PART 1:   ORGANOIDS VS TISSUE
###############################################################################################################################################################


#####################################################
#FIGURE 2A:  CORR-MATRIX.   TPM
#####################################################

organoid_rpm_corr<-cor(BladderOrganoid_tpm)

heatmap.2(organoid_rpm_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "rpm corr",key.xlab = "rpm corr",cexRow = 0.4,cexCol = 0.6)

#####################################################
#FIGURE 2A:   CORR-MATRIX.  z-score
#####################################################

organoid_zscore_corr<-cor(allSamples_zscore)

heatmap.2(organoid_zscore_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "z-score corr",key.xlab = "z-score corr",cexRow = 0.6,cexCol = 0.6)





#####################################################
#FIGURE 2B:   PCA   color by patients (below) or conditiions
#####################################################

table(sapply(strsplit(colnames(allSamples_zscore),"_"), function(x) x[2]))
#    Bo Murphy  Ohren  Pearl  Penny Sadie2 Sadie3 Sadie4 Sadie7    Tux 
#.    2      3      3      3     11      2      6      3      3      4 

#RUN PCA
organoid_pca<-prcomp(t(allSamples_zscore))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)


sample_names_vector[grep("tissue_Penny",sample_names_vector)]<-"darkblue"
sample_names_vector[grep("organoid_Penny",sample_names_vector)]<-"blue"
sample_names_vector[grep("Sadie2",sample_names_vector)]<-"red1"
sample_names_vector[grep("Sadie3",sample_names_vector)]<-"red2"
sample_names_vector[grep("Sadie4",sample_names_vector)]<-"red3"
sample_names_vector[grep("Sadie7",sample_names_vector)]<-"red4"
sample_names_vector[grep("Bo",sample_names_vector)]<-"darkorange"
sample_names_vector[grep("Murphy",sample_names_vector)]<-"darkgreen"
sample_names_vector[grep("Ohren",sample_names_vector)]<-"purple"
sample_names_vector[grep("Pearl",sample_names_vector)]<-"black"
sample_names_vector[grep("Tux",sample_names_vector)]<-"grey30"



plot(organoid_pca_data[,1],organoid_pca_data[,3],col=sample_names_vector,pch=16,
     xlab="PC1 (38%):  Penny tissue vs organoid",ylab="PC3 (12%): tux-tissue v organoid")
text(organoid_pca_data[,1],organoid_pca_data[,3],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)



plot(organoid_pca_data[,2],organoid_pca_data[,4],col=sample_names_vector,pch=16,
     xlab="PC2 (14%):  organoid vs tissue",ylab="PC4 (9%):")
text(organoid_pca_data[,2],organoid_pca_data[,4],rownames(organoid_pca_data),col=sample_names_vector,cex=0.5)


#####################################################
#heatmaps
#####################################################

#RANK PATHWAYS BASED ON VARIANCE:
genes_in_order<-names(sort(apply(BladderOrganoid_tpm,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(allSamples_zscore[genes_in_order[1:500],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.1,cexCol = 0.5,mar=c(5,10))





#####################################################
#FIGURE 3:   PATHWAY ANALYSIS
#####################################################

#####################
#LOAD PACKAGES
#####################


library(viper)
library(msigdbr)



#####################
#2 MAKE HALLMARK-OBJECT FOR PATHWAY ANALYSIS
#####################

#EXTRACT MOUSE-HALLMARKS:
HALLMARK_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
HALLMARK_gene_sets$gs_name<-gsub("HALLMARK_","",HALLMARK_gene_sets$gs_name)


#CONVERT MOUSE-HALLMARK DATA-FRAME TO A "REGULON"
msigdbr_list = split(x = HALLMARK_gene_sets$gene_symbol, f = HALLMARK_gene_sets$gs_name)
hallmark_regulon<-list()

for (tf in names(msigdbr_list)){
  sig_genes<-msigdbr_list[[tf]]
  
  hallmark_regulon[[tf]]<-list()
  tmp_tfmode<-rep(1,length(sig_genes))
  names(tmp_tfmode)<-sig_genes
  
  hallmark_regulon[[tf]]$'tfmode'<-tmp_tfmode
  hallmark_regulon[[tf]]$'likelihood'<-as.numeric(tmp_tfmode)
}

#####################
#RUN VIPER-Pathway analysis
#####################
organoid_Hallmarks <- viper(allSamples_zscore,hallmark_regulon,method="scale")



#####################
#VISUALIZE PATHWAY ANALYSIS:
#####################

#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(organoid_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:30],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.4,mar=c(5,10))

}


####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########      ########      ########     ########             #########           ##       ############      ##
##     ##     ##     ##     ##     ##           ##     ##            ##      ##         ####           ##          ####
##     ##     ##     ##    ##      ##           ##     ##            ##       ##       ##  ##          ##         ##  ##
##     ########      #######       ########     ########             ##        ##     ##    ##         ##        ##    ##
##     ##            ##    ##      ##           ##                   ##       ##     ##########        ##       ##########
##     ##            ##     ##     ##           ##                   ##      ##     ##        ##       ##      ##        ##
##     ##            ##      ##    ########     ##                   #########     ##          ##      ##     ##          ##
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#DATE:   1-23-25
#GOAL:   OLD COMBINE BATCHES OF DATA:    ROS_CFam_1
#PLAN:
setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 42 bladder cancer")


load("Combined_Bladder_batches.RData")


library("gplots")
library("RColorBrewer")



####################################################################################################################################################################################################################################################################
#0 PREPARE COMBINED BATCHES:  loaded above
#################################################################################################################################################################################################################################################################### 

{
load("Bladder40_Organoid_paired.RData")
samples40<-BladderOrganoid_tpm


load("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/04 10 paired organoid-tissue/BladderOrganoid_paired.RData")
samples14<-BladderOrganoid_tpm
new_colnames<-paste(paste(sapply(strsplit(colnames(BladderOrganoid_tpm),split="_"),function(x) x[2]),
                    sapply(strsplit(colnames(BladderOrganoid_tpm),split="_"),function(x) x[1]),sep="_"),
                    sapply(strsplit(colnames(BladderOrganoid_tpm),split="_"),function(x) x[3]),sep="_")

colnames(samples14)<-new_colnames

overlap_genes<-intersect(rownames(samples14),rownames(samples40))

Overall_Bladder_tpm<-cbind(samples40[overlap_genes,],samples14[overlap_genes,])

#FILTER OUT ALL ZERO ROWS:
zero_rows<-which(rowSums(Overall_Bladder_tpm)==0)
organoid_logRPM_zero<-Overall_Bladder_tpm[-zero_rows,]

Overall_Bladder_tpm<-organoid_logRPM_zero


Overall_Bladder_zscore<-t(scale(t(Overall_Bladder_tpm)))


#####################################################
#NORMALIZE TISSUE AND ORGANOID SEPARALY
#####################################################
organoid_zscore<-t(scale(t(Overall_Bladder_tpm[,grep("organoid",colnames(Overall_Bladder_tpm))]),center=T,scale=T))
tissue_zscore<-t(scale(t(Overall_Bladder_tpm[,grep("issue",colnames(Overall_Bladder_tpm))]),center=T,scale=T))

Overall_organ_tissue_zscore_raw<-cbind(organoid_zscore,tissue_zscore)

Overall_organ_tissue_zscore<-Overall_organ_tissue_zscore_raw[-which(is.na(rowSums(Overall_organ_tissue_zscore_raw))==T),]


write.csv(log10(Overall_Bladder_tpm+1),file="Overall_Bladder_logTPM.csv")

save(list=c("Overall_Bladder_tpm","Overall_Bladder_zscore","Overall_organ_tissue_zscore"),file="Combined_Bladder_batches.RData")




###############################################################################################################################################################
#PART 1:   ORGANOIDS VS TISSUE
###############################################################################################################################################################


#####################################################
#FIGURE 2A:  CORR-MATRIX.    Published analysis
#####################################################

organoid_rpm_corr<-cor(Overall_Bladder_tpm)

heatmap.2(organoid_rpm_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "rpm corr",key.xlab = "rpm corr",cexRow = 0.4,cexCol = 0.6)

#####################################################
#FIGURE 2A:   CORR-MATRIX.   Refined analysis
#####################################################

organoid_zscore_corr<-cor(Overall_Bladder_zscore)

heatmap.2(organoid_zscore_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "z-score corr",key.xlab = "z-score corr",cexRow = 0.3,cexCol = 0.6)


organoid_zscore_corr<-cor(Overall_organ_tissue_zscore)

heatmap.2(organoid_zscore_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "z-score corr",key.xlab = "z-score corr",cexRow = 0.3,cexCol = 0.6)





#####################################################
#PCA:  all
#####################################################

#RUN PCA
organoid_pca<-prcomp(t(Overall_Bladder_zscore))
organoid_pca_data<-organoid_pca$x
PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)

plot(organoid_pca_data[,1],organoid_pca_data[,2],col="black",pch=16,
     xlab="PC1 (29%):  tissue vs organoid",ylab="PC2 (15%): tissue v organoid", main="PCA:  all samples")
text(organoid_pca_data[,1],organoid_pca_data[,2],rownames(organoid_pca_data),col="black",cex=0.3)



#####################################################
#PCA:  organoids only
#####################################################

#RUN PCA
organoid_pca<-prcomp(t(Overall_organ_tissue_zscore[,grep("ganoid",colnames(Overall_organ_tissue_zscore))]))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)
table(sapply(strsplit(sample_names_vector,"_"), function(x) x[2]))
#Bo    boo    liz  molly Murphy  Ohren  Pearl  Penny  pilot  sadie Sadie2 Sadie3 Sadie4 Sadie7    Tux 
#2      3      1      1      3      3      3      5      1      2      2      6      3      3      3 

name2color<-c("darkblue","blue","red","darkred","purple","darkorange","yellow4","cyan4","pink3","darkgreen","black")
names(name2color)<-c("Bo","boo","liz","molly","Murphy","Ohren","Pearl","Penny","pilot","adie","Tux")

for (idx in 1:length(sample_names_vector)){
  sample_names_vector[grep(names(name2color)[idx],sample_names_vector)]<-name2color[idx]
}

tmp_names<-gsub("_rep","",gsub("organoid_","",rownames(organoid_pca_data)))


plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16,
     xlab="PC1 (24%):  tissue vs organoid",ylab="PC2 (17%): Tux vs all",main="PCA:  Organoids only")
text(organoid_pca_data[,1],organoid_pca_data[,2],tmp_names,col=sample_names_vector,cex=0.3)

plot(organoid_pca_data[,1],organoid_pca_data[,3],col=sample_names_vector,pch=16,
     xlab="PC1 (24%):  tissue vs organoid",ylab="PC2 (11%):",main="PCA:  Organoids only")
text(organoid_pca_data[,1],organoid_pca_data[,3],tmp_names,col=sample_names_vector,cex=0.3)


#####################################################
#PCA:  organoids only
#####################################################

#RUN PCA
organoid_pca<-prcomp(t(Overall_organ_tissue_zscore[,grep("issue",colnames(Overall_organ_tissue_zscore))]))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)
table(sapply(strsplit(sample_names_vector,"_"), function(x) x[2]))
#boo  Penny  pilot  sadie sadie3    Tux 
#2      6      1      1      2      1 

name2color<-c("blue","cyan4","pink3","darkgreen","black")
names(name2color)<-c("boo","Penny","pilot","adie","Tux")

for (idx in 1:length(sample_names_vector)){
  sample_names_vector[grep(names(name2color)[idx],sample_names_vector)]<-name2color[idx]
}

tmp_names<-gsub("_rep","",gsub("tissue_","",rownames(organoid_pca_data)))


plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16,
     xlab="PC1 (47%):  ",ylab="PC2 (12%): ",main="PCA:  Tissue Only")
text(organoid_pca_data[,1],organoid_pca_data[,2],tmp_names,col=sample_names_vector,cex=0.3)




#####################################################
#heatmaps
#####################################################

#RANK PATHWAYS BASED ON VARIANCE:
genes_in_order<-names(sort(apply(Overall_Bladder_tpm,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_zscore[genes_in_order[1:2000],]),col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.3,cexCol = 0.1,mar=c(5,10))





#####################################################
#FIGURE 3:   PATHWAY ANALYSIS
#####################################################

#####################
#LOAD PACKAGES
#####################


library(viper)
library(msigdbr)



#####################
#2 MAKE HALLMARK-OBJECT FOR PATHWAY ANALYSIS
#####################

#EXTRACT MOUSE-HALLMARKS:
HALLMARK_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
HALLMARK_gene_sets$gs_name<-gsub("HALLMARK_","",HALLMARK_gene_sets$gs_name)


#CONVERT MOUSE-HALLMARK DATA-FRAME TO A "REGULON"
msigdbr_list = split(x = HALLMARK_gene_sets$gene_symbol, f = HALLMARK_gene_sets$gs_name)
hallmark_regulon<-list()

for (tf in names(msigdbr_list)){
  sig_genes<-msigdbr_list[[tf]]
  
  hallmark_regulon[[tf]]<-list()
  tmp_tfmode<-rep(1,length(sig_genes))
  names(tmp_tfmode)<-sig_genes
  
  hallmark_regulon[[tf]]$'tfmode'<-tmp_tfmode
  hallmark_regulon[[tf]]$'likelihood'<-as.numeric(tmp_tfmode)
}

#####################
#RUN VIPER-Pathway analysis
#####################
organoid_Hallmarks <- viper(Overall_Bladder_zscore,hallmark_regulon,method="scale")



#####################
#VISUALIZE PATHWAY ANALYSIS:
#####################

#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(organoid_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:30],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.3,mar=c(5,10))



####################################################################################################################
#VIPER
####################################################################################################################


library(viper)
library(dorothea)

View(dorothea_hs)

allSamples_viper<-run_viper(Overall_Bladder_zscore,dorothea_hs)


#SELECTED GENES:

basal_luminal_genes<-c("SMAD3","FOXA1","GATA3","E2F1","FOXM1","TOP2A","ESR1","CENPF","E2F2","E2F3","ZEB1", "CDH1", "SNAI1","KMT2C","KMT2D")

#GENE EXPRESSION:
basal_luminal_overlap<-basal_luminal_genes[which(basal_luminal_genes %in% rownames(Overall_Bladder_zscore))]
heatmap.2(Overall_Bladder_zscore[basal_luminal_overlap,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Diff  expression",key.xlab = "z-score",cexRow = 0.8,cexCol = 0.3,mar=c(5,12))



#VIPER:
basal_luminal_overlap_viper<-basal_luminal_genes[which(basal_luminal_genes %in% rownames(allSamples_viper))]
heatmap.2(allSamples_viper[basal_luminal_overlap_viper,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "TF-Activity",key.xlab = "z-score",cexRow = 0.8,cexCol = 0.3,mar=c(5,12))

top_tfs<-names(sort(apply(allSamples_viper,1,var),decreasing=T))


heatmap.2(allSamples_viper[top_tfs[1:20],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "TF-Activity",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.3,mar=c(5,12))

write.csv(allSamples_viper,file="06 Excel - All-samples VIPER.csv")
  
}
  


####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########      ########      ########     ########             #########           ##       ############      ##
##     ##     ##     ##     ##     ##           ##     ##            ##      ##         ####           ##          ####
##     ##     ##     ##    ##      ##           ##     ##            ##       ##       ##  ##          ##         ##  ##
##     ########      #######       ########     ########             ##        ##     ##    ##         ##        ##    ##
##     ##            ##    ##      ##           ##                   ##       ##     ##########        ##       ##########
##     ##            ##     ##     ##           ##                   ##      ##     ##        ##       ##      ##        ##
##     ##            ##      ##    ########     ##                   #########     ##          ##      ##     ##          ##
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#DATE:   3-5-25
#GOAL:   NEW COMBINE CanFam3-1 (doesn't include OR.A-E and so this one is primary for the EXCEL File they want)
#PLAN:

setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 42 bladder cancer")


load("Combined_Bladder_CanFam3-1.RData")


library("gplots")
library("RColorBrewer")



####################################################################################################################################################################################################################################################################
#0 PREPARE COMBINED BATCHES:  loaded above
#################################################################################################################################################################################################################################################################### 

{
  load("Bladder40_Organoid_CanFam3-1.RData")
  samples40<-BladderOrganoid_tpm
  
  
  load("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/04 10 paired organoid-tissue/BladderOrganoid_paired_CanFam-3-1.RData")
  samples14<-BladderOrganoid_tpm
  new_colnames<-paste(paste(sapply(strsplit(colnames(BladderOrganoid_tpm),split="_"),function(x) x[2]),
                            sapply(strsplit(colnames(BladderOrganoid_tpm),split="_"),function(x) x[1]),sep="_"),
                      sapply(strsplit(colnames(BladderOrganoid_tpm),split="_"),function(x) x[3]),sep="_")
  
  colnames(samples14)<-new_colnames
  
  overlap_genes<-intersect(rownames(samples14),rownames(samples40))
  
  Overall_Bladder_tpm<-cbind(samples40[overlap_genes,],samples14[overlap_genes,])
  
  #FILTER OUT ALL ZERO ROWS:
  zero_rows<-which(rowSums(Overall_Bladder_tpm)==0)
  organoid_logRPM_zero<-Overall_Bladder_tpm[-zero_rows,]
  
  Overall_Bladder_tpm<-organoid_logRPM_zero
  
  
  Overall_Bladder_zscore<-t(scale(t(Overall_Bladder_tpm)))
  
  
  #####################################################
  #NORMALIZE TISSUE AND ORGANOID SEPARALY
  #####################################################
  organoid_zscore<-t(scale(t(Overall_Bladder_tpm[,grep("organoid",colnames(Overall_Bladder_tpm))]),center=T,scale=T))
  tissue_zscore<-t(scale(t(Overall_Bladder_tpm[,grep("issue",colnames(Overall_Bladder_tpm))]),center=T,scale=T))
  
  Overall_organ_tissue_zscore_raw<-cbind(organoid_zscore,tissue_zscore)
  
  Overall_organ_tissue_zscore<-Overall_organ_tissue_zscore_raw[-which(is.na(rowSums(Overall_organ_tissue_zscore_raw))==T),]
  
  
  write.csv(log10(Overall_Bladder_tpm+1),file="Overall_Bladder_logTPM.csv")
  
  save(list=c("Overall_Bladder_tpm","Overall_Bladder_zscore","Overall_organ_tissue_zscore"),file="Combined_Bladder_CanFam3-1.RData")
  
}


###############################################################################################################################################################
#PART 1:   ORGANOIDS VS TISSUE
###############################################################################################################################################################

{
#####################################################
#FIGURE 2A:  CORR-MATRIX.    Published analysis
#####################################################

organoid_rpm_corr<-cor(Overall_Bladder_tpm)

heatmap.2(organoid_rpm_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "rpm corr",key.xlab = "rpm corr",cexRow = 0.4,cexCol = 0.6)

#####################################################
#FIGURE 2A:   CORR-MATRIX.   Refined analysis
#####################################################

organoid_zscore_corr<-cor(Overall_Bladder_zscore)

heatmap.2(organoid_zscore_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "z-score corr",key.xlab = "z-score corr",cexRow = 0.3,cexCol = 0.6)


organoid_zscore_corr<-cor(Overall_organ_tissue_zscore)

heatmap.2(organoid_zscore_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "z-score corr",key.xlab = "z-score corr",cexRow = 0.3,cexCol = 0.6)





#####################################################
#PCA:  all
#####################################################

#RUN PCA
organoid_pca<-prcomp(t(Overall_Bladder_zscore))
organoid_pca_data<-organoid_pca$x
PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)

plot(organoid_pca_data[,1],organoid_pca_data[,2],col="black",pch=16,
     xlab="PC1 (29%):  tissue vs organoid",ylab="PC2 (15%): tissue v organoid", main="PCA:  all samples")
text(organoid_pca_data[,1],organoid_pca_data[,2],rownames(organoid_pca_data),col="black",cex=0.3)



#####################################################
#PCA:  organoids only
#####################################################

#RUN PCA
organoid_pca<-prcomp(t(Overall_organ_tissue_zscore[,grep("ganoid",colnames(Overall_organ_tissue_zscore))]))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)
table(sapply(strsplit(sample_names_vector,"_"), function(x) x[2]))
#Bo    boo    liz  molly Murphy  Ohren  Pearl  Penny  pilot  sadie Sadie2 Sadie3 Sadie4 Sadie7    Tux 
#2      3      1      1      3      3      3      5      1      2      2      6      3      3      3 

name2color<-c("darkblue","blue","red","darkred","purple","darkorange","yellow4","cyan4","pink3","darkgreen","black")
names(name2color)<-c("Bo","boo","liz","molly","Murphy","Ohren","Pearl","Penny","pilot","adie","Tux")

for (idx in 1:length(sample_names_vector)){
  sample_names_vector[grep(names(name2color)[idx],sample_names_vector)]<-name2color[idx]
}

tmp_names<-gsub("_rep","",gsub("organoid_","",rownames(organoid_pca_data)))


plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16,
     xlab="PC1 (24%):  tissue vs organoid",ylab="PC2 (17%): Tux vs all",main="PCA:  Organoids only")
text(organoid_pca_data[,1],organoid_pca_data[,2],tmp_names,col=sample_names_vector,cex=0.3)

plot(organoid_pca_data[,1],organoid_pca_data[,3],col=sample_names_vector,pch=16,
     xlab="PC1 (24%):  tissue vs organoid",ylab="PC2 (11%):",main="PCA:  Organoids only")
text(organoid_pca_data[,1],organoid_pca_data[,3],tmp_names,col=sample_names_vector,cex=0.3)


#####################################################
#PCA:  organoids only
#####################################################

#RUN PCA
organoid_pca<-prcomp(t(Overall_organ_tissue_zscore[,grep("issue",colnames(Overall_organ_tissue_zscore))]))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)


#MAKE COLOR VECTOR:
sample_names_vector<-rownames(organoid_pca_data)
table(sapply(strsplit(sample_names_vector,"_"), function(x) x[2]))
#boo  Penny  pilot  sadie sadie3    Tux 
#2      6      1      1      2      1 

name2color<-c("blue","cyan4","pink3","darkgreen","black")
names(name2color)<-c("boo","Penny","pilot","adie","Tux")

for (idx in 1:length(sample_names_vector)){
  sample_names_vector[grep(names(name2color)[idx],sample_names_vector)]<-name2color[idx]
}

tmp_names<-gsub("_rep","",gsub("tissue_","",rownames(organoid_pca_data)))


plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16,
     xlab="PC1 (47%):  ",ylab="PC2 (12%): ",main="PCA:  Tissue Only")
text(organoid_pca_data[,1],organoid_pca_data[,2],tmp_names,col=sample_names_vector,cex=0.3)




#####################################################
#heatmaps
#####################################################

#RANK PATHWAYS BASED ON VARIANCE:
genes_in_order<-names(sort(apply(Overall_Bladder_tpm,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_zscore[genes_in_order[1:2000],]),col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.3,cexCol = 0.1,mar=c(5,10))





#####################################################
#FIGURE 3:   PATHWAY ANALYSIS
#####################################################

#####################
#LOAD PACKAGES
#####################


library(viper)
library(msigdbr)



#####################
#2 MAKE HALLMARK-OBJECT FOR PATHWAY ANALYSIS
#####################

#EXTRACT MOUSE-HALLMARKS:
HALLMARK_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
HALLMARK_gene_sets$gs_name<-gsub("HALLMARK_","",HALLMARK_gene_sets$gs_name)


#CONVERT MOUSE-HALLMARK DATA-FRAME TO A "REGULON"
msigdbr_list = split(x = HALLMARK_gene_sets$gene_symbol, f = HALLMARK_gene_sets$gs_name)
hallmark_regulon<-list()

for (tf in names(msigdbr_list)){
  sig_genes<-msigdbr_list[[tf]]
  
  hallmark_regulon[[tf]]<-list()
  tmp_tfmode<-rep(1,length(sig_genes))
  names(tmp_tfmode)<-sig_genes
  
  hallmark_regulon[[tf]]$'tfmode'<-tmp_tfmode
  hallmark_regulon[[tf]]$'likelihood'<-as.numeric(tmp_tfmode)
}

#####################
#RUN VIPER-Pathway analysis
#####################
organoid_Hallmarks <- viper(Overall_Bladder_zscore,hallmark_regulon,method="scale")



#####################
#VISUALIZE PATHWAY ANALYSIS:
#####################

#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(organoid_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:30],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.3,mar=c(5,10))



####################################################################################################################
#VIPER
####################################################################################################################


library(viper)
library(dorothea)

View(dorothea_hs)

allSamples_viper<-run_viper(Overall_Bladder_zscore,dorothea_hs)


#SELECTED GENES:

basal_luminal_genes<-c("SMAD3","FOXA1","GATA3","E2F1","FOXM1","TOP2A","ESR1","CENPF","E2F2","E2F3","ZEB1", "CDH1", "SNAI1","KMT2C","KMT2D")

#GENE EXPRESSION:
basal_luminal_overlap<-basal_luminal_genes[which(basal_luminal_genes %in% rownames(Overall_Bladder_zscore))]
heatmap.2(Overall_Bladder_zscore[basal_luminal_overlap,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Diff  expression",key.xlab = "z-score",cexRow = 0.8,cexCol = 0.3,mar=c(5,12))



#VIPER:
basal_luminal_overlap_viper<-basal_luminal_genes[which(basal_luminal_genes %in% rownames(allSamples_viper))]
heatmap.2(allSamples_viper[basal_luminal_overlap_viper,],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "TF-Activity",key.xlab = "z-score",cexRow = 0.8,cexCol = 0.3,mar=c(5,12))

top_tfs<-names(sort(apply(allSamples_viper,1,var),decreasing=T))


heatmap.2(allSamples_viper[top_tfs[1:20],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "TF-Activity",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.3,mar=c(5,12))

write.csv(allSamples_viper,file="06 Excel - All-samples VIPER.csv")

}


####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########      ########      ########     ########             #########           ##       ############      ##
##     ##     ##     ##     ##     ##           ##     ##            ##      ##         ####           ##          ####
##     ##     ##     ##    ##      ##           ##     ##            ##       ##       ##  ##          ##         ##  ##
##     ########      #######       ########     ########             ##        ##     ##    ##         ##        ##    ##
##     ##            ##    ##      ##           ##                   ##       ##     ##########        ##       ##########
##     ##            ##     ##     ##           ##                   ##      ##     ##        ##       ##      ##        ##
##     ##            ##      ##    ########     ##                   #########     ##          ##      ##     ##          ##
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######    ##     ##   ########    #########   #######                                                               
##     ##          ##    ##          ##     ##   ##     ##   ##         ##                           
##     ##          ##   ##           ##     ##   ##     ##   ##         ##                          
##     ########    ##   ##   ####    ##     ##   ########    #######     #######                                                
##     ##          ##   ##      ##   ##     ##   ##    ##    ##                ##                          
##     ##          ##    ##     ##   ##     ##   ##     ##   ##                ##                   
##     ##          ##     #######     #######    ##      ##  ########    #######                                                 
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#DATE:   3-5-25
#GOAL:   COMBINE BATCHES (same as prev section but at samples OR.A, OR.E, OR.B)
#PLAN:

setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 42 bladder cancer")

load("BATCH-Combined_Bladder_CanFam3-1.RData")


library("gplots")
library("RColorBrewer")



####################################################################################################################################################################################################################################################################
#0 PREPARE COMBINED BATCHES:  loaded above
#################################################################################################################################################################################################################################################################### 

{

#######################################
#COMBINE DATA:  
#######################################
  
load("Bladder40_Organoid_CanFam3-1.RData")
samples40<-BladderOrganoid_tpm


load("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/04 10 paired organoid-tissue/BladderOrganoid_paired_CanFam-3-1.RData")
samples14<-BladderOrganoid_tpm

load("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/01 analysis - blca-VIPER/blca_organoid_tpm.RData")
colnames(blca_organoid_tpm)<-gsub("_scaledTPM","",colnames(blca_organoid_tpm))

#NEW NAMES:
tmp_converter<-read.csv("05 BATCHES - combine QC/new_sample_convert/ALL_samples_convert.csv",header=T,as.is=T)
old2new<-tmp_converter$Eugene.ID.to.USE
names(old2new)<-tmp_converter$old.name

new2batch<-tmp_converter$batch
names(new2batch)<-tmp_converter$Eugene.ID.to.USE

#15800 for CamFam3-1 (14,437 if I add OR.A-E samples)
overlap_genes<-intersect(intersect(rownames(samples14),rownames(samples40)),rownames(blca_organoid_tpm))

Overall_Bladder_tpm<-cbind(samples40[overlap_genes,],samples14[overlap_genes,],blca_organoid_tpm[overlap_genes,])
dim(Overall_Bladder_tpm)

#MISSING:  "organoid_Sadie2_urine_rep1"   was RO-70 and looked at original FASTQ files and wasn't there
Overall_Bladder_tpm_paperALL<-Overall_Bladder_tpm[,which(colnames(Overall_Bladder_tpm) %in% names(old2new))]
dim(Overall_Bladder_tpm_paperALL)

colnames(Overall_Bladder_tpm_paperALL)<-old2new[colnames(Overall_Bladder_tpm_paperALL)]

batches_ordered<-new2batch[colnames(Overall_Bladder_tpm_paperALL)]

#missing samples: RO-70 (no fastq), RO85-6 (no fastq)
unname(old2new)[-which(old2new %in% colnames(Overall_Bladder_tpm_paperALL))]

#######################################
#QUANTILE NORMALIZE 
#######################################


library(limma)

Overall_Bladder_logTPM_batch <- removeBatchEffect(log10(Overall_Bladder_tpm_paperALL+1),batches_ordered)




#######################################
#Z-score
#######################################


#FILTER OUT ALL ZERO ROWS:
zero_rows<-which(rowSums(Overall_Bladder_logTPM_batch)==0)
organoid_logRPM_zero<-Overall_Bladder_logTPM_batch[-zero_rows,]

Overall_Bladder_tpm<-organoid_logRPM_zero


Overall_Bladder_zscore<-t(scale(t(Overall_Bladder_tpm)))



save(list=c("Overall_Bladder_tpm","Overall_Bladder_zscore"),file="BATCH-Combined_Bladder_CanFam3-1.RData")


}

###############################################################################################################################################################
#1:   CHECK ALL SAMPLES
###############################################################################################################################################################


#####################################################
#FIGURE 2A:  CORR-MATRIX.    Published analysis
#####################################################

organoid_rpm_corr<-cor(Overall_Bladder_tpm)

heatmap.2(organoid_rpm_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "rpm corr",key.xlab = "rpm corr",cexRow = 0.4,cexCol = 0.6)

#####################################################
#FIGURE 2A:   CORR-MATRIX.   Refined analysis
#####################################################

organoid_zscore_corr<-cor(Overall_Bladder_zscore)

heatmap.2(organoid_zscore_corr,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "z-score corr",key.xlab = "z-score corr",cexRow = 0.3,cexCol = 0.6)




#####################################################
#PCA:  all
#####################################################

library(Seurat)

color_vector <- scales::hue_pal()(14)
name2color<-color_vector[14:1]
names(name2color)<-c("P01","P02","P03","P05","P06","P07","P08","P09","P10","P11","P12","P13","P14")



#MAKE COLOR VECTOR:
sample_names_vector<-colnames(Overall_Bladder_zscore)
table(sapply(strsplit(sample_names_vector,"_"), function(x) x[1]))
#P01 P02 P03 P05 P06 P07 P08 P09 P10 P11 
#1   1   1   1  18   5   1   2   3   3 


for (idx in 1:length(sample_names_vector)){
  sample_names_vector[grep(names(name2color)[idx],sample_names_vector)]<-name2color[idx]
}



#RUN PCA
organoid_pca<-prcomp(t(Overall_Bladder_zscore))
organoid_pca_data<-organoid_pca$x
PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)

plot(organoid_pca_data[,1],organoid_pca_data[,2],col=sample_names_vector,pch=16,
     xlab="PC1 (31%):  tissue vs organoid",ylab="PC2 (16%): tissue vs organoid", main="PCA:  all samples")
text(organoid_pca_data[,1],organoid_pca_data[,2],rownames(organoid_pca_data),col=sample_names_vector,cex=0.3)



#####################################################
#PCA:  organoids only
#####################################################

#RUN PCA
organoid_pca<-prcomp(t(Overall_Bladder_zscore[,grep("Org",colnames(Overall_Bladder_zscore))]))
organoid_pca_data<-organoid_pca$x

PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)

org_names_vector<-rownames(organoid_pca_data)
for (idx in 1:length(org_names_vector)){
  org_names_vector[grep(names(name2color)[idx],org_names_vector)]<-name2color[idx]
}


tmp_names<-gsub("_rep","",gsub("organoid_","",rownames(organoid_pca_data)))


plot(organoid_pca_data[,1],organoid_pca_data[,2],col=org_names_vector,pch=16,
     xlab="PC1 (22%):"  ,ylab="PC2 (17%):",main="PCA:  Organoids only")
text(organoid_pca_data[,1],organoid_pca_data[,2],tmp_names,col=org_names_vector,cex=0.3)

plot(organoid_pca_data[,1],organoid_pca_data[,3],col=sample_names_vector,pch=16,
     xlab="PC1 (22%):  ",ylab="PC2 (12%):",main="PCA:  Organoids only")
text(organoid_pca_data[,1],organoid_pca_data[,3],tmp_names,col=sample_names_vector,cex=0.3)






#####################################################
#heatmaps
#####################################################

#RANK PATHWAYS BASED ON VARIANCE:
genes_in_order<-names(sort(apply(Overall_Bladder_tpm,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_zscore[genes_in_order[1:2000],]),col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.3,cexCol = 0.1,mar=c(5,10))





#####################################################
#FIGURE 3:   PATHWAY ANALYSIS
#####################################################

#####################
#LOAD PACKAGES
#####################

#BiocManager::install("viper")

library(viper)
library(msigdbr)



#####################
#2 MAKE HALLMARK-OBJECT FOR PATHWAY ANALYSIS
#####################

#EXTRACT MOUSE-HALLMARKS:
HALLMARK_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
HALLMARK_gene_sets$gs_name<-gsub("HALLMARK_","",HALLMARK_gene_sets$gs_name)


#CONVERT MOUSE-HALLMARK DATA-FRAME TO A "REGULON"
msigdbr_list = split(x = HALLMARK_gene_sets$gene_symbol, f = HALLMARK_gene_sets$gs_name)
hallmark_regulon<-list()

for (tf in names(msigdbr_list)){
  sig_genes<-msigdbr_list[[tf]]
  
  hallmark_regulon[[tf]]<-list()
  tmp_tfmode<-rep(1,length(sig_genes))
  names(tmp_tfmode)<-sig_genes
  
  hallmark_regulon[[tf]]$'tfmode'<-tmp_tfmode
  hallmark_regulon[[tf]]$'likelihood'<-as.numeric(tmp_tfmode)
}

#####################
#RUN VIPER-Pathway analysis
#####################
organoid_Hallmarks <- viper(Overall_Bladder_zscore,hallmark_regulon,method="scale")



#####################
#VISUALIZE PATHWAY ANALYSIS:
#####################

#RANK PATHWAYS BASED ON VARIANCE:
Hallmarks_in_order<-names(sort(apply(organoid_Hallmarks,1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:30],],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.3,mar=c(5,10))



###############################################################################################################################################################
#2:   PAPER FIGURES
###############################################################################################################################################################

hist(Overall_Bladder_tpm,breaks=100)

#LOAD SAMPLE SETS:
fig1_samples<-read.csv("05b BATCHES - paper figs/new Chris directions/Fig01_all-samples.csv",as.is=T,header=T)[,1]
fig1_samples[-which(fig1_samples %in% colnames(Overall_Bladder_tpm))]#RO-85 and RO-86 not in fastq = "P12_Tis_Bladder_R1" "P12_Tis_Bladder_R2"
fig1_samples_filter<-fig1_samples[which(fig1_samples %in% colnames(Overall_Bladder_tpm))]


fig2_samples<-read.csv("05b BATCHES - paper figs/new Chris directions/Fig02_organoids.csv",as.is=T,header=T)[,1]
fig2_samples[-which(fig2_samples %in% colnames(Overall_Bladder_tpm))]

fig9_samples<-read.csv("05b BATCHES - paper figs/new Chris directions/Fig09_all-sadie.csv",as.is=T,header=T)[,1]
fig9_samples[-which(fig9_samples %in% colnames(Overall_Bladder_tpm))]#Missing RO-70 (no fastq) = "P06_Org_V2_Urine_R1"
fig9_samples_filter<-fig9_samples[which(fig9_samples %in% colnames(Overall_Bladder_tpm))]


fig10_samples<-read.csv("05b BATCHES - paper figs/new Chris directions/Fig10_sadie-organoids.csv",as.is=T,header=T)[,1]
fig10_samples_filter<-fig10_samples[which(fig10_samples %in% colnames(Overall_Bladder_tpm))]#Missing RO-70 (no fastq) = "P06_Org_V2_Urine_R1"

#GENE SETS:
fig4_genes<-c("GATA3", "TP63", "PAX8", "KRT5", "KRT7", "CD44", "UPK1A", "PCP4L1", "KRT20", "FOXA1","KRT75","KRT18","KRT14")
fig4_gene_overlap<-fig4_genes[which(fig4_genes %in% rownames(Overall_Bladder_tpm))]

fig5_genes<-c("PTGER2", "ERBB2", "IRGM", "IL6", "IL1B", "CSF2", "CCND1", "TP53", "TNF", "FOXM1", "HGF", "E2F3", "FOXO1", "NFKB1", "IFNG", "LET7", "RABL6", "JUN", "IL1A", "E2F1", "JUN", "IL22", "VEGFA", "MAP2K3", "RELA", "MAPK3", "EGFR", "IFI16", "E2F2", "ZFP36")
fig5_gene_overlap<-fig5_genes[which(fig5_genes %in% rownames(Overall_Bladder_tpm))]
  
fig6_genes<-c("SOX12", "RGMA", "PIK3CD", "LURAP1", "MAPK8IP3", "S100A13", "TGFB1", "SLC38A5", 
              "CLMP", "PPP1R1B", "ZNF503", "ANOS1", "CA9", "TGM1", "SULF2", "FOLR2", "FNBP1", 
              "TSHZ2", "CCNJL", "PTPRM", "LGR6", "MATN4", "SP8", "PGR", "CXCL8", "PAQR6", 
              "GJA1", "DHRS9", "ACKR3", "DKK3", "PDPN", "CDH13", "SPRY1", "SPP1", "APP", 
              "KRT17", "STOM")
fig6_gene_overlap<-fig6_genes[which(fig6_genes %in% rownames(Overall_Bladder_tpm))]




#########################################################
#FIG 1:  PCA:  tissue vs org
#########################################################

fig1_colors<-fig1_samples_filter

for (idx in 1:length(fig1_colors)){
  fig1_colors[grep(names(name2color)[idx],fig1_samples_filter)]<-name2color[idx]
}



#RUN PCA
organoid_pca<-prcomp(t(Overall_Bladder_zscore[,fig1_samples_filter]))
organoid_pca_data<-organoid_pca$x
PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)

plot(organoid_pca_data[,1],organoid_pca_data[,2],col=fig1_colors,pch=16,
     xlab="PC1 (33%):  tissue vs organoid",ylab="PC2 (16%): patient replicates", main="FIG 1:  PCA:  all samples")
text(organoid_pca_data[,1],organoid_pca_data[,2],rownames(organoid_pca_data),col=fig1_colors,cex=0.3)




#########################################################
#FIG 2:  PCA:  organoids
#########################################################
fig2_colors<-fig2_samples

for (idx in 1:length(fig2_colors)){
  fig2_colors[grep(names(name2color)[idx],fig2_samples)]<-name2color[idx]
}



#RUN PCA
organoid_pca<-prcomp(t(Overall_Bladder_zscore[,fig2_samples]))
organoid_pca_data<-organoid_pca$x
PoV <- round(organoid_pca$sdev^2/sum(organoid_pca$sdev^2),2)*100
names(PoV)<-colnames(organoid_pca$x)

plot(organoid_pca_data[,1],organoid_pca_data[,3],col=fig2_colors,pch=16,
     xlab="PC1 (25%):  patients",ylab="PC3 (14%): patient replicates", main="FIG 2:  PCA:  organoids")
text(organoid_pca_data[,1],organoid_pca_data[,3],rownames(organoid_pca_data),col=fig2_colors,cex=0.3)






#########################################################
#FIG 3:  HEATMAP:  top 200 var genes
#########################################################


#RANK PATHWAYS BASED ON VARIANCE:
genes_in_order<-names(sort(apply(Overall_Bladder_tpm[,fig1_samples_filter],1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_zscore[genes_in_order[1:2000],fig1_samples_filter]),RowSideColors = fig1_colors,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.5,cexCol = 0.1,mar=c(5,10))





#########################################################
#FIG 4:  HEATMAP Z-SCORE small number genes (GATA3)
#########################################################

#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_zscore[fig4_gene_overlap,fig1_samples_filter]),RowSideColors = fig1_colors,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 0.5,cexCol = 1,mar=c(5,10))



#########################################################
#FIG 5:  HEATMAP TPM small number genes (IL6)
#########################################################

Overall_Bladder_tpm[which(Overall_Bladder_tpm<0)]<-0

#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_tpm[fig5_gene_overlap,fig1_samples_filter]),RowSideColors = fig1_colors,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "logTPM",
          cexRow = 0.5,cexCol = 0.7,mar=c(5,10))






#########################################################
#FIG 6:  HEATMAP TPM small number genes (SOX2)
#########################################################

#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_tpm[fig6_gene_overlap,fig1_samples_filter]),RowSideColors = fig1_colors,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "logTPM",
          cexRow = 0.5,cexCol = 0.7,mar=c(5,10))


#########################################################
#FIG 7:  HEATMAP TPM small number genes (SOX2)
#########################################################

#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_tpm[fig6_gene_overlap,fig2_samples]),RowSideColors = fig2_colors,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "logTPM",
          cexRow = 0.5,cexCol = 0.7,mar=c(5,10))



#########################################################
#FIG 8:  HEATMAP PATHWAY:   tissue vs org
#########################################################


#PLOT TOP-20 PATHWAYS:
heatmap.2(organoid_Hallmarks[Hallmarks_in_order[1:30],fig1_samples_filter],ColSideColors = fig1_colors,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "enrichment score",key.xlab = "z-score",cexRow = 0.6,cexCol = 0.6,mar=c(5,10))




#########################################################
#FIG 9:  HEATMAP SADIE
#########################################################
fig9_colors<-fig9_samples_filter

for (idx in 1:length(fig2_colors)){
  fig9_colors[grep(names(name2color)[idx],fig9_samples_filter)]<-name2color[idx]
}



#RANK PATHWAYS BASED ON VARIANCE:
genes_in_order<-names(sort(apply(Overall_Bladder_tpm[,fig9_samples_filter],1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_zscore[genes_in_order[1:2000],fig9_samples_filter]),RowSideColors = fig9_colors,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),Rowv = NA,,
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
          cexRow = 1,cexCol = 0.1,mar=c(5,10))


#########################################################
#FIG 10:  HEATMAP SADIE  - organoid
#########################################################
fig10_colors<-fig10_samples_filter

for (idx in 1:length(fig2_colors)){
  fig10_colors[grep(names(name2color)[idx],fig10_samples_filter)]<-name2color[idx]
}



#RANK PATHWAYS BASED ON VARIANCE:
genes_in_order<-names(sort(apply(Overall_Bladder_tpm[,fig10_samples_filter],1,var),decreasing=T))


#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_zscore[genes_in_order[1:2000],fig10_samples_filter]),RowSideColors = fig10_colors,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),Rowv = NA,,
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "z-score",
            cexRow = 1,cexCol = 0.1,mar=c(5,10))

#########################################################
#FIG 11:  HEATMAP SADIE  - organoid
#########################################################




#PLOT TOP-20 PATHWAYS:
heatmap.2(t(Overall_Bladder_tpm[c("CD274","PDCD1"),fig10_samples_filter]),RowSideColors = fig10_colors,
          col=colorRampPalette(c("blue", "white", "red"))(n = 20),Rowv = NA,,
          density.info = "none",trace="none",key.title = "Gene Expression",key.xlab = "logTPM",
          cexRow = 0.5,cexCol = 1,mar=c(5,10))

####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########      ########    ########      ##     ##   #######    ######   ##########    #######                                  
##     ##     ##     ##         ##      ##     ##     ##   ##        ##            ##       ##                
##     ##    ##      ##         ##      ##     ##     ##   ##        ##            ##       ##                
##     #######       ########   ##      ##     ##     ##   #######    ######       ##        #######                         
##     ##    ##      ##         ##   ## ##     ##     ##   ##              ##      ##              ##          
##     ##     ##     ##         ##     ####    ##     ##   ##              ##      ##              ##          
##     ##      ##    ########    ######## ###   #######    ########   ######       ##        #######                                 
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#CHRIS WASNT VOLCANO PLOTS:    urine vs tissue


{
setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 42 bladder cancer")


###############################################################################################################################################################
#3:   CHRIS' requested Volcano plots:  can't use DESeq with batch corrected data (use integrated TPMs)
###############################################################################################################################################################

load("BATCH-Combined_Bladder_CanFam3-1.RData")
Overall_Bladder_tpm[which(Overall_Bladder_tpm<0)]<-0
hist(Overall_Bladder_tpm,breaks=100)

sample_maps<-read.csv("05c BATCHES - volcanos/FINAL_Eugene_final_figs_TCC.csv")

#CAN"T use DESEQ because i need counts (can't use with batch corrected data)
volcano_loop_ttest<-function(data,urine,tissue){
  
  
  data<-data[which(apply(data,1,var)>0),]
  
  #DEFINE Tissues:
  #data=Tab1_KI_RNAseq_logRPM
  #lineage1<-"Tab1KI_OIR"
  #lineage2<-"wt_OIR"
  
  #define column #'s for two lineages:
  lineage1_idx<-which(colnames(data) %in% urine)
  lineage2_idx<-which(colnames(data) %in% tissue)
  
  
  #MAKE A MATRIX with P-value & Fold change for every gene
  diff_Express_matrix<-matrix(NA,ncol=2,nrow=nrow(data))
  colnames(diff_Express_matrix)<-c("log-FC","log-Pvalue")
  rownames(diff_Express_matrix)<-rownames(data)
  
  #FOR LOOP:  calculates fold-change and p-value ~20,000 times:
  for (gene in rownames(diff_Express_matrix)){
    #extract values
    lineage1_values<-data[gene,lineage1_idx]
    lineage2_values<-data[gene,lineage2_idx]
    
    #calculate/store fold change
    logFoldChange<-mean(lineage1_values)-mean(lineage2_values)
    diff_Express_matrix[gene,"log-FC"]<-logFoldChange
    
    #calculate/store p-value
    lineage_ttest<-t.test(lineage1_values,lineage2_values)
    diff_Express_matrix[gene,"log-Pvalue"]<-log10(lineage_ttest$p.value)
    
  }
  
  return(diff_Express_matrix)
}


#######################################
#ALL SAMPLES:
#######################################
urine_samples<-sample_maps$Eugene.ID.to.USE[which(sample_maps$Type=="urine")]
tissue_samples<-sample_maps$Eugene.ID.to.USE[which(sample_maps$Type=="tissue")]

diff_Express_matrix<-volcano_loop_ttest(Overall_Bladder_tpm,urine=urine_samples,tissue=tissue_samples)

write.csv(diff_Express_matrix,file="ALL - urine vs tissue.csv")

color_vector<-diff_Express_matrix[,"log-FC"]
color_vector[which(diff_Express_matrix[,"log-FC"] < 0)]<-"darkred"
color_vector[which(diff_Express_matrix[,"log-FC"] >0)]<-"darkblue"



#PLOT VOLCANO PLOT:
plot(diff_Express_matrix[,"log-FC"],-diff_Express_matrix[,"log-Pvalue"],
     col=color_vector,cex=0.5,pch=16,xlab="log-Fold Change (urine blue vs tissue red)",ylab="-logPvalue")
abline(h=1.3, lty=2)
#Add labels for genes
top_genes_idx<-intersect(which(abs(diff_Express_matrix[,"log-FC"])>0.2),which(diff_Express_matrix[,"log-Pvalue"]< (-1.3)))
top_genes<-rownames(diff_Express_matrix)[top_genes_idx]
text(diff_Express_matrix[top_genes,"log-FC"]+0.03,-diff_Express_matrix[top_genes,"log-Pvalue"],top_genes,cex=0.4,col=color_vector[top_genes])


#######################################
#ALL SAMPLES:
#######################################
urine_samples<-c("P06_Org_V4_Urine_R1","P06_Org_V4_Urine_R2","P06_Org_V4_Urine_R3")
tissue_samples<-c("P06_Org_V4_Bladder_R1","P06_Org_V4_Bladder_R2","P06_Org_V4_Bladder_R3")


diff_Express_matrix<-volcano_loop_ttest(Overall_Bladder_tpm[,c(urine_samples,tissue_samples)],urine=urine_samples,tissue=tissue_samples)

write.csv(diff_Express_matrix,file="Sadie - urine vs tissue.csv")


color_vector<-diff_Express_matrix[,"log-FC"]
color_vector[which(diff_Express_matrix[,"log-FC"] < 0)]<-"darkred"
color_vector[which(diff_Express_matrix[,"log-FC"] >0)]<-"darkblue"



#PLOT VOLCANO PLOT:
plot(diff_Express_matrix[,"log-FC"],-diff_Express_matrix[,"log-Pvalue"],
     col=color_vector,cex=0.5,pch=16,xlab="log-Fold Change (urine blue vs tissue red)",ylab="-logPvalue",
     xlim=c(-1.2,1))
abline(h=1.3, lty=2)
#Add labels for genes
top_genes_idx<-intersect(which(abs(diff_Express_matrix[,"log-FC"])>0.45),which(diff_Express_matrix[,"log-Pvalue"]< (-1.3)))
top_genes<-rownames(diff_Express_matrix)[top_genes_idx]
text(diff_Express_matrix[top_genes,"log-FC"]+0.06,-diff_Express_matrix[top_genes,"log-Pvalue"],top_genes,cex=0.4,col=color_vector[top_genes])


###############################################################################################################################################################
#3:  IF NEED DESEQ: CHRIS' requested Volcano plots
###############################################################################################################################################################

# Load required packages
library(DESeq2)
library(viper)
library(msigdbr)
library(fdrtool)

#LOAD COUNTS WITH OLD NAMES
load("Bladder40_Organoid_CanFam3-1.RData")
hist(log10(BladderOrganoid_counts+1),breaks=100)

BladderOrganoid_counts_round<-round(BladderOrganoid_counts,0)
rm(BladderOrganoid_counts)
rm(BladderOrganoid_tpm)

volcano_loop_deseq<-function(data,urine,tissue){
  
  tmp_matrix<-data[which(apply(data,1,var)>100),]
  
  
  sample_info <- data.frame(
    row.names = c(urine,tissue),
    condition = factor(c(rep("urine",length(urine)),rep("tissue",length(tissue))))  # First 4 are group A, next 4 are group B
  )
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = tmp_matrix, colData = sample_info, design = ~condition)
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)
  
  # Extract results with shrinkage applied to fold changes
  res <- results(dds, contrast = c("condition", "urine", "tissue"))
  
  # Compute FDR-adjusted p-value (q-value)
  res$QValue <- p.adjust(res$pvalue, method = "BH")
  
  # Compute False Discovery Rate (FDR) using fdrtool
  valid_pvals <- !is.na(res$pvalue)  # Remove NA values
  fdr_results <- fdrtool(res$pvalue[valid_pvals], statistic = "pvalue", plot = FALSE)
  res$FDR <- NA  # Initialize with NA
  res$FDR[valid_pvals] <- fdr_results$lfdr  # Store FDR values
  
  # Create result dataframe
  result_df <- data.frame(
    Gene = rownames(res),
    Log2FoldChange = res$log2FoldChange,
    PValue = res$pvalue,
    QValue = res$QValue,  # Adjusted p-value (q-value)
    FDR = res$FDR         # False Discovery Rate
  )
  
  # Remove NAs (optional)
  result_df <- na.omit(result_df)
  
  
  return(result_df)
}



#######################################
#ALL SAMPLES:
#######################################

#EXTRACT SAMPLE NAMES WITHIN ONE BATCH:
sample_maps<-read.csv("05c BATCHES - volcanos/FINAL_Eugene_final_figs_TCC.csv")
sample_converter<-read.csv('05b BATCHES - paper figs/new Chris directions/ALL_samples_convert.csv')
new2old<-sample_converter$old.name
names(new2old)<-sample_converter$Eugene.ID.to.USE
urine_samples<-new2old[sample_maps$Eugene.ID.to.USE[which(sample_maps$Type=="urine")]]
tissue_samples<-new2old[sample_maps$Eugene.ID.to.USE[which(sample_maps$Type=="tissue")]]
tissue_samples<-tissue_samples[which(tissue_samples %in% colnames(BladderOrganoid_counts_round))]#filter out those not from batch
urine_samples<-urine_samples[which(urine_samples %in% colnames(BladderOrganoid_counts_round))]


diff_Express_matrix<-volcano_loop_deseq(BladderOrganoid_counts_round[,c(urine_samples,tissue_samples)],
                                        urine=urine_samples,tissue=tissue_samples)


write.csv(diff_Express_matrix,file="ALL - urine vs tissue.csv")

rownames(diff_Express_matrix)<-diff_Express_matrix$Gene

color_vector<-diff_Express_matrix[,"Log2FoldChange"]
names(color_vector)<-diff_Express_matrix$Gene
color_vector[which(diff_Express_matrix[,"Log2FoldChange"] < 0)]<-"darkred"
color_vector[which(diff_Express_matrix[,"Log2FoldChange"] >0)]<-"darkblue"



#PLOT VOLCANO PLOT:
plot(diff_Express_matrix[,"Log2FoldChange"],-log(diff_Express_matrix[,"QValue"]),
     col=color_vector,cex=0.5,pch=16,xlab="log-Fold Change (urine blue vs tissue red)",ylab="-logPvalue",
     ylim=c(0,22))
abline(h=1.3, lty=2)
#Add labels for genes
top_genes_idx<-intersect(which(abs(diff_Express_matrix[,"Log2FoldChange"])>1.5),which(log(diff_Express_matrix[,"QValue"])< (-5)))
top_genes<-rownames(diff_Express_matrix)[top_genes_idx]
text(diff_Express_matrix[top_genes,"Log2FoldChange"]+0.3,-log(diff_Express_matrix[top_genes,"QValue"]),
     top_genes,cex=0.4,col=color_vector[top_genes])

length(which(diff_Express_matrix$QValue>0.05))/length(diff_Express_matrix$QValue)

#######################################
#ALL SAMPLES:
#######################################

tissue_samples<-c("organoid_Sadie3_tumor_rep1","organoid_Sadie3_tumor_rep2","organoid_Sadie3_tumor_rep3")
urine_samples<-c("organoid_Sadie3_urine_rep1","organoid_Sadie3_urine_rep2","organoid_Sadie3_urine_rep3")



diff_Express_matrix<-volcano_loop_deseq(BladderOrganoid_counts_round[,c(urine_samples,tissue_samples)],
                                        urine=urine_samples,tissue=tissue_samples)


write.csv(diff_Express_matrix,file="Sadie - urine vs tissue.csv")
rownames(diff_Express_matrix)<-diff_Express_matrix$Gene

color_vector<-diff_Express_matrix[,"Log2FoldChange"]
names(color_vector)<-diff_Express_matrix$Gene
color_vector[which(diff_Express_matrix[,"Log2FoldChange"] < 0)]<-"darkred"
color_vector[which(diff_Express_matrix[,"Log2FoldChange"] >0)]<-"darkblue"



#PLOT VOLCANO PLOT:
plot(diff_Express_matrix[,"Log2FoldChange"],-log(diff_Express_matrix[,"PValue"]),
     col=color_vector,cex=0.5,pch=16,xlab="log-Fold Change (urine blue vs tissue red)",ylab="-logPvalue")
#Add labels for genes
top_genes_idx<-intersect(which(abs(diff_Express_matrix[,"Log2FoldChange"])>1),which(log(diff_Express_matrix[,"QValue"])< (-100)))
top_genes<-rownames(diff_Express_matrix)[top_genes_idx]
text(diff_Express_matrix[top_genes,"Log2FoldChange"]+0.3,-log(diff_Express_matrix[top_genes,"QValue"]),
     top_genes,cex=0.4,col=color_vector[top_genes])

}


####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##     ########      ########    ########      ##     ##   #######    ######   ##########    #######                                  
##     ##     ##     ##         ##      ##     ##     ##   ##        ##            ##       ##                
##     ##    ##      ##         ##      ##     ##     ##   ##        ##            ##       ##                
##     #######       ########   ##      ##     ##     ##   #######    ######       ##        #######                         
##     ##    ##      ##         ##   ## ##     ##     ##   ##              ##      ##              ##          
##     ##     ##     ##         ##     ####    ##     ##   ##              ##      ##              ##          
##     ##      ##    ########    ######## ###   #######    ########   ######       ##        #######                                 
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#DATE: 7-16-24
#GOAL:  Karin wants PCA or Heatmap comparisons of Human MIBC TCGA vs Canine batch corrected organoids
#
#PLAN:   #1 load tcga (already subtyped)
#        #2 copy/paste canine-subtyping code so formatted exactly like tcga
#        #3 try PCA plot but also try heatmaps for comparisons
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################


setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 42 bladder cancer")


#LOAD NON-BATCH CORRECTED:
#load("Bladder40_Organoid_CanFam3-1.RData")
#BladderOrganoid_logTPM<-log10(BladderOrganoid_tpm+1)
#hist(BladderOrganoid_logTPM,breaks=100)


#LOAD ALL ORGANOID DATA:
load("BATCH-Combined_Bladder_CanFam3-1.RData")
Overall_Bladder_tpm[which(Overall_Bladder_tpm<0)]<-0
hist(Overall_Bladder_tpm,breaks=100)
dim(Overall_Bladder_tpm)

#LOAD TCGA and CANINE-MIBC DATA FRAMEWORK: 
load("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2025-03-10 - MIBC - POH framework/canine-TCGA_bladder.RData")
load("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2025-03-10 - MIBC - POH framework/BLCA_gsea_methods.RData")



############################################################################################################################################################################
#0 PREP DATA:  MIBC subtype
############################################################################################################################################################################

ordered_type<-sapply(strsplit(colnames(Overall_Bladder_tpm),split="_"), function(x) x[2])
organodi_idx<-which(ordered_type=="Org")
organoid_names<-colnames(Overall_Bladder_tpm)[organodi_idx]

library(consensusMIBC)

tmp_class<-getConsensusClass(Overall_Bladder_tpm,gene_id="hgnc_symbol")
sort(table(tmp_class$consensusClass))
tmp_class$simpleClass<-gsub("LumNS","Lum",gsub("LumU","Lum",gsub("LumP","Lum",gsub("Stroma-rich","Ba/Sq",tmp_class$consensusClass))))
table(tmp_class$simpleClass)
#  Ba/Sq   Lum 
#    6    49 

############################################################################################################################################################################
#1 PCA plots
############################################################################################################################################################################


#############################################################
#4.1 Define Drug targets
#############################################################


#LOAD DRUGBANK TARGETS:  
load("/Users/eugenedouglass/Dropbox/00 Independent work/--- 09 TEACHING/03 PHRM8210 - Biostats/2024/02 New apps/01 Multi-Omic Heatmap app/data_chem/binary_drug_targets.RData")
drug_targets<-colnames(binary_drug_targets_matrix)
rm(binary_drug_targets_matrix)

#LOAD TME-TARGETS:
load("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2020-21/2023-08-01 - CellPhone Sequel/00 papers/00 original papers/02 cellchat/CellChat-master/data/CellChatDB.human.rda")
View(CellChatDB.human$interaction)
receptors<-unique(c(CellChatDB.human$interaction$ligand,CellChatDB.human$complex$subunit_1,CellChatDB.human$complex$subunit_2,CellChatDB.human$complex$subunit_3,CellChatDB.human$complex$subunit_4))
rm(CellChatDB.human)





#############################################################
#4.2 HUMAN PCA PLOT:
#############################################################

human_RNAseq<-BLCA_tissues_RNAseq$Human$RNAseq

type2color<-BLCA_gsea_methods$MIBCtype$simple2color


#CANINE DATA:
human2subtype<-BLCA_gsea_methods$MIBCtype$tcga$simpleClass
names(human2subtype)<-rownames(BLCA_gsea_methods$MIBCtype$tcga)

#HUMAN INTL:  pca = 11%
human_pca_intl<-prcomp(t(human_RNAseq[which(rownames(human_RNAseq) %in% drug_targets),]))
human_pca_intl_data<-human_pca_intl$x
human_percent<-round(human_pca_intl$sdev^2/sum(human_pca_intl$sdev^2),2)*100

ordered_colors<-type2color[human2subtype[rownames(human_pca_intl_data)]]



plot(human_pca_intl_data[,1],human_pca_intl_data[,2],col=ordered_colors,pch=16,
     xlab="PC1 (16%):  ",ylab="PC2 (9%): ",main="408 Human MIBC Samples:  all drug-targets")
abline(h=0,v=0,lty=2,lwd=0.5)


#############################################################
#4.3 Canine PCA PLOT:
#############################################################

human_RNAseq<-BLCA_tissues_RNAseq$Canine56$RNAseq

type2color<-BLCA_gsea_methods$MIBCtype$simple2color


#CANINE DATA:
human2subtype<-BLCA_gsea_methods$MIBCtype$canine$simpleClass
names(human2subtype)<-rownames(BLCA_gsea_methods$MIBCtype$canine)

#HUMAN INTL:  pca = 11%
human_pca_intl<-prcomp(t(human_RNAseq[which(rownames(human_RNAseq) %in% drug_targets),]))
human_pca_intl_data<-human_pca_intl$x
human_percent<-round(human_pca_intl$sdev^2/sum(human_pca_intl$sdev^2),2)*100

ordered_colors<-type2color[human2subtype[rownames(human_pca_intl_data)]]



plot(-human_pca_intl_data[,1],-human_pca_intl_data[,2],col=ordered_colors,pch=16,
     xlab="PC1 (23%):  ",ylab="PC2 (8%): ",main="56 Canine MIBC Samples:  all drug-targets")
abline(h=0,v=0,lty=2,lwd=0.5)
#text(-human_pca_intl_data[,1],-human_pca_intl_data[,2],rownames(human_pca_intl_data))



#############################################################
#4.3 Organoid PCA PLOT:
#############################################################

human_RNAseq<-Overall_Bladder_tpm

type2color<-BLCA_gsea_methods$MIBCtype$simple2color


#CANINE DATA:
human2subtype<-tmp_class$simpleClass
names(human2subtype)<-rownames(tmp_class)

#HUMAN INTL:  pca = 11%
human_pca_intl<-prcomp(t(human_RNAseq[which(rownames(human_RNAseq) %in% drug_targets),]))
human_pca_intl_data<-human_pca_intl$x
human_percent<-round(human_pca_intl$sdev^2/sum(human_pca_intl$sdev^2),2)*100

ordered_colors<-type2color[human2subtype[rownames(human_pca_intl_data)]]



plot(-human_pca_intl_data[,1],-human_pca_intl_data[,2],col=ordered_colors,pch=16,
     xlab="PC1 (23%):  ",ylab="PC2 (8%): ",main="55 Canine MIBC Organoids:  all drug-targets")
abline(h=0,v=0,lty=2,lwd=0.5)
text(-human_pca_intl_data[,1],-human_pca_intl_data[,2],rownames(human_pca_intl_data),col=ordered_colors,cex=0.3)


############################################################################################################################################################################
#1 HEATMAP
############################################################################################################################################################################

library("gplots")
library("RColorBrewer")


############################
#CORRELATION MATRIX:
############################
canine29_RNAseq<-Overall_Bladder_zscore
tcga_RNAseq<-BLCA_tissues_RNAseq$Human$RNAseq_zscore

overlap_genes<-intersect(rownames(canine29_RNAseq),rownames(tcga_RNAseq))



Canine_classification_genes<-read.csv("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2025-03-10 - MIBC - POH framework/data/canine29/classifier_genes.csv",as.is=T)

overlap_genes_canineClass<-overlap_genes[which(overlap_genes %in% Canine_classification_genes$Gene.Symbol)]


canine_human_corr_canineClass<-cor(tcga_RNAseq[overlap_genes_canineClass,],canine29_RNAseq[overlap_genes_canineClass,],use="pairwise.complete.obs")

############################
#PLOT
############################

#MAKE COLOR BARS CONVERTIONS:
human2subtype<-BLCA_gsea_methods$MIBCtype$tcga$simpleClass
names(human2subtype)<-rownames(BLCA_gsea_methods$MIBCtype$tcga)
human_colors<-type2color[human2subtype[rownames(canine_human_corr_canineClass)]]


dog2subtype<-tmp_class$simpleClass
names(dog2subtype)<-rownames(tmp_class)
dog_colors<-type2color[dog2subtype[colnames(canine_human_corr_canineClass)]]


heatmap.2(t(canine_human_corr_canineClass),ColSideColors = human_colors,RowSideColors = dog_colors,col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title = "corr",key.xlab = "corr",cexRow = 0.4,cexCol = 0.2)


org_idx<-grep("_Org_",colnames(canine_human_corr_canineClass))
heatmap.2(t(canine_human_corr_canineClass[,org_idx]),ColSideColors = human_colors,RowSideColors = dog_colors[org_idx],col=colorRampPalette(c("blue", "white", "red"))(n = 20),
          density.info = "none",trace="none",key.title ="corr",key.xlab = "corr",cexRow = 0.4,cexCol = 0.2)


