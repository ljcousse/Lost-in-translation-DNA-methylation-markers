###############################################
###############################################
##                                           ##
##  Targeted Bisulphite sequencing analysis  ##
##                                           ##
###############################################
###############################################


## Location qMSP amplicon (obtained from Kim Lommen)
# NDRG4; 58497439-58497541
# GREM1; 33009978-33010025 
# LY75; 160761053-160761124

######################
# Parameter settings #
######################

options(warn=1L)
# mincov<-4 # Minimal coverage of CpG in all samples
# minperc<-1 # Percentage of samples that have to meet the minimal Coverage
# !!!!! If minperc < 1 you have to filter out rows with NA values for boxplot etc.
# minmethreads<-4 # Minimum number of (hydroxy)methylated reads in all samples (B-samples)
output_dir <- "/data/homes/louisc/Project_KimSmits/"
ct <- 0.01 # Constant value for calculation M-values from methylation degrees


#############################
# Loading required packages #
#############################

library(BiSeq)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(limma)
library(biomaRt)
library(openxlsx)
library(stringr)
library(precrec)
library(scales)
library(parallel)
library(survival)
library(biomaRt)
library(survminer)
library(gridExtra)
library(ggplot2)
library(data.table)

######################
# plotGenes function #
######################

source("Plots_AIC_function.R")

tiff(paste(output_dir,"LegendFirstPanel.jpg",sep=""),width=2048+512+256,height=1024*2/3)
par(mar=c(6,8,2,40))
plot(1,1,type="n")
multiplier <- 2
cex_legend <- 1*multiplier
col_gene <- c("cadetblue","blue","navy","lightgrey")
par(xpd=T)
legend(par("usr")[2]+((par("usr")[2]-par("usr")[1])/100*2.5),par("usr")[4],
       c(expression(paste(Delta,"AIC between base model and")),
         "   methylation model","","","CpG in primer","","Amplicons","","","Amplicons without sufficient data","",
         expression(paste(Delta,"AIC=0; Base model equally")),
         "probable as methylation model",
         expression(paste(Delta,"AIC=-6; Base only 0.05x as ")),
         "probable as methylation model"),
       pch=c(18,18,18,NA,4,NA,NA,NA,NA,NA,NA,NA,NA,NA),
       lty=c(NA,NA,NA,NA,NA,NA,1,1,1,1,NA,1,NA,2,NA),
       lwd=c(NA,NA,NA,NA,NA,NA,3,3,3,3,NA,3,NA,3,NA),
       col=c(col_gene[1:3],NA,"black",NA,col_gene,NA,"grey",NA,"grey",NA),cex=cex_legend)
par(xpd=F)
dev.off()

#########################
#########################
#########################
###                   ###
###                   ###   
###   Preprocessing   ###
###                   ###
###                   ###
#########################
#########################
#########################

#######################
# Importing Data data #
#######################

# Methylation data is extracted by bismark and saved as ".cov" files. The BiSeq package is used to load in the data for both runs seperatly. Data for both runs will be merged later on.
# Set working directory
wd_ini <- getwd()
setwd(paste(output_dir,"meth",sep=""))

# Get filenames of cov files
files <- list.files()[grep("CpG",list.files())]
files <- files[grep("cov$",files)]
# files <- files[grepl("run_231",files)]
files[grepl("run_239",files)] <- files[grepl("run_239",files)][sort(as.numeric(gsub("_.*$","",gsub(".*_S","",files[grepl("run_239",files)]))),index.return=T)$ix] # resort files for run 239
# print(files)

# Get sample annotation
sample_annot <- read.xlsx(paste(output_dir,"SampleSheet MUMC run.xlsx",sep=""),1,startRow=17,colNames=T)
sample_annot <- sample_annot[,1:2]
sample_annot[19,2] <- "111N_2057T_M31"
sample_annot <- cbind(sample_annot[,1],t(matrix(unlist(strsplit(sample_annot[,2],"_")),nrow=3)))
sample_annot[29,4] <- ""
colnames(sample_annot) <- c("Sample_ID","CRC","ccRCC","melanoma")
sample_annot <- data.frame(sample_annot,stringsAsFactors = F)
sample_annot$Sample_ID <- paste(rep("run231",nrow(sample_annot)),sample_annot$Sample_ID,sep="_")

# Get additional sample anotation melanoma follow up
sample_annot_2 <- read.xlsx(paste(output_dir,"Clinical follow-up data additional melanoma samples - Location-III.xlsx",sep=""),1,colNames=T)
sample_annot_2 <- sample_annot_2[-c(which(sample_annot_2$Sample.ID%in%c(21,19))),]
sample_annot_add <- cbind(paste(rep("run239",nrow(sample_annot_2)),sample_annot_2$Sample.ID,sep="_"),
                          rep("",nrow(sample_annot_2)),
                          rep("",nrow(sample_annot_2)),
                          sample_annot_2$Sample.nr.)
colnames(sample_annot_add) <- c("Sample_ID","CRC","ccRCC","melanoma")
sample_annot <- rbind(sample_annot,sample_annot_add)
sample_annot_prefilter <- sample_annot

# Remove run231_30 & run231_44 due to to low aligned reads (<1000)
sample_annot <- sample_annot[!(sample_annot$Sample_ID%in%c("run231_30","run231_44")),] 
files <- files[!(grepl("run_231_30_",files) | grepl("run_231_44_",files))]

# Get primer annotation
output_dir <- "/data/homes/louisc/Project_KimSmits/"
primer_annot <- read.xlsx(paste(output_dir,"SupplementaryTable1.xlsx",sep=""))
primer_annot[,1] <- gsub("FY","LY",primer_annot[,1])
primer_annot[,1] <- gsub(" ","",primer_annot[,1])
colnames(primer_annot) <- primer_annot[1,]

primer_annot_GREM1 <- primer_annot[grepl("GREM1_",primer_annot$name),]
primer_annot_NDRG4 <- primer_annot[grepl("NDRG4_",primer_annot$name),]
primer_annot_LY75 <- primer_annot[grepl("LY75_",primer_annot$name),]
for (i in 4:7){
  primer_annot_GREM1[,i] <- as.numeric(primer_annot_GREM1[,i])
  primer_annot_NDRG4[,i] <- as.numeric(primer_annot_NDRG4[,i])
  primer_annot_LY75[,i] <- as.numeric(primer_annot_LY75[,i])
}

primer_annot_GREM1
primer_annot_NDRG4
primer_annot_LY75

# CRC 
CRC_annot <- read.xlsx(paste(output_dir,"Klinische data location III CRC, ccRCC, melanoma- compleet 17022022.xlsx",sep=""),
                       1,startRow=1,colNames=T)
CRC_annot$Sample_ID <- paste(str_extract(CRC_annot[,1],"^[0-9]+"),CRC_annot[,2],sep="")

sum(sample_annot$CRC%in%CRC_annot$Sample_ID)
CRC_annot <- CRC_annot[CRC_annot$Sample_ID%in%sample_annot$CRC,]

# ccRCC
ccRCC_annot <- read.xlsx(paste(output_dir,"Klinische data location III CRC, ccRCC, melanoma- compleet 17022022.xlsx",sep=""),
                         2,startRow=1,colNames=T)
ccRCC_annot$Sample_ID <- paste(ccRCC_annot[,1],ccRCC_annot[,2],sep="")

sum(sample_annot$ccRCC%in%ccRCC_annot$Sample_ID)
ccRCC_annot <- ccRCC_annot[ccRCC_annot$Sample_ID%in%sample_annot$ccRCC,]

# melanoma
melanoma_annot <- read.xlsx(paste(output_dir,"Klinische data location III CRC, ccRCC, melanoma- compleet 17022022_edited.xlsx",sep=""),
                         3,startRow=1,colNames=T)
melanoma_annot[6,2] <- "P32" 
melanoma_annot$Sample_ID <- melanoma_annot[,2]
melanoma_annot$run <- "run231"

melanoma_annot_2 <- read.xlsx(paste(output_dir,"Clinical follow-up data additional melanoma samples - Location-III.xlsx",sep=""),1,colNames=T)
melanoma_annot_2[,3] <- tolower(melanoma_annot_2[,3])
melanoma_annot_2 <- cbind(melanoma_annot_2[,1:8],rep(NA,nrow(melanoma_annot_2)),melanoma_annot_2[,9:ncol(melanoma_annot_2)],melanoma_annot_2[,2])
melanoma_annot_2$Date.last.contact[is.na(melanoma_annot_2$Date.last.contact)] <- "?"
melanoma_annot_2$run <- "run239"
colnames(melanoma_annot_2) <- colnames(melanoma_annot)

melanoma_annot <- rbind(melanoma_annot,melanoma_annot_2)

sum(sample_annot$melanoma%in%melanoma_annot$Sample_ID)
melanoma_annot <- melanoma_annot[melanoma_annot$Sample_ID%in%sample_annot$melanoma,]

# make sure only annotated samples are included
sample_annot$CRC[!(sample_annot$CRC%in%CRC_annot$Sample_ID)] <- ""
sum(sample_annot$CRC%in%CRC_annot$Sample_ID)
length(CRC_annot$Sample_ID)

sample_annot$ccRCC[!(sample_annot$ccRCC%in%ccRCC_annot$Sample_ID)] <- ""
sum(sample_annot$ccRCC%in%ccRCC_annot$Sample_ID)
length(ccRCC_annot$Sample_ID)

sample_annot$melanoma[!(sample_annot$melanoma%in%melanoma_annot$Sample_ID)] <- ""
sum(sample_annot$melanoma%in%melanoma_annot$Sample_ID)
length(melanoma_annot$Sample_ID)

write.table(sample_annot,file="sample_annot.txt",col.names = T,row.names = F,quote=F,sep="\t")

# Load RRBS amplicon aware
amplicons <- c(primer_annot_LY75$name,primer_annot_GREM1$name,primer_annot_NDRG4$name)
amplicon_allrrbs_rows <- NULL
for (i in 1:length(amplicons)){
  print(amplicons[i])
  files_amplicon_i <- files[grepl(paste0(amplicons[i],"_"),files)]
  
  if (length(files_amplicon_i)>0){
    sample_amplicon_i <- paste0(gsub("_.+.cov$","",gsub(".+_run_","run",files_amplicon_i)),rep("_",length(files_amplicon_i)),
                                gsub("\\.run.+cov$","",gsub(".+_run_[0-9]+_[0-9]+_S","",files_amplicon_i)))
      
    # sort files in order of sample annot
    files_amplicon_i <- files_amplicon_i[order(sample_amplicon_i)]
    sample_amplicon_i <- sample_amplicon_i[order(sample_amplicon_i)]
    
    sample_annot_amplicon_i <- sample_annot[sample_annot$Sample_ID%in%sample_amplicon_i,]
    
    # Check if annotation has right length
    if (length(files_amplicon_i)!=length(sample_annot_amplicon_i[,1])){
      stop("colData has wrong dimensions")
    }

    # Check if annotation is in right order
    if (sum(sample_amplicon_i==sample_annot_amplicon_i[,1])!=length(sample_amplicon_i)){
      stop("colData in wrong order")
    }
    
    rrbs_amplicon_i <- readBismark(files_amplicon_i,sample_annot_amplicon_i[,1])
    
    totalReads_amplicon_i <- matrix(as.integer(rep(0,nrow(totalReads(rrbs_amplicon_i))*nrow(sample_annot))),
                                    ncol=nrow(sample_annot))
    colnames(totalReads_amplicon_i) <- sample_annot$Sample_ID
    methReads_amplicon_i <- totalReads_amplicon_i
    
    totalReads_amplicon_i[,colnames(totalReads_amplicon_i)%in%colnames(totalReads(rrbs_amplicon_i))] <- totalReads(rrbs_amplicon_i)
    methReads_amplicon_i[,colnames(methReads_amplicon_i)%in%colnames(methReads(rrbs_amplicon_i))] <- methReads(rrbs_amplicon_i)
    
    rrbs_amplicon_i <- BSraw(metadata = list(), rowRanges = rowRanges(rrbs_amplicon_i),
                             methReads_amplicon_i, totalReads_amplicon_i)
    
    if (i == 1){
      allrrbs <- rrbs_amplicon_i
    } else {
      allrrbs <- c(allrrbs,rrbs_amplicon_i)
    }
    
    amplicon_allrrbs_rows <- c(amplicon_allrrbs_rows,nrow(rrbs_amplicon_i))
  } else {
    amplicon_allrrbs_rows <- c(amplicon_allrrbs_rows,0)
  }
}

col_mask_CRC <- sample_annot$CRC!=""
col_mask_ccRCC <- sample_annot$ccRCC!=""
col_mask_melanoma <- sample_annot$melanoma!=""

amplicons_allrrbs <- rep(amplicons,amplicon_allrrbs_rows)

# Change names Manon & FPassay
primer_annot[primer_annot[,1]=="GREM1_Manon",1] <- "GREM1_Original"
primer_annot[primer_annot[,1]=="NDRG4_FPassay",1] <- "NDRG4_Original"

primer_annot_GREM1[primer_annot_GREM1[,1]=="GREM1_Manon",1] <- "GREM1_Original"
primer_annot_NDRG4[primer_annot_NDRG4[,1]=="NDRG4_FPassay",1] <- "NDRG4_Original"

amplicons_allrrbs[amplicons_allrrbs=="GREM1_Manon"] <- "GREM1_Original"
amplicons_allrrbs[amplicons_allrrbs=="NDRG4_FPassay"] <- "NDRG4_Original"

amplicons[amplicons=="GREM1_Manon"]  <- "GREM1_Original"
amplicons[amplicons=="NDRG4_FPassay"]  <- "NDRG4_Original"

###########################
# Filter qualitative loci #
###########################

rownames(allrrbs) <- amplicons_allrrbs

locs_allrrbs <- rowRanges(allrrbs)
locs_allrrbs$name <- amplicons_allrrbs

row_NDRG4 <- data.frame(locs_allrrbs)[,1]==16
row_GREM1 <- data.frame(locs_allrrbs)[,1]==15
row_LY75 <- data.frame(locs_allrrbs)[,1]==2

# Remove NDRG4_92 due to design error
totalReads(allrrbs)[grepl("NDRG4_92$",rownames(allrrbs)),] <- 1

# coverage <10 considered missing value
totalReads(allrrbs)[row_NDRG4,col_mask_CRC][totalReads(allrrbs)[row_NDRG4,col_mask_CRC]<10] <- NA
totalReads(allrrbs)[row_GREM1,col_mask_ccRCC][totalReads(allrrbs)[row_GREM1,col_mask_ccRCC]<10] <- NA
totalReads(allrrbs)[row_LY75,col_mask_melanoma][totalReads(allrrbs)[row_LY75,col_mask_melanoma]<10] <- NA

# coverage <10 considered missing value
perc_non_missing_val <- 0.9
mask_filter <- rep(FALSE,nrow(allrrbs))
mask_filter[row_NDRG4] <- rowSums(!is.na(totalReads(allrrbs)[row_NDRG4,col_mask_CRC]))>perc_non_missing_val*sum(col_mask_CRC) 
mask_filter[row_GREM1] <- rowSums(!is.na(totalReads(allrrbs)[row_GREM1,col_mask_ccRCC]))>perc_non_missing_val*sum(col_mask_ccRCC)  
mask_filter[row_LY75] <- rowSums(!is.na(totalReads(allrrbs)[row_LY75,col_mask_melanoma]))>perc_non_missing_val*sum(col_mask_melanoma)
sum(mask_filter)

# coverage <10 considered missing value
rrbs <- allrrbs[mask_filter]
amplicons_rrbs <- amplicons_allrrbs[mask_filter]

locs <- rowRanges(rrbs)
locs$name <- amplicons_rrbs

table(data.frame(locs)[,1])
# CRC: NDRG4 (chr16:58462846-58513628)
head(data.frame(locs)[data.frame(locs)[,1]==16,])
# ccRCC: GREM1 (chr15:32718004-32745106)
head(data.frame(locs)[data.frame(locs)[,1]==15,])
# melanoma: LY75 (GRch38: GRch37: chr2:15980335-159904756)
head(data.frame(locs)[data.frame(locs)[,1]==2,])

bs_NDRG4 <- rrbs[data.frame(locs)[,1]==16,]
locs_NDRG4 <- locs[data.frame(locs)[,1]==16,]
bs_LY75 <- rrbs[data.frame(locs)[,1]==2,]
locs_LY75 <- locs[data.frame(locs)[,1]==2,]
bs_GREM1 <- rrbs[data.frame(locs)[,1]==15,]
locs_GREM1 <- locs[data.frame(locs)[,1]==15,]

bs_NDRG4.rel <- rawToRel(bs_NDRG4)
bs_LY75.rel <- rawToRel(bs_LY75)
bs_GREM1.rel <- rawToRel(bs_GREM1)

length(unique(paste(data.frame(locs_NDRG4)[,1],data.frame(locs_NDRG4)[,2],sep="_")))
length(unique(paste(data.frame(locs_GREM1)[,1],data.frame(locs_GREM1)[,2],sep="_")))
length(unique(paste(data.frame(locs_LY75)[,1],data.frame(locs_LY75)[,2],sep="_")))

locs_NDRG4$av_cov <- rowMeans(totalReads(bs_NDRG4)[,col_mask_CRC],na.rm=T)
locs_GREM1$av_cov <- rowMeans(totalReads(bs_GREM1)[,col_mask_ccRCC],na.rm=T)
locs_LY75$av_cov <- rowMeans(totalReads(bs_LY75)[,col_mask_melanoma],na.rm=T)

# rowSums((totalReads(bs_NDRG4)[,col_mask_CRC]<10) & !is.na(totalReads(bs_NDRG4)[,col_mask_CRC]))
# colSums((totalReads(bs_NDRG4)[,col_mask_CRC]<10) & !is.na(totalReads(bs_NDRG4)[,col_mask_CRC]))
# rowSums((totalReads(bs_GREM1)[,col_mask_ccRCC]<10) & !is.na(totalReads(bs_GREM1)[,col_mask_ccRCC]))
# colSums((totalReads(bs_GREM1)[,col_mask_ccRCC]<10) & !is.na(totalReads(bs_GREM1)[,col_mask_ccRCC]))
# rowSums((totalReads(bs_LY75)[,col_mask_melanoma]<10) & !is.na(totalReads(bs_LY75)[,col_mask_melanoma]))
# colSums((totalReads(bs_LY75)[,col_mask_melanoma]<10) & !is.na(totalReads(bs_LY75)[,col_mask_melanoma]))

# sum((totalReads(bs_NDRG4)[,col_mask_CRC]<10) & !is.na(totalReads(bs_NDRG4)[,col_mask_CRC]))
# sum((totalReads(bs_NDRG4)[,col_mask_CRC]<10) & !is.na(totalReads(bs_NDRG4)[,col_mask_CRC]))
# sum((totalReads(bs_GREM1)[,col_mask_ccRCC]<10) & !is.na(totalReads(bs_GREM1)[,col_mask_ccRCC]))
# sum((totalReads(bs_GREM1)[,col_mask_ccRCC]<10) & !is.na(totalReads(bs_GREM1)[,col_mask_ccRCC]))
# sum((totalReads(bs_LY75)[,col_mask_melanoma]<10) & !is.na(totalReads(bs_LY75)[,col_mask_melanoma]))
# sum((totalReads(bs_LY75)[,col_mask_melanoma]<10) & !is.na(totalReads(bs_LY75)[,col_mask_melanoma]))

# Number of CpGs
length(locs_GREM1)
length(locs_NDRG4)
length(locs_LY75)

# Unique CpGs
length(unique(paste(data.frame(locs_GREM1)[,1],data.frame(locs_GREM1)[,2],sep="_")))
length(unique(paste(data.frame(locs_NDRG4)[,1],data.frame(locs_NDRG4)[,2],sep="_")))
length(unique(paste(data.frame(locs_LY75)[,1],data.frame(locs_LY75)[,2],sep="_")))

# Average coverage
summary(locs_GREM1$av_cov)
summary(locs_NDRG4$av_cov)
summary(locs_LY75$av_cov)

########################
# Output betas to file #
########################

# Set working directory
setwd(output_dir)

rrbs.rel <- rawToRel(rrbs)
dat_meth_rates <- cbind(data.frame(locs),methLevel(rrbs.rel))

colnames(dat_meth_rates)[1] <- c("chr")
colnames(dat_meth_rates)[6] <- c("amplicon")

write.table(dat_meth_rates,file="methylation_rates.txt",row.names = F,col.names = T,sep="\t",quote = F)

#############################
# Get transcript annotation #
#############################

# Import gene structure 
if (!file.exists("Genes_struct.Rda")){
  ensembl <- useEnsembl(biomart="ensembl",GRCh=37, dataset="hsapiens_gene_ensembl")
  genes_struct <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id','ensembl_transcript_id','ensembl_exon_id','chromosome_name',
                                       'start_position','end_position',"exon_chrom_start","exon_chrom_end",'transcription_start_site')
                        ,mart = ensembl,filter="external_gene_name",values=c("NDRG4","LY75","GREM1"))
  save(genes_struct,file="Genes_struct.Rda")
  print(head(genes_struct))
} else {
  load("Genes_struct.Rda")
  print(head(genes_struct))
}

######################################
######################################
######################################
###                                ###
###                                ###   
###   Inter amplicon differences   ###
###                                ###
###                                ###
######################################
######################################
######################################


######################################
## Analysis CRC (NDRG4): DIAGNOSTIC ##
######################################

betas_NDRG4 <- methLevel(bs_NDRG4.rel)[,sample_annot$CRC!=""]
rownames(betas_NDRG4) <- paste(data.frame(locs_NDRG4)[,1],data.frame(locs_NDRG4)[,2],data.frame(locs_NDRG4)[,6],sep="_")
colnames(betas_NDRG4) <- sample_annot$CRC[sample_annot$CRC!=""]
print(head(betas_NDRG4))

dat_NDRG4 <- cbind(data.frame(locs_NDRG4)[,c(1:2,6)],betas_NDRG4)
colnames(dat_NDRG4)[1:3] <- c("CHR","POS","AMPLICON")
write.table(dat_NDRG4,file=paste(output_dir,"Meth_NDRG4.txt",sep=""),row.names = F,col.names = T,sep="\t",quote = F)

CRC_annot_sort <- CRC_annot[sort(CRC_annot$Sample_ID,index.return=T)$ix,]
CRC_annot_sort <- CRC_annot_sort[CRC_annot_sort$Sample_ID%in%colnames(betas_NDRG4),]
betas_NDRG4_sort <- betas_NDRG4[,sort(colnames(betas_NDRG4),index.return=T)$ix]
betas_NDRG4_sort <- betas_NDRG4_sort[,colnames(betas_NDRG4_sort)%in%CRC_annot_sort$Sample_ID]
sum(colnames(betas_NDRG4_sort)==CRC_annot_sort$Sample_ID)

# Correction for age and sex
colnames(CRC_annot_sort) <- c("Sample","Tissue","Gender","Age","Meth","Sample_ID")
CRC_annot_sort$Gender <- as.factor(CRC_annot_sort$Gender)
CRC_annot_sort$Meth <- as.factor(CRC_annot_sort$Meth)
CRC_annot_sort$Tissue <- as.factor(CRC_annot_sort$Tissue)
print(CRC_annot_sort$Gender)
print(CRC_annot_sort$Age)
print(CRC_annot_sort$Tissue)
print(CRC_annot_sort$Meth)
CRC_annot_sort$Patient_ID <- gsub("[TN]$","",CRC_annot_sort$Sample_ID)
write.table(CRC_annot_sort,file="CRC_annot_sort.txt",col.names = T,row.names = F,quote = F,sep="\t")

# Mvalues
Mval_NDRG4_sort <- log2((betas_NDRG4_sort+ct)/(1-betas_NDRG4_sort+ct))

# ROC and logit per CpG
dat_NDRG4_ROC_logit_perCpG <- data.frame(ID=rownames(Mval_NDRG4_sort),
                                         Chr=as.numeric(gsub("_.+$","",rownames(Mval_NDRG4_sort))),
                                         Pos=as.numeric(gsub("_.+$","",gsub("^[0-9]+_","",rownames(Mval_NDRG4_sort)))),
                                         Amplicon=as.character(gsub("^[0-9]+_[0-9]+_","",rownames(Mval_NDRG4_sort))),
                                         stringsAsFactors = F)
CRC_annot_sort_BS <- CRC_annot_sort

# cox regression confounders
formula_CF <- formula(paste("Tissue ~ Gender+Age+Patient_ID",sep=""))
cn_formula_CF <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_CF))),"+",fixed=TRUE)))
model_CF <- glm(formula_CF, 
                data = CRC_annot_sort_BS, 
                family = binomial(link="logit"))
print(summary(model_CF))

# cox regression msp
formula_MSP <- formula(paste("Tissue ~ Meth+Patient_ID",sep=""))
cn_formula_MSP <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_MSP))),"+",fixed=TRUE)))
model_MSP <- glm(formula_MSP, 
                data = CRC_annot_sort_BS, 
                family = binomial(link="logit"))
print(summary(model_MSP))

print(sprintf("AIC for baseline model: %s",AIC(model_CF)))
print(sprintf("Delta AIC for MSP amplicon: %s",AIC(model_MSP)-AIC(model_CF)))

dat_NDRG4_ROC_logit_perCpG$logit <- NA
dat_NDRG4_ROC_logit_perCpG$AIC <- NA
dat_NDRG4_ROC_logit_perCpG$AIC_CF <- AIC(model_CF)
dat_NDRG4_ROC_logit_perCpG$est_logit <- NA
dat_NDRG4_ROC_logit_perCpG$delta_AIC <- NA
for (i in 1:nrow(Mval_NDRG4_sort)){
  CRC_annot_sort_BS[paste("Meth_CpG_",i,sep="")] <- Mval_NDRG4_sort[i,]
  
  # logit 
  formula_i <- formula(paste("Tissue~Meth_CpG_",i,"+Patient_ID",sep=""))
  cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
  model_BS_i <- glm(formula_i, 
                    data = CRC_annot_sort_BS, 
                    family = binomial(link="logit"))
  dat_NDRG4_ROC_logit_perCpG$logit[i] <- -log10(summary(model_BS_i)$coefficients[2,4])
  dat_NDRG4_ROC_logit_perCpG$est_logit[i] <- summary(model_BS_i)$coefficients[2,1]
  
  print(summary(model_BS_i))
  
  # ROC BS
  dat_NDRG4_ROC_logit_perCpG$AIC[i] <- AIC(model_BS_i)
  dat_NDRG4_ROC_logit_perCpG$delta_AIC[i] <- dat_NDRG4_ROC_logit_perCpG$AIC[i] - dat_NDRG4_ROC_logit_perCpG$AIC_CF[i]
  
}                 
write.table(dat_NDRG4_ROC_logit_perCpG,"NDRG4_AIC_table.txt",col.names = T,row.names = F,quote = F, sep="\t")

########################################
# Analysis melanoma (LY75): PROGNOSTIC #
########################################

betas_LY75 <- methLevel(bs_LY75.rel)[,sample_annot$melanoma!=""]
rownames(betas_LY75) <- paste(data.frame(locs_LY75)[,1],data.frame(locs_LY75)[,2],data.frame(locs_LY75)[,6],sep="_")
colnames(betas_LY75) <- sample_annot$melanoma[sample_annot$melanoma!=""]
print(head(betas_LY75))

dat_LY75 <- cbind(data.frame(locs_LY75)[,c(1:2,6)],betas_LY75)
colnames(dat_LY75)[1:3] <- c("CHR","POS","AMPLICON")
write.table(dat_LY75,file=paste(output_dir,"Meth_LY75.txt",sep=""),row.names = F,col.names = T,sep="\t",quote = F)

melanoma_annot <- melanoma_annot[!(melanoma_annot$Date.last.contact=="?"),]

melanoma_annot_sort <- melanoma_annot[sort(melanoma_annot$Sample_ID,index.return=T)$ix,]
melanoma_annot_sort <- melanoma_annot_sort[melanoma_annot_sort$Sample_ID%in%colnames(betas_LY75),]
betas_LY75_sort <- betas_LY75[,sort(colnames(betas_LY75),index.return=T)$ix]
betas_LY75_sort <- betas_LY75_sort[,colnames(betas_LY75_sort)%in%melanoma_annot_sort$Sample_ID]
melanoma_annot_sort$NK_death <- melanoma_annot_sort$MK_death
melanoma_annot_sort$NK_death <- as.numeric(melanoma_annot_sort$NK_death)
sum(colnames(betas_LY75_sort)==melanoma_annot_sort$Sample_ID)
write.table(melanoma_annot_sort,file="melanoma_annot_sort.txt",col.names = T,row.names = F,quote = F,sep="\t")

# Correctie in model voor leeftijd, geslacht en Breslow
colnames(melanoma_annot_sort) <- c("Study_Nr","Sample","Gender","Age","Location","Date_diag",
                                   "Breslow","Ulceration","T","Date_last_followup","Outcome","Meth","Sample_ID","Run")
melanoma_annot_sort$Gender <- toupper(melanoma_annot_sort$Gender)
melanoma_annot_sort$Gender <- as.factor(melanoma_annot_sort$Gender)
melanoma_annot_sort$Run <- as.factor(melanoma_annot_sort$Run)
melanoma_annot_sort$Meth <- as.factor(melanoma_annot_sort$Meth)
melanoma_annot_sort$Breslow <- as.numeric(melanoma_annot_sort$Breslow)
melanoma_annot_sort$Date_diag <- as.Date(melanoma_annot_sort$Date_diag,origin="1899-12-30")
melanoma_annot_sort$Date_last_followup <- as.numeric(melanoma_annot_sort$Date_last_followup)
melanoma_annot_sort$Date_last_followup <- as.Date(melanoma_annot_sort$Date_last_followup,origin="1899-12-30")
melanoma_annot_sort$Followup_time <- melanoma_annot_sort$Date_last_followup - melanoma_annot_sort$Date_diag
melanoma_annot_sort$MK_death <- melanoma_annot_sort$Outcome=="DOD"

betas_LY75_sort <- betas_LY75_sort[,melanoma_annot_sort$Followup_time>0]
melanoma_annot_sort <- melanoma_annot_sort[melanoma_annot_sort$Followup_time>0,]

print(melanoma_annot_sort$Gender)
print(melanoma_annot_sort$Age)
print(melanoma_annot_sort$Meth)

# Mvalues
Mval_LY75_sort <- log2((betas_LY75_sort+ct)/(1-betas_LY75_sort+ct))

# ROC and logit per CpG
dat_LY75_ROC_logit_perCpG <- data.frame(ID=rownames(Mval_LY75_sort),
                                         Chr=as.numeric(gsub("_.+$","",rownames(Mval_LY75_sort))),
                                         Pos=as.numeric(gsub("_.+$","",gsub("^[0-9]+_","",rownames(Mval_LY75_sort)))),
                                         Amplicon=as.character(gsub("^[0-9]+_[0-9]+_","",rownames(Mval_LY75_sort))),
                                         stringsAsFactors = F)
melanoma_annot_sort_BS <- melanoma_annot_sort

# cox regression confounders
# formula_CF <- formula(paste("Surv(Followup_time, MK_death) ~ Gender+Age+Breslow",sep=""))
# formula_CF <- formula(paste("Surv(Followup_time, MK_death) ~ Gender+Age+Run",sep=""))
formula_CF <- formula(paste("Surv(Followup_time, MK_death) ~ Gender+Age",sep=""))
cn_formula_CF <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_CF))),"+",fixed=TRUE)))
cox_CF <- coxph(formula_CF, data = melanoma_annot_sort_BS)
print(summary(cox_CF))

# cox regression msp
# formula_MSP <- formula(paste("Surv(Followup_time, MK_death) ~ Meth+Run",sep=""))
formula_MSP <- formula(paste("Surv(Followup_time, MK_death) ~ Meth",sep=""))
cn_formula_MSP <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_MSP))),"+",fixed=TRUE)))
cox_MSP <- coxph(formula_MSP, data = melanoma_annot_sort_BS)
print(summary(cox_MSP))

print(sprintf("AIC for baseline model: %s",AIC(cox_CF)))
print(sprintf("Delta AIC for MSP amplicon: %s",AIC(cox_MSP)-AIC(cox_CF)))

dat_LY75_ROC_logit_perCpG$logit <- NA
dat_LY75_ROC_logit_perCpG$AIC <- NA
dat_LY75_ROC_logit_perCpG$AIC_CF <- AIC(cox_CF)
dat_LY75_ROC_logit_perCpG$est_logit <- NA
dat_LY75_ROC_logit_perCpG$delta_AIC <- NA
for (i in 1:nrow(Mval_LY75_sort)){
  melanoma_annot_sort_BS[paste("Meth_CpG_",i,sep="")] <- Mval_LY75_sort[i,]
  
  # cox regression
  # formula_i <- formula(paste("Surv(Followup_time, MK_death) ~ Meth_CpG_",i,"+Gender+Age+Breslow",sep=""))
  # formula_i <- formula(paste("Surv(Followup_time, MK_death) ~ Meth_CpG_",i,"+Run",sep=""))
  formula_i <- formula(paste("Surv(Followup_time, MK_death) ~ Meth_CpG_",i,sep=""))
  cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
  cox_BS_i <- coxph(formula_i, data = melanoma_annot_sort_BS)
  print(summary(cox_BS_i))
  dat_LY75_ROC_logit_perCpG$logit[i] <- -log10(summary(cox_BS_i)$coefficients[1,5])
  dat_LY75_ROC_logit_perCpG$est_logit[i] <- summary(cox_BS_i)$coefficients[1,1]
  
  # ROC BS
  dat_LY75_ROC_logit_perCpG$AIC[i] <- AIC(cox_BS_i)
  dat_LY75_ROC_logit_perCpG$delta_AIC[i] <- dat_LY75_ROC_logit_perCpG$AIC[i] - dat_LY75_ROC_logit_perCpG$AIC_CF[i]
  
}                               
write.table(dat_LY75_ROC_logit_perCpG,"LY75_AIC_table.txt",col.names = T,row.names = F,quote = F, sep="\t")

######################################
# Analysis ccRCC (GREM1): PROGNOSTIC #
######################################

betas_GREM1 <- methLevel(bs_GREM1.rel)[,sample_annot$ccRCC!=""]
betas_GREM1 <- betas_GREM1[sort(data.frame(locs_GREM1)[,2],index.return=T)$ix,]
locs_GREM1 <- locs_GREM1[sort(data.frame(locs_GREM1)[,2],index.return=T)$ix,]
rownames(betas_GREM1) <- paste(data.frame(locs_GREM1)[,1],data.frame(locs_GREM1)[,2],data.frame(locs_GREM1)[,6],sep="_")
colnames(betas_GREM1) <- sample_annot$ccRCC[sample_annot$ccRCC!=""]
print(head(betas_GREM1))

dat_GREM1 <- cbind(data.frame(locs_GREM1)[,c(1:2,6)],betas_GREM1)
colnames(dat_GREM1)[1:3] <- c("CHR","POS","AMPLICON")
write.table(dat_GREM1,file=paste(output_dir,"Meth_GREM1.txt",sep=""),row.names = F,col.names = T,sep="\t",quote = F)

ccRCC_annot_sort <- ccRCC_annot[sort(ccRCC_annot$Sample_ID,index.return=T)$ix,]
ccRCC_annot_sort <- ccRCC_annot_sort[ccRCC_annot_sort$Sample_ID%in%colnames(betas_GREM1),]
betas_GREM1_sort <- betas_GREM1[,sort(colnames(betas_GREM1),index.return=T)$ix]
betas_GREM1_sort <- betas_GREM1_sort[,colnames(betas_GREM1_sort)%in%ccRCC_annot_sort$Sample_ID]
sum(colnames(betas_GREM1_sort)==ccRCC_annot_sort$Sample_ID)
write.table(ccRCC_annot_sort,file="ccRCC_annot_sort.txt",col.names = T,row.names = F,quote = F,sep="\t")


# Correctie in model voor leeftijd, geslacht, Stage, Size
colnames(ccRCC_annot_sort) <- c("Sample","Tissue","Prognosis","Gender","Age","Date_diag","Grade_old","Size","Stage",
                                   "Date_last_followup","NK_death","Meth","Grade_new","Lymph_invasion","Necrosis","Sample_ID")
ccRCC_annot_sort$Gender <- as.factor(ccRCC_annot_sort$Gender)
ccRCC_annot_sort$Meth <- as.factor(ccRCC_annot_sort$Meth)
ccRCC_annot_sort$Date_diag <- as.Date(ccRCC_annot_sort$Date_diag,origin="1899-12-30")
ccRCC_annot_sort$Date_last_followup <- as.Date(ccRCC_annot_sort$Date_last_followup,origin="1899-12-30")
ccRCC_annot_sort$Followup_time <- ccRCC_annot_sort$Date_last_followup - ccRCC_annot_sort$Date_diag
print(ccRCC_annot_sort$Gender)
print(ccRCC_annot_sort$Age)
print(ccRCC_annot_sort$Meth)


# Mvalues
Mval_GREM1_sort <- log2((betas_GREM1_sort+ct)/(1-betas_GREM1_sort+ct))

# ROC and logit per CpG
dat_GREM1_ROC_logit_perCpG <- data.frame(ID=rownames(Mval_GREM1_sort),
                                        Chr=as.numeric(gsub("_.+$","",rownames(Mval_GREM1_sort))),
                                        Pos=as.numeric(gsub("_.+$","",gsub("^[0-9]+_","",rownames(Mval_GREM1_sort)))),
                                        Amplicon=as.character(gsub("^[0-9]+_[0-9]+_","",rownames(Mval_GREM1_sort))),
                                        stringsAsFactors = F)
ccRCC_annot_sort_BS <- ccRCC_annot_sort

# cox regression confounders
# formula_CF <- formula(paste("Surv(Followup_time, NK_death) ~ Gender+Age+Stage+Size",sep=""))
formula_CF <- formula(paste("Surv(Followup_time, NK_death) ~ Gender+Age",sep=""))
cn_formula_CF <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_CF))),"+",fixed=TRUE)))
cox_CF <- coxph(formula_CF, data = ccRCC_annot_sort_BS)
print(summary(cox_CF))

# cox regression msp
formula_MSP <- formula(paste("Surv(Followup_time, NK_death) ~ Meth",sep=""))
cn_formula_MSP <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_MSP))),"+",fixed=TRUE)))
cox_MSP <- coxph(formula_MSP, data = ccRCC_annot_sort_BS)
print(summary(cox_MSP))

print(sprintf("AIC for baseline model: %s",AIC(cox_CF)))
print(sprintf("Delta AIC for MSP amplicon: %s",AIC(cox_MSP)-AIC(cox_CF)))

dat_GREM1_ROC_logit_perCpG$logit <- NA
dat_GREM1_ROC_logit_perCpG$AIC <- NA
dat_GREM1_ROC_logit_perCpG$AIC_CF <- AIC(cox_CF)
dat_GREM1_ROC_logit_perCpG$est_logit <- NA
dat_GREM1_ROC_logit_perCpG$delta_AIC <- NA
for (i in 1:nrow(Mval_GREM1_sort)){
  ccRCC_annot_sort_BS[paste("Meth_CpG_",i,sep="")] <- Mval_GREM1_sort[i,]
  
  # cox regression incl Meth
  # formula_i <- formula(paste("Surv(Followup_time, NK_death) ~ Meth_CpG_",i,"+Gender+Age+Stage+Size",sep=""))
  formula_i <- formula(paste("Surv(Followup_time, NK_death) ~ Meth_CpG_",i,sep=""))
  cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
  cox_BS_i <- coxph(formula_i, data = ccRCC_annot_sort_BS)
  print(summary(cox_BS_i))
  dat_GREM1_ROC_logit_perCpG$logit[i] <- -log10(summary(cox_BS_i)$coefficients[1,5])
  dat_GREM1_ROC_logit_perCpG$est_logit[i] <- summary(cox_BS_i)$coefficients[1,1]

  # ROC BS
  dat_GREM1_ROC_logit_perCpG$AIC[i] <- AIC(cox_BS_i)
  dat_GREM1_ROC_logit_perCpG$delta_AIC[i] <- dat_GREM1_ROC_logit_perCpG$AIC[i] - dat_GREM1_ROC_logit_perCpG$AIC_CF[i]
}              
write.table(dat_GREM1_ROC_logit_perCpG,"GREM1_AIC_table.txt",col.names = T,row.names = F,quote = F, sep="\t")

####################
# Overview objects #
####################

nrow(dat_NDRG4_ROC_logit_perCpG)
length(unique(paste(dat_NDRG4_ROC_logit_perCpG$Chr,dat_NDRG4_ROC_logit_perCpG$Pos,sep="_")))
nrow(dat_GREM1_ROC_logit_perCpG)
length(unique(paste(dat_GREM1_ROC_logit_perCpG$Chr,dat_GREM1_ROC_logit_perCpG$Pos,sep="_")))
nrow(dat_LY75_ROC_logit_perCpG)
length(unique(paste(dat_LY75_ROC_logit_perCpG$Chr,dat_LY75_ROC_logit_perCpG$Pos,sep="_")))

head(dat_NDRG4_ROC_logit_perCpG)
head(dat_GREM1_ROC_logit_perCpG)
head(dat_LY75_ROC_logit_perCpG)

#####################
# Primer_annot_list #
#####################

primer_annot$len_fwd <- nchar(primer_annot$fwd_primer)
primer_annot$len_rev <- nchar(primer_annot$rev_primer)
primer_annot$CpG_fwd <- grepl("[CY][GR]",primer_annot$fwd_primer)
primer_annot$CpG_rev <- grepl("[CY][GR]",primer_annot$rev_primer)

primer_annot_list <- list(primer_annot[grepl("GREM1_",primer_annot$name),],
                          primer_annot[grepl("NDRG4_",primer_annot$name),],
                          primer_annot[grepl("LY75_",primer_annot$name),])
names(primer_annot_list) <- c("GREM1","NDRG4","LY75")

for (i in names(primer_annot_list)){
  for (j in 4:9){
    primer_annot_list[[i]][,j] <- as.numeric(primer_annot_list[[i]][,j])
  }
}

for (gene in names(primer_annot_list)){
  primer_annot_list[[gene]]$meth_CpG1_fwd <- NA
  primer_annot_list[[gene]]$meth_CpG1_rev <- NA
  primer_annot_list[[gene]]$meth_CpG2_fwd <- NA
  primer_annot_list[[gene]]$meth_CpG2_rev <- NA
  primer_annot_list[[gene]]$meth_CpG3_fwd <- NA
  primer_annot_list[[gene]]$meth_CpG3_rev <- NA
  locs_g <- data.frame(get(paste0("locs_",gene)))
  betas_g <-  get(paste0("betas_",gene))
  for (i in 1:nrow(primer_annot_list[[gene]])){
    mask_CpG_in_amplicon_fwd <- (locs_g$start>primer_annot_list[[gene]]$begin[i]) & 
      (locs_g$start<=(primer_annot_list[[gene]]$begin[i]+primer_annot_list[[gene]]$len_fwd[i]-1)) &
      (locs_g$name==primer_annot_list[[gene]]$name[i])
    mask_CpG_in_amplicon_rev <- (locs_g$start<primer_annot_list[[gene]]$end[i]) & 
      (locs_g$start>=(primer_annot_list[[gene]]$end[i]-primer_annot_list[[gene]]$len_rev[i]+1)) &
      (locs_g$name==primer_annot_list[[gene]]$name[i])

    if (sum(mask_CpG_in_amplicon_fwd)>0){primer_annot_list[[gene]]$CpG_fwd[i] <- TRUE}
    if (sum(mask_CpG_in_amplicon_rev)>0){primer_annot_list[[gene]]$CpG_rev[i] <- TRUE}

    # Sense
    if(sum(mask_CpG_in_amplicon_fwd)==1){
      print(paste0(gene," (",i,", ",primer_annot_list[[gene]]$name[i],"); Sense: ",sum(mask_CpG_in_amplicon_fwd)))
      primer_annot_list[[gene]]$meth_CpG1_fwd[i] <- mean(betas_g[paste(locs_g$seqnames[mask_CpG_in_amplicon_fwd],
                                                       locs_g$start[mask_CpG_in_amplicon_fwd],
                                                       locs_g$name[mask_CpG_in_amplicon_fwd],sep="_"),],na.rm=T)
    } else if (sum(mask_CpG_in_amplicon_fwd)>1){
      print(paste0(gene," (",i,", ",primer_annot_list[[gene]]$name[i],"); Sense: ",sum(mask_CpG_in_amplicon_fwd)))
      for (j in 1:sum(mask_CpG_in_amplicon_fwd)){
        primer_annot_list[[gene]][i,colnames(primer_annot_list[[gene]])==paste0("meth_CpG",j,"_fwd")] <- mean(betas_g[paste(locs_g$seqnames[mask_CpG_in_amplicon_fwd],
                                                                                                            locs_g$start[mask_CpG_in_amplicon_fwd],
                                                                                                            locs_g$name[mask_CpG_in_amplicon_fwd],sep="_"),][j,],na.rm=T)
      }
    } else if ((sum(mask_CpG_in_amplicon_fwd)==0) & (primer_annot_list[[gene]]$CpG_fwd[i])){
      print(paste0(gene," (",i,", ",primer_annot_list[[gene]]$name[i],"); Sense: ",sum(mask_CpG_in_amplicon_fwd)))
    }

    # Antisense
    if(sum(mask_CpG_in_amplicon_rev)==1){
      print(paste0(gene," (",i,", ",primer_annot_list[[gene]]$name[i],"); Antisense: ",sum(mask_CpG_in_amplicon_rev)))
      primer_annot_list[[gene]]$meth_CpG1_rev[i] <- mean(betas_g[paste(locs_g$seqnames[mask_CpG_in_amplicon_rev],
                                                        locs_g$start[mask_CpG_in_amplicon_rev],
                                                        locs_g$name[mask_CpG_in_amplicon_rev],sep="_"),],na.rm=T)
    } else if (sum(mask_CpG_in_amplicon_rev)>1){
      print(paste0(gene," (",i,", ",primer_annot_list[[gene]]$name[i],"); Antisense: ",sum(mask_CpG_in_amplicon_rev)))
      for (j in sum(mask_CpG_in_amplicon_rev):1){
        primer_annot_list[[gene]][i,colnames(primer_annot_list[[gene]])==paste0("meth_CpG",j,"_rev")] <- mean(betas_g[paste(locs_g$seqnames[mask_CpG_in_amplicon_rev],
                                                                                                             locs_g$start[mask_CpG_in_amplicon_rev],
                                                                                                             locs_g$name[mask_CpG_in_amplicon_rev],sep="_"),][j,],na.rm=T)
      }
    } else if ((sum(mask_CpG_in_amplicon_rev)==0) & (primer_annot_list[[gene]]$CpG_rev[i])){
      print(paste0(gene," (",i,", ",primer_annot_list[[gene]]$name[i],"); Antisense: ",sum(mask_CpG_in_amplicon_rev)))
    }
  }
}

cbind(rbind(primer_annot_list[["GREM1"]],primer_annot_list[["NDRG4"]],primer_annot_list[["LY75"]])[,1],rowSums(!is.na(rbind(primer_annot_list[["GREM1"]],primer_annot_list[["NDRG4"]],primer_annot_list[["LY75"]])[,12:17])))

primer_annot_list[["GREM1"]]
primer_annot_list[["NDRG4"]]
primer_annot_list[["LY75"]]

write.table(rbind(primer_annot_list[["GREM1"]],primer_annot_list[["NDRG4"]],primer_annot_list[["LY75"]]),
            file="supplementary_table_1.txt",col.names=T,row.names=F,quote=F,sep="\t")

###########################
# Average delta AIC table #
###########################

avAIC_table <- NULL
for (k in 1:nrow(primer_annot_GREM1)){
  mask_amplicon <- grepl(paste0(primer_annot_GREM1$name[k],"$"),rownames(Mval_GREM1_sort))
  avAIC_table <- rbind(avAIC_table,
                       c(primer_annot_GREM1$name[k],sum(mask_amplicon),
                         mean(dat_GREM1_ROC_logit_perCpG$delta_AIC[mask_amplicon]),
                         min(dat_GREM1_ROC_logit_perCpG$delta_AIC[mask_amplicon]),
                         mean(totalReads(bs_GREM1)[grepl(primer_annot_GREM1$name[k],locs_GREM1$name),col_mask_ccRCC],na.rm=T)))
}
for (k in 1:nrow(primer_annot_LY75)){
  mask_amplicon <- grepl(paste0(primer_annot_LY75$name[k],"$"),rownames(Mval_LY75_sort))
  avAIC_table <- rbind(avAIC_table,
                       c(primer_annot_LY75$name[k],sum(mask_amplicon),
                         mean(dat_LY75_ROC_logit_perCpG$delta_AIC[mask_amplicon]),
                         min(dat_LY75_ROC_logit_perCpG$delta_AIC[mask_amplicon]),
                         mean(totalReads(bs_LY75)[grepl(primer_annot_LY75$name[k],locs_LY75$name),col_mask_melanoma],na.rm=T)))
}
for (k in 1:nrow(primer_annot_NDRG4)){
  mask_amplicon <- grepl(paste0(primer_annot_NDRG4$name[k],"$"),rownames(Mval_NDRG4_sort))
  avAIC_table <- rbind(avAIC_table,
                       c(primer_annot_NDRG4$name[k],sum(mask_amplicon),
                         mean(dat_NDRG4_ROC_logit_perCpG$delta_AIC[mask_amplicon]),
                         min(dat_NDRG4_ROC_logit_perCpG$delta_AIC[mask_amplicon]),
                         mean(totalReads(bs_NDRG4)[grepl(primer_annot_NDRG4$name[k],locs_NDRG4$name),col_mask_CRC],na.rm=T)))
}
colnames(avAIC_table) <- c("AmpliconID","Am_CpG","avAIC","minAIC","Am_cov")
avAIC_table <- data.frame(avAIC_table)
for (i in 2:4){
  avAIC_table[,i] <- as.numeric(avAIC_table[,i])
}

for (i in 1:nrow(avAIC_table)){
  gene_i <- gsub('_.*$','',avAIC_table$AmpliconID[i])
  avAIC_table$CpG_in_FWD[i] <- primer_annot_list[[gene_i]]$CpG_fwd[primer_annot_list[[gene_i]]$name==avAIC_table$AmpliconID[i]]
  avAIC_table$CpG_in_REV[i] <- primer_annot_list[[gene_i]]$CpG_rev[primer_annot_list[[gene_i]]$name==avAIC_table$AmpliconID[i]]
  avAIC_table$n_CpG_FWD[i] <- sum(!is.na(primer_annot_list[[gene_i]][primer_annot_list[[gene_i]]$name==avAIC_table$AmpliconID[i],grepl("CpG[1-3]_fwd$",colnames(primer_annot_list[[gene_i]]))]))
  avAIC_table$n_CpG_REV[i]  <- sum(!is.na(primer_annot_list[[gene_i]][primer_annot_list[[gene_i]]$name==avAIC_table$AmpliconID[i],grepl("CpG[1-3]_rev$",colnames(primer_annot_list[[gene_i]]))]))
}
          
write.table(avAIC_table,"avAIC_table.txt",col.names = T,row.names = F,quote = F,sep="\t")
rownames(avAIC_table) <- avAIC_table$AmpliconID

###############
# Data correl #
###############

dat_correl <- data.frame(ID=c(dat_GREM1_ROC_logit_perCpG$ID,dat_NDRG4_ROC_logit_perCpG$ID,dat_LY75_ROC_logit_perCpG$ID),
                         deltaAIC=c(dat_GREM1_ROC_logit_perCpG$delta_AIC,dat_NDRG4_ROC_logit_perCpG$delta_AIC,dat_LY75_ROC_logit_perCpG$delta_AIC),
                         CG_in_Am_FWD=avAIC_table[c(dat_GREM1_ROC_logit_perCpG$Amplicon,dat_NDRG4_ROC_logit_perCpG$Amplicon,dat_LY75_ROC_logit_perCpG$Amplicon),]$n_CpG_FWD,
                         CG_in_Am_REV=avAIC_table[c(dat_GREM1_ROC_logit_perCpG$Amplicon,dat_NDRG4_ROC_logit_perCpG$Amplicon,dat_LY75_ROC_logit_perCpG$Amplicon),]$n_CpG_REV)

dat_correl$CG_qual <- as.factor((dat_correl$CG_in_Am_FWD!=0) + (dat_correl$CG_in_Am_REV!=0))
dat_correl$CG_quant <- as.factor(dat_correl$CG_in_Am_FWD + dat_correl$CG_in_Am_REV)

dat_correl$dist_CpG <- NA
dat_correl$av_meth_primer_CpG <- NA
dat_correl$av_meth_CpG <- NA
for (i in 1:nrow(dat_correl)){
  chr_i <- unlist(strsplit(dat_correl$ID[i],"_"))[1]
  pos_i <- unlist(strsplit(dat_correl$ID[i],"_"))[2]
  gene_i <- unlist(strsplit(dat_correl$ID[i],"_"))[3]
  amplicon_i <- paste(unlist(strsplit(dat_correl$ID[i],"_"))[3:length(unlist(strsplit(dat_correl$ID[i],"_")))],collapse="_")

  primer_annot_i <- primer_annot_list[[gene_i]][primer_annot_list[[gene_i]]$name==amplicon_i,]

  dat_i <- get(paste0("dat_",gene_i,"_ROC_logit_perCpG"))
  dat_i <- dat_i[dat_i$Amplicon==amplicon_i,]

  betas_i <- rowMeans(get(paste0("betas_",gene_i)),na.rm=T)

  dat_correl$av_meth_CpG[i] <- betas_i[dat_correl$ID[i]]

  CpGs_in_primer_i <- dat_i$ID[(dat_i$Pos[dat_i$Amplicon==amplicon_i]>=primer_annot_i$begin) & (dat_i$Pos[dat_i$Amplicon==amplicon_i]<=(primer_annot_i$begin+primer_annot_i$len_fwd-1)) |
                                  (dat_i$Pos[dat_i$Amplicon==amplicon_i]>=(primer_annot_i$end-primer_annot_i$len_rev+1)) & (dat_i$Pos[dat_i$Amplicon==amplicon_i]<=primer_annot_i$end)]
  CpGs_in_primer_pos_i <- dat_i$Pos[(dat_i$Pos[dat_i$Amplicon==amplicon_i]>=primer_annot_i$begin) & (dat_i$Pos[dat_i$Amplicon==amplicon_i]<=(primer_annot_i$begin+primer_annot_i$len_fwd-1)) |
                                  (dat_i$Pos[dat_i$Amplicon==amplicon_i]>=(primer_annot_i$end-primer_annot_i$len_rev+1)) & (dat_i$Pos[dat_i$Amplicon==amplicon_i]<=primer_annot_i$end)]

  if (length(CpGs_in_primer_i)!=0){
    dat_correl$dist_CpG[i] <- min(abs(as.numeric(pos_i)-as.numeric(CpGs_in_primer_pos_i)))
    dat_correl$av_meth_primer_CpG[i] <- betas_i[CpGs_in_primer_i[which.min(abs(as.numeric(pos_i)-as.numeric(CpGs_in_primer_pos_i)))]]
  }
}
dat_correl$diff_meth <- abs(dat_correl$av_meth_CpG-0.5)
rownames(dat_correl) <- dat_correl$ID

# dat_correl[grepl("GREM1_53",rownames(dat_correl)),]
# dat_correl[grepl("NDRG4_Original",rownames(dat_correl)),]
# dat_correl[grepl("LY75_9_S",rownames(dat_correl)),]

################################
# Plot CRC (NDRG4): DIAGNOSTIC #
################################

mask_primers_NDRG4 <- (dat_correl[dat_NDRG4_ROC_logit_perCpG$ID,]$dist==0) & (!is.na(dat_correl[dat_NDRG4_ROC_logit_perCpG$ID,]$dist==0))

tiff(paste(output_dir,"NDRG4_AIC.jpg",sep=""),width=2048+512+256,height=1024*2/3)
par(mfrow=c(2,1))
plotPlots_AIC(dat_gene=dat_NDRG4_ROC_logit_perCpG,
              locs_gene=locs_NDRG4, gene="NDRG4",
              MSP_amplicon=c(58497439,58497541),
              primer_annot_gene=primer_annot_NDRG4,
              gene_obj=genes_struct,
              chrom=as.character(data.frame(locs_NDRG4)[1,1]),
              transcript_ids=unique(genes_struct$ensembl_transcript_id[genes_struct$external_gene_name=="NDRG4"]),
              elevate_amplicon=rbind(c(4,2),c(8,2),c(10,1),c(12,1),c(13,2)),
              primer_mask=mask_primers_NDRG4,
              max_y=36)
dev.off()

####################################
# Plot melanoma (LY75): PROGNOSTIC #
####################################

mask_primers_LY75 <- (dat_correl[dat_LY75_ROC_logit_perCpG$ID,]$dist==0) & (!is.na(dat_correl[dat_LY75_ROC_logit_perCpG$ID,]$dist==0))

tiff(paste(output_dir,"LY75_AIC.jpg",sep=""),width=2048+512+256,height=1024*2/3)
par(mfrow=c(2,1))
plotPlots_AIC(dat_gene=dat_LY75_ROC_logit_perCpG,
              locs_gene=locs_LY75, gene="LY75",
              MSP_amplicon=c(160761053,160761124),
              primer_annot_gene=primer_annot_LY75,
              gene_obj=genes_struct,
              chrom=as.character(data.frame(locs_LY75)[1,1]),
              transcript_ids=unique(genes_struct$ensembl_transcript_id[genes_struct$external_gene_name=="LY75"]),
              # elevate_amplicon=c(6,2),
              primer_mask=mask_primers_LY75,
              max_y=16)
dev.off()

##################################
# Plot ccRCC (GREM1): PROGNOSTIC #
##################################

mask_primers_GREM1 <- (dat_correl[dat_GREM1_ROC_logit_perCpG$ID,]$dist==0) & (!is.na(dat_correl[dat_GREM1_ROC_logit_perCpG$ID,]$dist==0))

tiff(paste(output_dir,"GREM1_AIC.jpg",sep=""),width=2048+512+256,height=1024*2/3)
par(mfrow=c(2,1))
plotPlots_AIC(dat_gene=dat_GREM1_ROC_logit_perCpG,
          locs_gene=locs_GREM1, gene="GREM1",
          MSP_amplicon=c(33009978,33010025),
          primer_annot_gene=primer_annot_GREM1,
          gene_obj=genes_struct,
          chrom=as.character(data.frame(locs_GREM1)[1,1]),
          transcript_ids=unique(genes_struct$ensembl_transcript_id[genes_struct$external_gene_name=="GREM1"]),
          elevate_amplicon=rbind(c(4,1),c(8,2),c(10,1)),
          primer_mask=mask_primers_GREM1,
          max_y=25)
dev.off()

################
# Checks Johan #
################

jpeg("NDRG4_av_cov.jpg")
plot(data.frame(locs_NDRG4)$start,data.frame(locs_NDRG4)$av_cov)
dev.off()

jpeg("GREM1_av_cov.jpg")
plot(data.frame(locs_GREM1)$start,data.frame(locs_GREM1)$av_cov)
dev.off()

jpeg("LY75_av_cov.jpg")
plot(data.frame(locs_LY75)$start,data.frame(locs_LY75)$av_cov)
dev.off()


# #################
# # Checks Johan #
# ################
# 
# # GREM1
# dat_GREM1_notlocated <- cbind(data.frame(locs_GREM1)[,1:2],
#                               rowMeans(totalReads(bs_GREM1)),
#                               rowMeans(methReads(bs_GREM1)),
#                               rep(F,length(locs_GREM1)))
# for (i in 1:nrow(primer_annot_GREM1)){
#   dat_GREM1_notlocated[(dat_GREM1_ROC_logit_perCpG$Pos>primer_annot_GREM1$begin[i]) &
#                          (dat_GREM1_ROC_logit_perCpG$Pos<primer_annot_GREM1$end[i]),5] <- T
# }
# colnames(dat_GREM1_notlocated) <- c("chr","pos","coverage","meth_cov","located")
# head(dat_GREM1_notlocated)
# as.character(dat_GREM1_notlocated$chr[1])
# min(dat_GREM1_notlocated$pos)
# max(dat_GREM1_notlocated$pos)
#                               
# boxplot(coverage~located,dat_GREM1_notlocated)
# sum(!dat_GREM1_notlocated$located)
# mean(dat_GREM1_notlocated$coverage[dat_GREM1_notlocated$located])
# mean(dat_GREM1_notlocated$coverage[!dat_GREM1_notlocated$located])
# t.test(dat_GREM1_notlocated$coverage[dat_GREM1_notlocated$located],
#        dat_GREM1_notlocated$coverage[!dat_GREM1_notlocated$located])
# 
# # NDRG4
# dat_NDRG4_notlocated <- cbind(data.frame(locs_NDRG4)[,1:2],
#                               rowMeans(totalReads(bs_NDRG4)),
#                               rowMeans(methReads(bs_NDRG4)),
#                               rep(F,length(locs_NDRG4)))
# for (i in 1:nrow(primer_annot_NDRG4)){
#   dat_NDRG4_notlocated[(dat_NDRG4_ROC_logit_perCpG$Pos>primer_annot_NDRG4$begin[i]) &
#                          (dat_NDRG4_ROC_logit_perCpG$Pos<primer_annot_NDRG4$end[i]),5] <- T
# }
# colnames(dat_NDRG4_notlocated) <- c("chr","pos","coverage","meth_cov","located")
# head(dat_NDRG4_notlocated)
# as.character(dat_NDRG4_notlocated$chr[1])
# min(dat_NDRG4_notlocated$pos)
# max(dat_NDRG4_notlocated$pos)
# 
# boxplot(coverage~located,dat_NDRG4_notlocated)
# sum(!dat_NDRG4_notlocated$located)
# mean(dat_NDRG4_notlocated$coverage[dat_NDRG4_notlocated$located])
# mean(dat_NDRG4_notlocated$coverage[!dat_NDRG4_notlocated$located])
# t.test(dat_NDRG4_notlocated$coverage[dat_NDRG4_notlocated$located],
#        dat_NDRG4_notlocated$coverage[!dat_NDRG4_notlocated$located])


######################################
######################################
######################################
###                                ###
###                                ###
###   Intra amplicon differences   ###
###                                ###
###                                ###
######################################
######################################
######################################

source("/data/homes/louisc/Project_KimSmits/Intra_amplicon_diff_functions.R")

# Create directories
if(!dir.exists("AIC_plots")){dir.create("AIC_plots")}
if(!dir.exists("AIC_plots/CpG_plots")){dir.create("AIC_plots/CpG_plots")}

# all comparisons
cex_AIC_plots <- 1.25
seed_perm <- 42
fdr_treshold_toplot <- 0.10
n_bp_add <- -1

## including primer CpGs
####################

#NDRG4
pval_list_NDRG4 <- perform_AIC_analysis(primer_annot_NDRG4,Mval_NDRG4_sort,CRC_annot_sort,"NDRG4")

#GREM1
pval_list_GREM1 <- perform_AIC_analysis(primer_annot_GREM1,Mval_GREM1_sort,ccRCC_annot_sort,"GREM1")

#LY75
pval_list_LY75 <- perform_AIC_analysis(primer_annot_LY75,Mval_LY75_sort,melanoma_annot_sort,"LY75")

# Full table
part_pval_GREM1 <- cbind(rep(names(pval_list_GREM1)[pval_list_GREM1!="No CpGs"],times=as.numeric(unlist(lapply(pval_list_GREM1,nrow)))),
                         do.call(rbind,pval_list_GREM1[pval_list_GREM1!="No CpGs"]))
colnames(part_pval_GREM1) <- c("Amplicon",colnames(part_pval_GREM1[2:ncol(part_pval_GREM1)]))
rownames(part_pval_GREM1) <- 1:nrow(part_pval_GREM1)
part_pval_NDRG4 <- cbind(rep(names(pval_list_NDRG4)[pval_list_NDRG4!="No CpGs"],times=as.numeric(unlist(lapply(pval_list_NDRG4,nrow)))),
                         do.call(rbind,pval_list_NDRG4[pval_list_NDRG4!="No CpGs"]))
colnames(part_pval_NDRG4) <- c("Amplicon",colnames(part_pval_NDRG4[2:ncol(part_pval_NDRG4)]))
rownames(part_pval_NDRG4) <- 1:nrow(part_pval_NDRG4)
part_pval_LY75 <- cbind(rep(names(pval_list_LY75)[pval_list_LY75!="No CpGs"],times=as.numeric(unlist(lapply(pval_list_LY75,nrow)))),
                        do.call(rbind,pval_list_LY75[pval_list_LY75!="No CpGs"]))
colnames(part_pval_LY75) <- c("Amplicon",colnames(part_pval_LY75[2:ncol(part_pval_LY75)]))
rownames(part_pval_LY75) <- 1:nrow(part_pval_LY75)
                
pval_total_list <- rbind(part_pval_GREM1,part_pval_NDRG4,part_pval_LY75)
rownames(pval_total_list) <- 1:nrow(pval_total_list)
write.table(pval_total_list,file="pval_total_list.txt",col.names=T,row.names=F,quote=F,sep="\t")

unique(pval_total_list$Amplicon[pval_total_list$fdr<0.10])
unique(pval_total_list$Amplicon[pval_total_list$pval<0.10])


# # checks
# dat_check_GREM1 <- dat_GREM1_ROC_logit_perCpG[,c("ID","AIC")]
# dat_check_GREM1_pvals <- do.call(rbind,pval_list_GREM1)[!duplicated(do.call(rbind,pval_list_GREM1)$CpG_1),c("CpG_1","AIC_1","pval","fdr")]
# dat_check_GREM1 <- cbind(dat_check_GREM1,matrix(rep(NA,4*nrow(dat_check_GREM1)),nrow=nrow(dat_check_GREM1)))
# dat_check_GREM1[dat_check_GREM1$ID%in%dat_check_GREM1_pvals$CpG_1,3:6] <- dat_check_GREM1_pvals 
# daT_check_GREM1
# 
# dat_check_NDRG4 <- dat_NDRG4_ROC_logit_perCpG[,c("ID","AIC")]
# dat_check_NDRG4_pvals <- do.call(rbind,pval_list_NDRG4)[!duplicated(do.call(rbind,pval_list_NDRG4)$CpG_1),c("CpG_1","AIC_1","pval","fdr")]
# dat_check_NDRG4 <- cbind(dat_check_NDRG4,matrix(rep(NA,4*nrow(dat_check_NDRG4)),nrow=nrow(dat_check_NDRG4)))
# dat_check_NDRG4[dat_check_NDRG4$ID%in%dat_check_NDRG4_pvals$CpG_1,3:6] <- dat_check_NDRG4_pvals 
# dat_check_NDRG4
# 
# dat_check_LY75 <- dat_LY75_ROC_logit_perCpG[,c("ID","AIC")]
# dat_check_LY75_pvals <- do.call(rbind,pval_list_LY75)[!duplicated(do.call(rbind,pval_list_LY75)$CpG_1),c("CpG_1","AIC_1","pval","fdr")]
# dat_check_LY75 <- cbind(dat_check_LY75,matrix(rep(NA,4*nrow(dat_check_LY75)),nrow=nrow(dat_check_LY75)))
# dat_check_LY75[dat_check_LY75$ID%in%dat_check_LY75_pvals$CpG_1,3:6] <- dat_check_LY75_pvals 
# dat_check_LY75

## excluding primer CpGs
####################

# unlink("AIC_plots_red/*")
if(!dir.exists("AIC_plots_red")){dir.create("AIC_plots_red")}

#source("/data/homes/louisc/Project_KimSmits/Intra_amplicon_diff_functions.R")

pval_cutoff_red <- 0.05

reduced_AIC_analysis("GREM1",pval_list_GREM1,pval_cutoff_red,dat_GREM1_ROC_logit_perCpG)

reduced_AIC_analysis("NDRG4",pval_list_NDRG4,pval_cutoff_red,dat_NDRG4_ROC_logit_perCpG)

pval_cutoff_red <- 0.10

reduced_AIC_analysis("LY75",pval_list_LY75,pval_cutoff_red,dat_LY75_ROC_logit_perCpG)

####################################
####################################
####################################
###                              ###
###                              ###
###   Technical considerations   ###
###                              ###
###                              ###
####################################
####################################
####################################

################################
# Results for overlapping CpGs #
################################

res_dupl_CpG <- NULL

dupl_CpG_GREM1 <- gsub("_[NDRG|GREM|LY].+$","",rownames(Mval_GREM1_sort))
dupl_CpG_GREM1 <- unique(dupl_CpG_GREM1[duplicated(dupl_CpG_GREM1)])
for (i in 1:length(dupl_CpG_GREM1)){
  dat_dupl_CpG_GREM1 <- Mval_GREM1_sort[grep(dupl_CpG_GREM1[i],rownames(Mval_GREM1_sort)),]
  beta_dupl_CpG_GREM1 <- betas_GREM1_sort[grep(dupl_CpG_GREM1[i],rownames(betas_GREM1_sort)),]
  
  res_dupl_CpG <- rbind(res_dupl_CpG,
                              c("GREM1",dupl_CpG_GREM1[i],
                                gsub(paste0(dupl_CpG_GREM1[i],"_"),"",rownames(dat_dupl_CpG_GREM1)),
                                data.frame(locs_GREM1)$av_cov[data.frame(locs_GREM1)$seqnames==gsub("_.+$","",dupl_CpG_GREM1[i]) & data.frame(locs_GREM1)$start==gsub("^.+_","",dupl_CpG_GREM1[i])],
                                t.test(dat_dupl_CpG_GREM1[1,],dat_dupl_CpG_GREM1[2,])$p.val,
                                mean(beta_dupl_CpG_GREM1[1,],na.rm=T),mean(beta_dupl_CpG_GREM1[2,],na.rm=T),
                                mean(beta_dupl_CpG_GREM1[1,],na.rm=T)-mean(beta_dupl_CpG_GREM1[2,],na.rm=T)))
}
dupl_CpG_LY75 <- gsub("_[NDRG|GREM|LY].+$","",rownames(Mval_LY75_sort))
dupl_CpG_LY75 <- unique(dupl_CpG_LY75[duplicated(dupl_CpG_LY75)])
for (i in 1:length(dupl_CpG_LY75)){
  dat_dupl_CpG_LY75 <- Mval_LY75_sort[grep(dupl_CpG_LY75[i],rownames(Mval_LY75_sort)),]
  beta_dupl_CpG_LY75 <- betas_LY75_sort[grep(dupl_CpG_LY75[i],rownames(betas_LY75_sort)),]
  
  res_dupl_CpG <- rbind(res_dupl_CpG,
                        c("LY75",dupl_CpG_LY75[i],
                          gsub(paste0(dupl_CpG_LY75[i],"_"),"",rownames(dat_dupl_CpG_LY75)),
                          data.frame(locs_LY75)$av_cov[data.frame(locs_LY75)$seqnames==gsub("_.+$","",dupl_CpG_LY75[i]) & data.frame(locs_LY75)$start==gsub("^.+_","",dupl_CpG_LY75[i])],
                          t.test(dat_dupl_CpG_LY75[1,],dat_dupl_CpG_LY75[2,])$p.val,
                          mean(beta_dupl_CpG_LY75[1,],na.rm=T),mean(beta_dupl_CpG_LY75[2,],na.rm=T),
                          mean(beta_dupl_CpG_LY75[1,],na.rm=T)-mean(beta_dupl_CpG_LY75[2,],na.rm=T)))
}
dupl_CpG_NDRG4 <- gsub("_[NDRG|GREM|LY].+$","",rownames(Mval_NDRG4_sort))
dupl_CpG_NDRG4 <- unique(dupl_CpG_NDRG4[duplicated(dupl_CpG_NDRG4)])
for (i in 1:length(dupl_CpG_NDRG4)){
  dat_dupl_CpG_NDRG4 <- Mval_NDRG4_sort[grep(dupl_CpG_NDRG4[i],rownames(Mval_NDRG4_sort)),]
  beta_dupl_CpG_NDRG4 <- betas_NDRG4_sort[grep(dupl_CpG_NDRG4[i],rownames(betas_NDRG4_sort)),]
  
  res_dupl_CpG <- rbind(res_dupl_CpG,
                        c("NDRG4",dupl_CpG_NDRG4[i],
                          gsub(paste0(dupl_CpG_NDRG4[i],"_"),"",rownames(dat_dupl_CpG_NDRG4)),
                          data.frame(locs_NDRG4)$av_cov[data.frame(locs_NDRG4)$seqnames==gsub("_.+$","",dupl_CpG_NDRG4[i]) & data.frame(locs_NDRG4)$start==gsub("^.+_","",dupl_CpG_NDRG4[i])],
                          t.test(dat_dupl_CpG_NDRG4[1,],dat_dupl_CpG_NDRG4[2,])$p.val,
                          mean(beta_dupl_CpG_NDRG4[1,],na.rm=T),mean(beta_dupl_CpG_NDRG4[2,],na.rm=T),
                          mean(beta_dupl_CpG_NDRG4[1,],na.rm=T)-mean(beta_dupl_CpG_NDRG4[2,],na.rm=T)))
}
res_dupl_CpG <- data.frame(res_dupl_CpG)
colnames(res_dupl_CpG) <- c("Gene","CpG","Amplicon1","Amplicon2","A1_av_cov","A2_av_cov","pval","A1_av_meth","A2_av_meth","diff_meth")
res_dupl_CpG$pval <- as.numeric(res_dupl_CpG$pval)
res_dupl_CpG$diff_meth <- as.numeric(res_dupl_CpG$diff_meth)
write.table(res_dupl_CpG,file="res_dupl_CpG.txt",col.names = T,row.names = F,quote = F,sep="\t")

res_dupl_CpG$A1_fwd_meth <- NA
res_dupl_CpG$A1_rev_meth <- NA
res_dupl_CpG$A2_fwd_meth <- NA
res_dupl_CpG$A2_rev_meth <- NA
for (i in 1:nrow(res_dupl_CpG)){
  meth_A1_fwd_i <- primer_annot_list[[res_dupl_CpG$Gene[i]]][primer_annot_list[[res_dupl_CpG$Gene[i]]]$name==res_dupl_CpG$Amplicon1[i],grepl("meth_CpG[1-3]_fwd",colnames(primer_annot_list[[res_dupl_CpG$Gene[i]]]))]
  meth_A1_rev_i <- primer_annot_list[[res_dupl_CpG$Gene[i]]][primer_annot_list[[res_dupl_CpG$Gene[i]]]$name==res_dupl_CpG$Amplicon1[i],grepl("meth_CpG[1-3]_rev",colnames(primer_annot_list[[res_dupl_CpG$Gene[i]]]))]
  meth_A2_fwd_i <- primer_annot_list[[res_dupl_CpG$Gene[i]]][primer_annot_list[[res_dupl_CpG$Gene[i]]]$name==res_dupl_CpG$Amplicon2[i],grepl("meth_CpG[1-3]_fwd",colnames(primer_annot_list[[res_dupl_CpG$Gene[i]]]))]  
  meth_A2_rev_i <- primer_annot_list[[res_dupl_CpG$Gene[i]]][primer_annot_list[[res_dupl_CpG$Gene[i]]]$name==res_dupl_CpG$Amplicon2[i],grepl("meth_CpG[1-3]_rev",colnames(primer_annot_list[[res_dupl_CpG$Gene[i]]]))]

  if (sum(!is.na(meth_A1_fwd_i))>0){
    meth_A1_fwd_i <- meth_A1_fwd_i[!is.na(meth_A1_fwd_i)]
    res_dupl_CpG$A1_fwd_meth[i] <- paste(round(meth_A1_fwd_i,3),collapse=",")
  }
  if (sum(!is.na(meth_A1_rev_i))>0){
    meth_A1_rev_i <- meth_A1_rev_i[!is.na(meth_A1_rev_i)]
    res_dupl_CpG$A1_rev_meth[i] <- paste(round(meth_A1_rev_i,3),collapse=",")
  }
  if (sum(!is.na(meth_A2_fwd_i))>0){
    meth_A2_fwd_i <- meth_A2_fwd_i[!is.na(meth_A2_fwd_i)]
    res_dupl_CpG$A2_fwd_meth[i] <- paste(round(meth_A2_fwd_i,3),collapse=",")
  }
  if (sum(!is.na(meth_A2_rev_i))>0){
    meth_A2_rev_i <- meth_A2_rev_i[!is.na(meth_A2_rev_i)]
    res_dupl_CpG$A2_rev_meth[i] <- paste(round(meth_A2_rev_i,3),collapse=",")
  }

}

res_dupl_CpG$A1_deltaAIC <- NA
res_dupl_CpG$A2_deltaAIC <- NA
res_dupl_CpG$A1_dist_to_primer_CpG <- NA
res_dupl_CpG$A2_dist_to_primer_CpG <- NA
res_dupl_cpG$A1_av_meth <- NA
res_dupl_cpG$A2_av_meth <- NA
for (i in 1:nrow(res_dupl_CpG)){
  dat_i <- get(paste0("dat_",res_dupl_CpG$Gene[i],"_ROC_logit_perCpG"))

  res_dupl_CpG$A1_deltaAIC[i] <- dat_i$delta_AIC[dat_i$ID==paste(res_dupl_CpG$CpG[i],res_dupl_CpG$Amplicon1[i],sep="_")]
  res_dupl_CpG$A2_deltaAIC[i] <- dat_i$delta_AIC[dat_i$ID==paste(res_dupl_CpG$CpG[i],res_dupl_CpG$Amplicon2[i],sep="_")]

  res_dupl_CpG$A1_dist_to_primer_CpG[i] <- dat_correl$dist_CpG[dat_correl$ID==paste(res_dupl_CpG$CpG[i],res_dupl_CpG$Amplicon1[i],sep="_")]
  res_dupl_CpG$A2_dist_to_primer_CpG[i] <- dat_correl$dist_CpG[dat_correl$ID==paste(res_dupl_CpG$CpG[i],res_dupl_CpG$Amplicon2[i],sep="_")]
}

write.table(res_dupl_CpG,file="res_dupl_CpG.txt",col.names = T,row.names = F,quote = F,sep="\t")


################################
# Primer-specific CpG analysis #
################################

# Select amplicons that have primer-specific CpGs (including primer-spec CpGs)
amplicons_psCpGs <- unique(gsub("^[0-9]+_[0-9]+_","",rownames(dat_correl[(dat_correl$dist_CpG==0) & !is.na(dat_correl$dist_CpG),])))
dat_diff_psCpGs <- NULL
for (i in 1:length(amplicons_psCpGs)){
  dat_diff_psCpGs <- rbind(dat_diff_psCpGs,
                           c(amplicons_psCpGs[i],sum(grepl(paste0(amplicons_psCpGs[i],"$"),rownames(dat_correl)) & (dat_correl$dist_CpG==0) & (!is.na(dat_correl$dist_CpG))),
                             mean(dat_correl$deltaAIC[grepl(paste0(amplicons_psCpGs[i],"$"),rownames(dat_correl)) & (dat_correl$dist_CpG==0) & (!is.na(dat_correl$dist_CpG))]),
                             mean(dat_correl$deltaAIC[grepl(paste0(amplicons_psCpGs[i],"$"),rownames(dat_correl)) & (dat_correl$dist_CpG!=0) & (!is.na(dat_correl$dist_CpG))])))
}
dat_diff_psCpGs <- data.frame(dat_diff_psCpGs)
colnames(dat_diff_psCpGs) <- c("amplicon","n_CpGs","av_dAIC_in_primer","av_dAIC_out_primer")
for (i in 2:4){
  dat_diff_psCpGs[,i] <- as.numeric(dat_diff_psCpGs[,i])
}
dat_diff_psCpGs$diff <- dat_diff_psCpGs$av_dAIC_in_primer - dat_diff_psCpGs$av_dAIC_out_primer
write.table(dat_diff_psCpGs,file="dat_diff_psCpGs.txt",col.names=T, row.names=F,sep="\t",quote=F)

# Select target specific CpGs
dat_correl_not_in_primer <- dat_correl[(dat_correl$dist != 0) | is.na(dat_correl$dist),]
dim(dat_correl_not_in_primer)

sum(((dat_correl$dist != 0) | is.na(dat_correl$dist)) & (dat_correl$CG_qual==0))
sum(((dat_correl$dist != 0) | is.na(dat_correl$dist)) & (dat_correl$CG_qual==1))
sum(((dat_correl$dist != 0) | is.na(dat_correl$dist)) & (dat_correl$CG_qual==2))

# deltaAIC ifo quantitative presence of CpG in primer sequence
jpeg("/data/homes/louisc/Project_KimSmits/violin_deltaAIC_CpG_quant.jpg")
g <- ggplot(dat_correl_not_in_primer, aes(x=deltaAIC, y=CG_quant)) + 
  geom_violin(trim=FALSE) + 
  theme(axis.text.x= element_text(size=14),
        axis.text.y= element_text(size=14))
print(g)
dev.off()

# deltaAIC ifo presence of CpG in primer sequence
jpeg("/data/homes/louisc/Project_KimSmits/violin_deltaAIC_CpG_qual.jpg")
g <- ggplot(dat_correl_not_in_primer, aes(x=deltaAIC, y=CG_qual)) + 
  geom_violin(trim=FALSE) + 
  theme(axis.text.x= element_text(size=14),
        axis.text.y= element_text(size=14))
print(g)
dev.off()

# Statistical testing
t.test(dat_correl_not_in_primer$deltaAIC[dat_correl_not_in_primer$CG_qual==0],
       dat_correl_not_in_primer$deltaAIC[dat_correl_not_in_primer$CG_qual==1])

t.test(dat_correl_not_in_primer$deltaAIC[dat_correl_not_in_primer$CG_qual==0],
       dat_correl_not_in_primer$deltaAIC[dat_correl_not_in_primer$CG_qual==2])
 

# Select target specific CpGs that have at least one primer-specific CpGs in one of their primers
dat_correl_filt <- dat_correl[(dat_correl$CG_qual!=0) & (dat_correl$dist!=0) & (!is.na(dat_correl$dist)),]
dim(dat_correl_filt)

# deltaAIC ifo deviation of 50% average DNA methylation
jpeg("/data/homes/louisc/Project_KimSmits/correl_deltaAIC_diffmeth.jpg")
g <- ggplot(dat_correl_filt[!(dat_correl_filt$dist==0),], aes(x=diff_meth, y=deltaAIC)) + 
  geom_point()+ 
  geom_smooth(method=lm)
print(g)
dev.off()

cor.test(dat_correl_filt[!(dat_correl_filt$dist==0),]$diff_meth,
         dat_correl_filt[!(dat_correl_filt$dist==0),]$deltaAIC)
                         
# deltaAIC ifo distance to nearest primer
jpeg("/data/homes/louisc/Project_KimSmits/correl_deltaAIC_dist.jpg")
g <- ggplot(dat_correl_filt[!(dat_correl_filt$dist==0),], aes(x=dist_CpG, y=deltaAIC)) + 
  geom_point()+ 
  geom_smooth(method=lm)  + 
  theme(axis.text.x= element_text(size=14),
        axis.text.y= element_text(size=14))
print(g)
dev.off()         

cor.test(dat_correl_filt[!(dat_correl_filt$dist==0),]$dist_CpG,
         dat_correl_filt[!(dat_correl_filt$dist==0),]$deltaAIC)


# deltaAIC ifo difference in DNA methylation of nearest primer CpG
jpeg("/data/homes/louisc/Project_KimSmits/correl_deltaAIC_diffmeth_nearestprimerCpG.jpg")
g <- ggplot(dat_correl_filt[!(dat_correl_filt$dist==0),], aes(x=abs(av_meth_primer_CpG-0.5), y=deltaAIC)) + 
  geom_point()+ 
  geom_smooth(method=lm) + 
  theme(axis.text.x= element_text(size=14),
        axis.text.y= element_text(size=14))
print(g)
dev.off()         

cor.test(abs(dat_correl_filt[!(dat_correl_filt$dist==0),]$av_meth_primer_CpG-0.5),
         dat_correl_filt[!(dat_correl_filt$dist==0),]$deltaAIC)


# deltaAIC ifo average methylation rate
jpeg("/data/homes/louisc/Project_KimSmits/correl_deltaAIC_avMeth.jpg")
g <- ggplot(dat_correl, aes(x=av_meth_CpG, y=deltaAIC)) + 
  geom_point()+
  geom_smooth(method=lm)
print(g)
dev.off()

cor.test(dat_correl$av_meth_CpG,dat_correl$deltaAIC)


(paste(res_dupl_CpG$CpG,res_dupl_CpG$Amplicon1,sep="_") %in% dat_correl$ID[dat_correl$dist_CpG==0 & !is.na(dat_correl$dist_CpG)]) | (paste(res_dupl_CpG$CpG,res_dupl_CpG$Amplicon2,sep="_") %in% dat_correl$ID[dat_correl$dist_CpG==0 & !
is.na(dat_correl$dist_CpG)])


###############################
###############################
###############################
###                         ###
###                         ###   
###   Read based analysis   ###
###                         ###
###                         ###
###############################
###############################
###############################


######################################################
# Extract read based data from wgbs_tools .pat files #
######################################################

setwd("/data/homes/louisc/Project_KimSmits")

CpG_info_wgbstools <- data.frame(fread("/data/homes/louisc/Project_KimSmits/wgbs_tools/references/hg37_87/CpG.bed.gz",header=F))
colnames(CpG_info_wgbstools) <- c("chr","pos","id")

if(!dir.exists("waterfall_plots")){dir.create("waterfall_plots")}

all_files <- list.files("output_wgbstools")
all_samples <- unique(c(paste0(paste0("run",gsub("_.*$","",gsub("^.*run_","",all_files[grepl("run_239",all_files)]))),"_",gsub("_.*$","",gsub("^.*_S","",all_files[grepl("run_239",all_files)]))),
                        paste0("run",gsub("_S.*$","",gsub("^.*run_","",all_files[!grepl("run_239",all_files)])))))

ov_rb_anal <- data.frame(matrix(rep(NA,length(amplicons)*length(all_samples)),ncol=length(all_samples)))
colnames(ov_rb_anal) <- all_samples
rownames(ov_rb_anal) <- amplicons

ccRCC_annot_sort$surv <- Surv(ccRCC_annot_sort$Followup_time,ccRCC_annot_sort$NK_death)
melanoma_annot_sort$surv <- Surv(melanoma_annot_sort$Followup_time,melanoma_annot_sort$MK_death)

# Extract data from .pat files
unlink("waterfall_plots/sink_rb_analysis.log")
sink("waterfall_plots/sink_rb_analysis.log")
counter <- 1
db_rb_anal <- list()
for (amplicon in amplicons){
  print(sprintf("Starting read based analysis for %s",amplicon))

  if (amplicon%in%c("LY75_6_S","GREM1_60","NDRG4_92")){
      print(sprintf("Amplicon %s skipped due to various reasons (coverage, primer design)",amplicon))
      next()
  }

  #if(!dir.exists(paste0("waterfall_plots/",amplicon))){dir.create(paste0("waterfall_plots/",amplicon))}

  files <- list.files("output_wgbstools")
  if (amplicon=="GREM1_Original"){
    files <- files[grepl("GREM1_Manon",files) & grepl("pat.gz$",files)]
  } else if (amplicon=="NDRG4_Original"){
    files <- files[grepl("NDRG4_FPassay",files) & grepl("pat.gz$",files)]
  } else {
    files <- files[grepl(amplicon,files) & grepl("pat.gz$",files)]
  }

  start_CpG_ampl <- CpG_info_wgbstools$id[CpG_info_wgbstools$chr==primer_annot$chromosome[primer_annot$name==amplicon] &
                    CpG_info_wgbstools$pos>as.numeric(primer_annot$begin[primer_annot$name==amplicon])][1]

  # original one, changed due to fact that LY75_9_S has a CpG in bisulf data that is not present in wgbstools annotation => mutation?
  #n_CpG_ampl <- sum(grepl(amplicon,rownames(dat_correl)))
  if (grepl("LY75",amplicon)){
    mask_n_CpG <- c(as.numeric(gsub("^.+_","",gsub("_[A-Z]+.+$","",rownames(dat_correl)[grepl(amplicon,rownames(dat_correl))])))-1)%in%head(CpG_info_wgbstools$pos[CpG_info_wgbstools$id>=start_CpG_ampl],50)
  } else {
    mask_n_CpG <- as.numeric(gsub("^.+_","",gsub("_[A-Z]+.+$","",rownames(dat_correl)[grepl(amplicon,rownames(dat_correl))])))%in%head(CpG_info_wgbstools$pos[CpG_info_wgbstools$id>=start_CpG_ampl],50)
  }
  n_CpG_ampl <- sum(mask_n_CpG)
  print(sprintf("%s CpGs in amplicon",n_CpG_ampl))

  db_rb_anal[[counter]] <- list()
  names(db_rb_anal)[counter] <- amplicon

  counter_j <- 1
  for (j in files){
    if(grepl("run_239",j)){
      sample_j_ori <- paste0(paste0("run",gsub("_.*$","",gsub("^.*run_","",j))),"_",gsub("_.*$","",gsub("^.*_S","",j)))
    } else {
      sample_j_ori <- paste0("run",gsub("_S.*$","",gsub("^.*run_","",j)))
    }
    wh_sample_j <- grep(sample_j_ori,sample_annot$Sample_ID)
    # revert to actual sample name (i.e. based on pathology)?
    index_col_ov <- colnames(ov_rb_anal)==sample_j_ori 

    # skip samples that were excluded
    if(length(wh_sample_j)==0){      
        print(sprintf("%s sample skipped for %s",sample_j_ori,amplicon))
        next()
    }

    # skip samples that did not contain any data for the respective cancer
    feat_oi_j <- NA
    if(grepl("LY75",amplicon)){
      if(sample_annot$melanoma[wh_sample_j]==""){
        print(sprintf("%s sample skipped for %s",sample_j_ori,amplicon))
        next()
      } else  {
        sample_j <- sample_annot$melanoma[wh_sample_j]
        feat_oi_j <- melanoma_annot_sort$surv[melanoma_annot_sort$Sample_ID==sample_j]
      }
    } else if (grepl("GREM1",amplicon)){
      if(sample_annot$ccRCC[wh_sample_j]==""){
        print(sprintf("%s sample skipped for %s",sample_j_ori,amplicon))
        next()
      } else  {
        sample_j <- sample_annot$ccRCC[wh_sample_j]
        feat_oi_j <- ccRCC_annot_sort$surv[ccRCC_annot_sort$Sample_ID==sample_j]
      }
    } else if (grepl("NDRG4",amplicon)){
      if(sample_annot$CRC[wh_sample_j]==""){
        print(sprintf("%s sample skipped for %s",sample_j_ori,amplicon))
        next()
      } else  {
        sample_j <- sample_annot$CRC[wh_sample_j]
        feat_oi_j <- CRC_annot_sort$Tissue[CRC_annot_sort$Sample_ID==sample_j]
      }
    }

    print(sprintf("%s sample (%s) started, formerly known as %s.",sample_j_ori,wh_sample_j,sample_j))

    rb_info <- fread(paste0("output_wgbstools/",j),header=F)
    rb_info <- data.frame(rb_info)
    colnames(rb_info) <- c("chr","pos","allele","count")

    rb_info_filt <- rb_info[rb_info$pos==start_CpG_ampl,]
    
    # if (sum(grepl(amplicon,rownames(dat_correl)))!=n_CpG_ampl){
    #   print(sprintf("Non matching CpG dim, skipping %s",amplicon))
    # }

    mask <- nchar(rb_info_filt$allele)==n_CpG_ampl & !grepl("\\.",rb_info_filt$allele)
    rb_info_compl <- rb_info_filt[mask,]

    if(sum(rb_info_compl$count)<10){
      print(sprintf("Sample %s (%s) skipped due too low coverage",sample_j,sample_j_ori))
      next()
    }

    print(sprintf("%s out of %s reads retained after filtering (%s pct)",
                  sum(rb_info_filt$count[mask]),sum(rb_info$count),
                  round(100*sum(rb_info_filt$count[mask])/sum(rb_info$count),digits=1)))


    db_rb_anal[[counter]][[counter_j]] <- rb_info_compl
    names(db_rb_anal[[counter]])[counter_j] <- sample_j_ori

    ov_rb_anal[rownames(ov_rb_anal)==amplicon,index_col_ov] <- round(100*sum(rb_info_filt$count[mask])/sum(rb_info$count),digits=1)

    counter_j <- counter_j + 1
  }
  counter <- counter + 1 
}
sink(NULL)

write.table(ov_rb_anal,file="overview_rb_analysis.txt",col.names=T,row.names=T,quote=F,sep="\t")

add_cn <- NULL
for (i in 1:length(colnames(ov_rb_anal))){
  if (sum(sample_annot$Sample_ID==colnames(ov_rb_anal)[i])!=0){
      add_cn <- cbind(add_cn,unlist(sample_annot[sample_annot$Sample_ID==colnames(ov_rb_anal)[i],-1]))
  } else {
      add_cn <- cbind(add_cn,c("","",""))
  }
}
colnames(add_cn) <- colnames(ov_rb_anal)

write.table(add_cn,file="overview_rb_analysis.txt",row.names=T,col.name=F,quote=F,sep="\t",append=T)

# rb_info <- fread("output_wgbs/GREM1_15_run_231_17_S17.run_231_17_S17_L001_R1_001_val_1_bismark_bt2_pe.pat.gz",header=F)
# rb_info <- data.frame(rb_info)
# colnames(rb_info) <- c("chr","pos","allele","count")

##################
# Order of plots #
##################

# Orders of plots - diagnostic
order_NDRG4 <- CRC_annot_sort$Sample_ID[sort(as.character(CRC_annot_sort$Tissue),index.return=T)$ix]
for (i in 1:length(order_NDRG4)){
  order_NDRG4[i] <- sample_annot$Sample_ID[sample_annot$CRC==order_NDRG4[i]]
}
order <- list(order_GREM1,order_LY75,order_NDRG4)
# Orders of plots - prognostic
ccRCC_annot_sort$surv <- Surv(ccRCC_annot_sort$Followup_time, ccRCC_annot_sort$NK_death)
order_GREM1 <- gsub(" ","",as.character(sort(ccRCC_annot_sort$surv)))
for_ex_aequo <- 1
for (i in 1:length(order_GREM1)){
  mask_i <- gsub(" ","",as.character(ccRCC_annot_sort$surv))==order_GREM1[i]
  if (sum(mask_i)>1){
    order_GREM1[i] <- ccRCC_annot_sort$Sample_ID[mask_i][for_ex_aequo]
    for_ex_aequo <- for_ex_aequo + 1
  } else {
    order_GREM1[i] <- ccRCC_annot_sort$Sample_ID[mask_i]
    for_ex_aequo <- 1
  }
  order_GREM1[i] <- sample_annot$Sample_ID[sample_annot$ccRCC==order_GREM1[i]]
}

melanoma_annot_sort$surv <- Surv(melanoma_annot_sort$Followup_time, melanoma_annot_sort$NK_death)
order_LY75<- c(sort(as.numeric(melanoma_annot_sort$Followup_time)+c(1*10^9,0)[melanoma_annot_sort$NK_death+1]))
for_ex_aequo <- 1
for (i in 1:length(order_LY75)){
  mask_i <- as.numeric(melanoma_annot_sort$Followup_time)+c(1*10^9,0)[melanoma_annot_sort$NK_death+1]==order_LY75[i]
  if (sum(mask_i)>1){
    order_LY75[i] <- melanoma_annot_sort$Sample_ID[mask_i][for_ex_aequo]
    for_ex_aequo <- for_ex_aequo + 1
  } else {
    order_LY75[i] <- melanoma_annot_sort$Sample_ID[mask_i]
    for_ex_aequo <- 1
  }
  order_LY75[i] <- sample_annot$Sample_ID[sample_annot$melanoma==order_LY75[i]]
}
# order_LY75<- sort(as.numeric(melanoma_annot_sort$Followup_time))
# for (i in 1:length(order_LY75)){
#   order_LY75[i] <- melanoma_annot_sort$Sample_ID[as.numeric(melanoma_annot_sort$Followup_time)==order_LY75[i]]
#   order_LY75[i] <- sample_annot$Sample_ID[sample_annot$melanoma==order_LY75[i]]
# }

order_NDRG4
order_GREM1
order_LY75

##############
# MESA plots #
##############

# unlink("waterfall_plots/*.jpg")
CpG_in_primer <- rownames(dat_correl)[dat_correl$dist_CpG==0 & !is.na(dat_correl$dist_CpG)]
dat_plot_rb <- list()
counter <- 1
color_cometh <- FALSE
plot_paper <- FALSE
for (amplicon in amplicons){
  print(sprintf("Plotting for %s",amplicon))

  jpeg(paste0("waterfall_plots/",amplicon,".jpg"),width=256*10.5,height=256*6)
  par(mfrow=c(4,7))

  if(grepl("LY75",amplicon)){
    order_ampl <- order_LY75
  } else if (grepl("GREM1",amplicon)){
    order_ampl <- order_GREM1
  } else if (grepl("NDRG4",amplicon)){
    order_ampl <- order_NDRG4
  }

  dat_plot_rb[[counter]] <- list()
  names(dat_plot_rb)[counter] <- amplicon

  counter_j <- 1
  for (sample_j_ori in order_ampl){

    wh_sample_j <- grep(sample_j_ori,sample_annot$Sample_ID)

    # skip samples that did not contain any data for the respective cancer
    feat_oi_j <- NA
    if(grepl("LY75",amplicon)){
      if(sample_annot$melanoma[wh_sample_j]==""){
        print(sprintf("%s sample skipped for %s",sample_j_ori,amplicon))
        next()
      } else  {
        sample_j <- sample_annot$melanoma[wh_sample_j]
        feat_oi_j <- melanoma_annot_sort$surv[melanoma_annot_sort$Sample_ID==sample_j]
      }
    } else if (grepl("GREM1",amplicon)){
      if(sample_annot$ccRCC[wh_sample_j]==""){
        print(sprintf("%s sample skipped for %s",sample_j_ori,amplicon))
        next()
      } else  {
        sample_j <- sample_annot$ccRCC[wh_sample_j]
        feat_oi_j <- ccRCC_annot_sort$surv[ccRCC_annot_sort$Sample_ID==sample_j]
      }
    } else if (grepl("NDRG4",amplicon)){
      if(sample_annot$CRC[wh_sample_j]==""){
        print(sprintf("%s sample skipped for %s",sample_j_ori,amplicon))
        next()
      } else  {
        sample_j <- sample_annot$CRC[wh_sample_j]
        feat_oi_j <- CRC_annot_sort$Tissue[CRC_annot_sort$Sample_ID==sample_j]
      }
    }

    start_CpG_ampl <- CpG_info_wgbstools$id[CpG_info_wgbstools$chr==primer_annot$chromosome[primer_annot$name==amplicon] &
                        CpG_info_wgbstools$pos>as.numeric(primer_annot$begin[primer_annot$name==amplicon])][1]
    
    # original one, changed due to fact that LY75_9_S has a CpG in bisulf data that is not present in wgbstools annotation => mutation?
    #n_CpG_ampl <- sum(grepl(amplicon,rownames(dat_correl)))
    if (grepl("LY75",amplicon)){
      mask_n_CpG <- c(as.numeric(gsub("^.+_","",gsub("_[A-Z]+.+$","",rownames(dat_correl)[grepl(amplicon,rownames(dat_correl))])))-1)%in%head(CpG_info_wgbstools$pos[CpG_info_wgbstools$id>=start_CpG_ampl],50)
    } else {
      mask_n_CpG <- as.numeric(gsub("^.+_","",gsub("_[A-Z]+.+$","",rownames(dat_correl)[grepl(amplicon,rownames(dat_correl))])))%in%head(CpG_info_wgbstools$pos[CpG_info_wgbstools$id>=start_CpG_ampl],50)
    }
    n_CpG_ampl <- sum(mask_n_CpG)

    rb_info_compl <- db_rb_anal[[amplicon]][[sample_j_ori]]

    if(sum(rb_info_compl$count)<10){
      print(sprintf("Sample %s skipped due too low coverage",sample_j))
      next()
    }

    if(grepl("LY75",amplicon)){
      betas_ampl <- betas_LY75[grepl(paste0(amplicon,"$"),rownames(betas_LY75)),]
    } else if (grepl("GREM1",amplicon)){
      betas_ampl <- betas_GREM1[grepl(paste0(amplicon,"$"),rownames(betas_GREM1)),]
    } else if (grepl("NDRG4",amplicon)){
      betas_ampl <- betas_NDRG4[grepl(paste0(amplicon,"$"),rownames(betas_NDRG4)),]
    }

    CpG_in_primer_ampl <- which(rownames(betas_ampl)%in%CpG_in_primer)

    n_cpg_read <- NULL
    for (i in 1:nrow(rb_info_compl)){
      meth_allele <- unlist(strsplit(rb_info_compl$allele[i],""))
      if (length(CpG_in_primer_ampl)!=0){
        meth_allele[CpG_in_primer_ampl] <- "T"
      }
      n_cpg_read <- c(n_cpg_read,sum(grepl("C",meth_allele)))
    }
    order_cpgs <- order(n_cpg_read*10^9+rb_info_compl$count,decreasing=T)
    # cbind(n_cpg_read[order_cpgs],rb_info_compl$count[order_cpgs])
    rb_info_compl <- rb_info_compl[order_cpgs,] 

    pos <- CpG_info_wgbstools$pos[which(CpG_info_wgbstools$id==rb_info_compl$pos[1]):(which(CpG_info_wgbstools$id==rb_info_compl$pos[1])+n_CpG_ampl-1)]

    if (plot_paper){
      plot((pos[1]+pos[length(pos)])/2,1,
          xlim=c(pos[1]-5,pos[length(pos)]+5),
          ylim=c(0,sum(rb_info_compl$count)/1000),
          xaxt="n",
          xlab="",ylab="",
          type="n", main=sprintf("%s: %s",sample_j,feat_oi_j),cex.main=3,cex.axis=2)
    } else {
      plot((pos[1]+pos[length(pos)])/2,1,
         xlim=c(pos[1]-5,pos[length(pos)]+5),
         ylim=c(0,sum(rb_info_compl$count)/1000),
         xlab="pos",ylab="DNA meth (per 1000 reads)",
         type="n", main=sprintf("%s (%s): %s (%s CpGs)",sample_j,sample_j_ori,feat_oi_j,n_CpG_ampl),cex=2)
    }

    prev_y <- 0
    dat_plot_obj <- data.frame(matrix(rep(0,nrow(rb_info_compl)*length(unlist(strsplit(rb_info_compl$allele[i],"")))),nrow=nrow(rb_info_compl)))
    rownames(dat_plot_obj) <- cumsum(rb_info_compl$count)
    colnames(dat_plot_obj) <- 1:length(unlist(strsplit(rb_info_compl$allele[i],"")))
    for (i in 1:nrow(rb_info_compl)){
      # print(sprintf("Allele %s, prev y: %s.",rb_info_compl$allele[i],prev_y))
      meth_allele <- unlist(strsplit(rb_info_compl$allele[i],""))

      for (t in 1:length(meth_allele)){
        if(color_cometh){
          if (length(CpG_in_primer_ampl)==0){
            cometh_i <- sum(meth_allele=="C")>1
          } else {
            cometh_i <- sum(meth_allele[-CpG_in_primer_ampl]=="C")>1
          }
          for (t in 1:length(meth_allele)){
            if (meth_allele[t]=="C"){
              dat_plot_obj[i,t] <- 1
              if(t %in% CpG_in_primer_ampl){
                rect(pos[t]-0.4,prev_y,
                    pos[t]+0.4,prev_y + (rb_info_compl$count[i]/1000),
                    col="#ffed29",border=NA)
              } else {
                if (cometh_i){
                rect(pos[t]-0.4,prev_y,
                    pos[t]+0.4,prev_y + (rb_info_compl$count[i]/1000),
                    col="#46B1E1",border=NA)
                } else {
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info_compl$count[i]/1000),
                      col="#CC5500",border=NA)
                }
              }
            }
          }
        } else {
          for (t in 1:length(meth_allele)){
            if (meth_allele[t]=="C"){
              dat_plot_obj[i,t] <- 1
              if(t %in% CpG_in_primer_ampl){
                rect(pos[t]-0.4,prev_y,
                    pos[t]+0.4,prev_y + (rb_info_compl$count[i]/1000),
                    col="#ffed29",border=NA)
              } else {
                rect(pos[t]-0.4,prev_y,
                    pos[t]+0.4,prev_y + (rb_info_compl$count[i]/1000),
                    col="#CC5500",border=NA)
              }
            }
          }
        }
      }
      prev_y <- prev_y + (rb_info_compl$count[i]/1000)
    }
    dat_plot_rb[[amplicon]][[counter_j]] <- dat_plot_obj
    names(dat_plot_rb[[counter]])[counter_j] <- sample_j_ori

    dat_plot_obj[,CpG_in_primer_ampl] <- 0
    
    peak_height_i <- NULL
    for (t in 1:ncol(dat_plot_obj)){
      if(sum(dat_plot_obj[,t]==0)==0){
        peak_height_i <- c(peak_height_i,rownames(dat_plot_obj)[nrow(dat_plot_obj)])
      } else if (which(dat_plot_obj[,t]==0)[1]==1){
        peak_height_i <- c(peak_height_i,0)
      } else{      
        peak_height_i <- c(peak_height_i,rownames(dat_plot_obj)[which(dat_plot_obj[,t]==0)[1]-1])
      }
    }
    #print(peak_height_i)
    text(pos[which.max(as.numeric(peak_height_i))],sum(rb_info_compl$count)/1000*0.9,labels="*",cex=4)

    counter_j <- counter_j + 1
  }
  counter <- counter + 1
  dev.off()
}

###########################################
# Compare wgbs_tools data with BiSeq data #
###########################################

unlink("comparison_WGBStools_BiSeq.log")
sink("comparison_WGBStools_BiSeq.log")
for (amplicon in amplicons){
  methrate <- NULL
  print(sprintf("Checking for %s",amplicon))

  start_CpG_ampl <- CpG_info_wgbstools$id[CpG_info_wgbstools$chr==primer_annot$chromosome[primer_annot$name==amplicon] &
                    CpG_info_wgbstools$pos>as.numeric(primer_annot$begin[primer_annot$name==amplicon])][1]

  for (sample_j_ori in names(db_rb_anal[[amplicon]])){

    if(grepl("LY75",amplicon)){
      betas_ampl <- betas_LY75[grepl(paste0(amplicon,"$"),rownames(betas_LY75)),]
      sample_j <- sample_annot$melanoma[sample_annot$Sample_ID==sample_j_ori]
    } else if (grepl("GREM1",amplicon)){
      betas_ampl <- betas_GREM1[grepl(paste0(amplicon,"$"),rownames(betas_GREM1)),]
      sample_j <- sample_annot$ccRCC[sample_annot$Sample_ID==sample_j_ori]
    } else if (grepl("NDRG4",amplicon)){
      betas_ampl <- betas_NDRG4[grepl(paste0(amplicon,"$"),rownames(betas_NDRG4)),]
      sample_j <- sample_annot$CRC[sample_annot$Sample_ID==sample_j_ori]
    }
    CpG_in_primer_ampl <- which(rownames(betas_ampl)%in%CpG_in_primer)

    if(sample_j==""){
      print(sprintf("%s for %s skipped because not annotated in sample annotation.",sample_j_ori,amplicon))
      next()
    }

    dat_sample_j <- db_rb_anal[[amplicon]][[sample_j_ori]]

    # Check rb dna meth counts vs biseq
    methrate <- NULL
    for (i in 1:length(unlist(strsplit(dat_sample_j$allele[1],"")))){
      if (i == 1){
        methrate <- c(methrate,sum(dat_sample_j$count[grepl("^C.*$",dat_sample_j$allele)])/sum(dat_sample_j$count))
      } else if (i == nrow(dat_sample_j)){
        methrate <- c(methrate,sum(dat_sample_j$count[grepl("^.*C$",dat_sample_j$allele)])/sum(dat_sample_j$count))
      } else 
        methrate <- c(methrate,sum(dat_sample_j$count[grepl(paste0("^",paste(rep(".",i-1),collapse=""),"C.*$"),dat_sample_j$allele)])/sum(dat_sample_j$count))
    }

    if (grepl("LY75",amplicon)){
      mask_n_CpG <- c(as.numeric(gsub("^.+_","",gsub("_[A-Z]+.+$","",rownames(dat_correl)[grepl(paste0(amplicon,"$"),rownames(dat_correl))])))-1)%in%head(CpG_info_wgbstools$pos[CpG_info_wgbstools$id>=start_CpG_ampl],50)
    } else {
      mask_n_CpG <- as.numeric(gsub("^.+_","",gsub("_[A-Z]+.+$","",rownames(dat_correl)[grepl(paste0(amplicon,"$"),rownames(dat_correl))])))%in%head(CpG_info_wgbstools$pos[CpG_info_wgbstools$id>=start_CpG_ampl],50)
    }

    if(sum(is.na(betas_ampl[mask_n_CpG,sample_j]))==length(betas_ampl[mask_n_CpG,sample_j])){
      print(sprintf("%s for %s skipped because too much NA's.",sample_j,amplicon))
      next()
    }
    if(length(betas_ampl[mask_n_CpG,sample_j])<=2){
      print(sprintf("%s for %s skipped because 2 CpGs.",sample_j,amplicon))
      next()
    }
    if(cor(betas_ampl[mask_n_CpG,sample_j],methrate,use = "complete.obs")<0.90){
      print(sprintf("%s for %s corr BiSeq & WGBStools: %s (p=%s)",sample_j,amplicon,
              round(cor(betas_ampl[mask_n_CpG,sample_j],methrate,use = "complete.obs"),digits=2),
              round(cor.test(betas_ampl[mask_n_CpG,sample_j],methrate,use = "complete.obs")$p.value,digits=3)))
    }

    #sum(dat_sample_j$count[grepl("^C.*$",dat_sample_j$allele)])/sum(dat_sample_j$count)  
  }
}
sink(NULL)

########################################
# Predictive value read based DNA meth #
########################################

# Extract relative peak height
meth_state <- data.frame(matrix(rep(NA,length(amplicons)*nrow(sample_annot)),nrow=nrow(sample_annot)))
colnames(meth_state) <- amplicons
rownames(meth_state) <- sample_annot$Sample_ID
meth_state2 <- meth_state
for (amplicon in amplicons){
  print(sprintf("Checking DNA meth for %s",amplicon))
  for (sample_j_ori in names(dat_plot_rb[[amplicon]])){
    
    if(grepl("LY75",amplicon)){
      betas_ampl <- betas_LY75[grepl(paste0(amplicon,"$"),rownames(betas_LY75)),]
      sample_j <- sample_annot$melanoma[sample_annot$Sample_ID==sample_j_ori]
    } else if (grepl("GREM1",amplicon)){
      betas_ampl <- betas_GREM1[grepl(paste0(amplicon,"$"),rownames(betas_GREM1)),]
      sample_j <- sample_annot$ccRCC[sample_annot$Sample_ID==sample_j_ori]
    } else if (grepl("NDRG4",amplicon)){
      betas_ampl <- betas_NDRG4[grepl(paste0(amplicon,"$"),rownames(betas_NDRG4)),]
      sample_j <- sample_annot$CRC[sample_annot$Sample_ID==sample_j_ori]
    }
    CpG_in_primer_ampl <- which(rownames(betas_ampl)%in%CpG_in_primer)

    dat_plot_rb_i <- dat_plot_rb[[amplicon]][[sample_j_ori]]
    dat_plot_rb_i[,CpG_in_primer_ampl] <- 0
    
    # Comethylation based
    if (sum(rowSums(dat_plot_rb_i)>1)==0){
      meth_state[sample_j_ori,amplicon] <- 0
      meth_state2[sample_j_ori,amplicon] <- 0
    } else {
      meth_state[sample_j_ori,amplicon] <- as.numeric(tail(rownames(dat_plot_rb_i)[rowSums(dat_plot_rb_i)>1],1))/as.numeric(tail(rownames(dat_plot_rb_i),1))

      # zone_oi <- dat_plot_rb_i[rowSums(dat_plot_rb_i)>1,]


      # # co-meth density
      zone_oi <- dat_plot_rb_i[rowSums(dat_plot_rb_i)>1,]
      # meth density
      # zone_oi <- dat_plot_rb_i[rowSums(dat_plot_rb_i)>0,]

      if (nrow(zone_oi)==1){
        zone_oi_incr <- as.numeric(rownames(zone_oi)[1])
        meth_state2[sample_j_ori,amplicon] <- sum(rowSums(zone_oi)*zone_oi_incr)/sum((ncol(zone_oi)-length(CpG_in_primer_ampl))*zone_oi_incr)
      } else {
        zone_oi_incr <- as.numeric(c(rownames(zone_oi)[1],as.numeric(rownames(zone_oi)[2:nrow(zone_oi)])-as.numeric(rownames(zone_oi)[1:(nrow(zone_oi)-1)])))
        meth_state2[sample_j_ori,amplicon] <- sum(rowSums(zone_oi)*zone_oi_incr)/sum((ncol(zone_oi)-length(CpG_in_primer_ampl))*zone_oi_incr)
      }
    }
  }
}
write.table(meth_state,file="meth_state.txt",col.names=T,row.names=T,quote=F,sep="\t")
write.table(meth_state2,file="meth_state2.txt",col.names=T,row.names=T,quote=F,sep="\t")


# Calculate AIC

avAIC_table_ext <- avAIC_table
avAIC_table_ext$rbAIC <- NA
avAIC_table_ext$rbAIC_sur <- NA
avAIC_table_ext$rbAIC_hybr <- NA

# RB AIC NDRG4
formula_CF <- formula(paste("Tissue ~ Gender+Age+Patient_ID",sep=""))
model_CF <- glm(formula_CF, 
                data = CRC_annot_sort_BS, 
                family = binomial(link="logit"))
print(summary(model_CF))

AIC_NDRG4_rb <- NULL
amplicons_NDRG4 <- amplicons[grepl("NDRG4",amplicons)]
for (amplicon in amplicons_NDRG4){
  if (amplicon=="NDRG4_92"){next()}
  # print(amplicon)
  CRC_annot_sort$meth_state <- NA
  for (i in 1:nrow(CRC_annot_sort)){
    CRC_annot_sort$meth_state[i] <- meth_state[sample_annot$Sample_ID[sample_annot$CRC==CRC_annot_sort$Sample_ID[i]],amplicon]
  }
  formula_new <- formula(paste("Tissue ~ meth_state+Patient_ID",sep=""))
  model_new <- glm(formula_new, 
                  data = CRC_annot_sort, 
                  family = binomial(link="logit"))
  # print(summary(model_new)$coefficients["meth_state",4])
  # print(AIC(model_new) - AIC(model_CF))
  CRC_annot_sort$meth_state2 <- NA
  for (i in 1:nrow(CRC_annot_sort)){
    CRC_annot_sort$meth_state2[i] <- meth_state2[sample_annot$Sample_ID[sample_annot$CRC==CRC_annot_sort$Sample_ID[i]],amplicon]
  }
  formula_new <- formula(paste("Tissue ~ meth_state2+Patient_ID",sep=""))
  model_new2 <- glm(formula_new, 
                  data = CRC_annot_sort, 
                  family = binomial(link="logit"))
  formula_new <- formula(paste("Tissue ~ meth_state+meth_state2+Patient_ID",sep=""))
  model_new3 <- glm(formula_new, 
                  data = CRC_annot_sort, 
                  family = binomial(link="logit"))
  # print(summary(model_new2)$coefficients["meth_state2",4])
  # print(AIC(model_new2) - AIC(model_CF))
  AIC_NDRG4_rb <- rbind(AIC_NDRG4_rb,c(amplicon,
                                       AIC(model_new) - AIC(model_CF),
                                       AIC(model_new2) - AIC(model_CF),
                                       AIC(model_new3) - AIC(model_CF),
                                       (as.numeric(primer_annot$end[primer_annot$name==amplicon])+as.numeric(primer_annot$begin[primer_annot$name==amplicon]))/2))
  avAIC_table_ext$rbAIC[avAIC_table_ext$AmpliconID==amplicon] <- AIC(model_new) - AIC(model_CF)
  avAIC_table_ext$rbAIC_sur[avAIC_table_ext$AmpliconID==amplicon] <- AIC(model_new2) - AIC(model_CF)
  avAIC_table_ext$rbAIC_hybr[avAIC_table_ext$AmpliconID==amplicon] <- AIC(model_new3) - AIC(model_CF)
}
AIC_NDRG4_rb

# RB AIC LY75
formula_CF <- formula(paste("Surv(Followup_time, MK_death) ~ Gender+Age",sep=""))
cox_CF <- coxph(formula_CF, data = melanoma_annot_sort_BS)
print(summary(cox_CF))

AIC_LY75_rb <- NULL
amplicons_LY75 <- amplicons[grepl("LY75",amplicons)]
for (amplicon in amplicons_LY75){
  if (amplicon=="LY75_6_S"){next()}
  melanoma_annot_sort$meth_state <- NA
  for (i in 1:nrow(melanoma_annot_sort)){
    melanoma_annot_sort$meth_state[i] <- meth_state[sample_annot$Sample_ID[sample_annot$melanoma==melanoma_annot_sort$Sample_ID[i]],amplicon]
  }
  formula_new <- formula(paste("Surv(Followup_time, NK_death) ~ meth_state",sep=""))
  cox_new <- coxph(formula_new, 
                  data = melanoma_annot_sort)
  #print(summary(cox_new))
   melanoma_annot_sort$meth_state2 <- NA
  for (i in 1:nrow(melanoma_annot_sort)){
    melanoma_annot_sort$meth_state2[i] <- meth_state2[sample_annot$Sample_ID[sample_annot$melanoma==melanoma_annot_sort$Sample_ID[i]],amplicon]
  }
  formula_new <- formula(paste("Surv(Followup_time, NK_death) ~ meth_state2",sep=""))
  cox_new2 <- coxph(formula_new, 
                  data = melanoma_annot_sort)
  formula_new <- formula(paste("Surv(Followup_time, NK_death) ~ meth_state + meth_state2",sep=""))
  cox_new3<- coxph(formula_new, 
                  data = melanoma_annot_sort)
  AIC_LY75_rb <- rbind(AIC_LY75_rb,c(amplicon,
                                     AIC(cox_new) - AIC(cox_CF),
                                     AIC(cox_new2) - AIC(cox_CF),
                                     AIC(cox_new3) - AIC(cox_CF),
                                     (as.numeric(primer_annot$end[primer_annot$name==amplicon])+as.numeric(primer_annot$begin[primer_annot$name==amplicon]))/2))
  avAIC_table_ext$rbAIC[avAIC_table_ext$AmpliconID==amplicon] <- AIC(cox_new) - AIC(cox_CF)
  avAIC_table_ext$rbAIC_sur[avAIC_table_ext$AmpliconID==amplicon] <- AIC(cox_new2) - AIC(cox_CF)
  avAIC_table_ext$rbAIC_hybr[avAIC_table_ext$AmpliconID==amplicon] <- AIC(cox_new3) - AIC(cox_CF)
}
AIC_LY75_rb

# RB AIC GREM1
formula_CF <- formula(paste("Surv(Followup_time, NK_death) ~ Gender+Age",sep=""))
cox_CF <- coxph(formula_CF, data = ccRCC_annot_sort_BS)
print(summary(cox_CF))

AIC_GREM1_rb <- NULL
amplicons_GREM1 <- amplicons[grepl("GREM1",amplicons)]
for (amplicon in amplicons_GREM1){
  if (amplicon=="GREM1_60"){next()}
  ccRCC_annot_sort$meth_state <- NA
  for (i in 1:nrow(ccRCC_annot_sort)){
    ccRCC_annot_sort$meth_state[i] <- meth_state[sample_annot$Sample_ID[sample_annot$ccRCC==ccRCC_annot_sort$Sample_ID[i]],amplicon]
  }
  formula_new <- formula(paste("Surv(Followup_time, NK_death) ~ meth_state",sep=""))
  cox_new <- coxph(formula_new, 
                  data = ccRCC_annot_sort)
  ccRCC_annot_sort$meth_state2 <- NA
  for (i in 1:nrow(ccRCC_annot_sort)){
    ccRCC_annot_sort$meth_state2[i] <- meth_state2[sample_annot$Sample_ID[sample_annot$ccRCC==ccRCC_annot_sort$Sample_ID[i]],amplicon]
  }
  formula_new <- formula(paste("Surv(Followup_time, NK_death) ~ meth_state2",sep=""))
  cox_new2<- coxph(formula_new, 
                  data = ccRCC_annot_sort)
  #print(summary(cox_new))
  formula_new <- formula(paste("Surv(Followup_time, NK_death) ~ meth_state + meth_state2",sep=""))
  cox_new3<- coxph(formula_new, 
                  data = ccRCC_annot_sort)
  AIC_GREM1_rb <- rbind(AIC_GREM1_rb,c(amplicon,
                                       AIC(cox_new) - AIC(cox_CF),
                                       AIC(cox_new2) - AIC(cox_CF),
                                       AIC(cox_new3) - AIC(cox_CF),
                                       (as.numeric(primer_annot$end[primer_annot$name==amplicon])+as.numeric(primer_annot$begin[primer_annot$name==amplicon]))/2))
  avAIC_table_ext$rbAIC[avAIC_table_ext$AmpliconID==amplicon] <- AIC(cox_new) - AIC(cox_CF)
  avAIC_table_ext$rbAIC_sur[avAIC_table_ext$AmpliconID==amplicon] <- AIC(cox_new2) - AIC(cox_CF)
  avAIC_table_ext$rbAIC_hybr[avAIC_table_ext$AmpliconID==amplicon] <- AIC(cox_new3) - AIC(cox_CF)
}
AIC_GREM1_rb

write.table(avAIC_table_ext,file="avAIC_table_ext.txt",col.names=T,row.names=F,quote=F,sep="\t")


#######################
# Supplementary Fig 3 #
#######################

if(!dir.exists("suppl_fig_3")){dir.create("suppl_fig_3")}

files_s3 <- list.files("output_wgbstools")

preamplicons_s3 <- files_s3[grepl("run_231_30",files_s3) & grepl("pat.gz$",files_s3) & grepl("pat.gz$",files_s3) & grepl("NDRG4",files_s3)]
amplicons_s3 <- unique(gsub("_run.*$","",preamplicons_s3))

files_s3 <- files_s3[(grepl("run_231_30",files_s3) | grepl("run_231_45",files_s3)) & grepl("pat.gz$",files_s3) ]

# unlink("suppl_fig_3/*.jpg")
color_cometh <- FALSE
for (amplicon in amplicons_s3){
  print(sprintf("Plotting for %s",amplicon))

  files_j <- files_s3[(grepl("run_231_30",files_s3) | grepl("run_231_45",files_s3)) &  grepl(paste0(amplicon,"_"),files_s3)]

  jpeg(paste0("suppl_fig_3/",amplicon,"_S3.jpg"),width=256*3,height=256*1.5)
  par(mfrow=c(1,2))
  par(mar=c(2,3,3,2))

  for (j in 1:length(files_j)){
    rb_info <- fread(paste0("output_wgbstools/",files_j[j]),header=F)
    rb_info <- data.frame(rb_info)
    colnames(rb_info) <- c("chr","pos","allele","count")

    start_CpG_ampl <- CpG_info_wgbstools$id[CpG_info_wgbstools$chr==primer_annot$chromosome[primer_annot$name==amplicon] &
                    CpG_info_wgbstools$pos>as.numeric(primer_annot$begin[primer_annot$name==amplicon])][1]

    mask_n_CpG <- as.numeric(gsub("^.+_","",gsub("_[A-Z]+.+$","",rownames(dat_correl)[grepl(amplicon,rownames(dat_correl))])))%in%head(CpG_info_wgbstools$pos[CpG_info_wgbstools$id>=start_CpG_ampl],50)
    n_CpG_ampl <- sum(mask_n_CpG)

    betas_ampl <- betas_NDRG4[grepl(paste0(amplicon,"$"),rownames(betas_NDRG4)),]
    CpG_in_primer_ampl <- which(rownames(betas_ampl)%in%CpG_in_primer)

    n_cpg_read <- NULL
    for (i in 1:nrow(rb_info)){
      meth_allele <- unlist(strsplit(rb_info$allele[i],""))
      if (length(CpG_in_primer_ampl)!=0){
        meth_allele[CpG_in_primer_ampl] <- "T"
      }
      n_cpg_read <- c(n_cpg_read,sum(grepl("C",meth_allele)))
    }
    order_cpgs <- order(n_cpg_read*10^9+rb_info$count,decreasing=T)
    rb_info <- rb_info[order_cpgs,] 

    pos <- CpG_info_wgbstools$pos[which(CpG_info_wgbstools$id==rb_info$pos[1]):(which(CpG_info_wgbstools$id==rb_info$pos[1])+n_CpG_ampl-1)]

    if (j==1){
      plot((pos[1]+pos[length(pos)])/2,1,
            xlim=c(pos[1]-5,pos[length(pos)]+5),
            ylim=c(0,sum(rb_info$count)),
            xaxt="n",
            xlab="",ylab="",
            type="n", main="178T",cex.main=1.5,cex.axis=1.5)

      prev_y <- 0
      for (i in 1:nrow(rb_info)){
        # print(sprintf("Allele %s, prev y: %s.",rb_info$allele[i],prev_y))
        meth_allele <- unlist(strsplit(rb_info$allele[i],""))

        for (t in 1:length(meth_allele)){
          if(color_cometh){
            if (length(CpG_in_primer_ampl)==0){
              cometh_i <- sum(meth_allele=="C")>1
            } else {
              cometh_i <- sum(meth_allele[-CpG_in_primer_ampl]=="C")>1
            }
            for (t in 1:length(meth_allele)){
              if (meth_allele[t]=="C"){
                if(t %in% CpG_in_primer_ampl){
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info$count[i]),
                      col="#ffed29",border=NA)
                } else {
                  if (cometh_i){
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info$count[i]),
                      col="#46B1E1",border=NA)
                  } else {
                    rect(pos[t]-0.4,prev_y,
                        pos[t]+0.4,prev_y + (rb_info$count[i]),
                        col="#CC5500",border=NA)
                  }
                }
              }
            }
          } else {
            for (t in 1:length(meth_allele)){
              if (meth_allele[t]=="C"){
                if(t %in% CpG_in_primer_ampl){
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info$count[i]),
                      col="#ffed29",border=NA)
                } else {
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info$count[i]),
                      col="#CC5500",border=NA)
                }
              }
              
            }
            
          }
        }
        prev_y <- prev_y + (rb_info$count[i])
      }
    } else {
      plot((pos[1]+pos[length(pos)])/2,1,
            xlim=c(pos[1]-5,pos[length(pos)]+5),
            ylim=c(0,sum(rb_info$count)/1000),
            xaxt="n",
            xlab="",ylab="",
            type="n", main="178N",cex.main=1.5,cex.axis=1.5)

      prev_y <- 0
      for (i in 1:nrow(rb_info)){
        # print(sprintf("Allele %s, prev y: %s.",rb_info$allele[i],prev_y))
        meth_allele <- unlist(strsplit(rb_info$allele[i],""))

        for (t in 1:length(meth_allele)){
          if(color_cometh){
            if (length(CpG_in_primer_ampl)==0){
              cometh_i <- sum(meth_allele=="C")>1
            } else {
              cometh_i <- sum(meth_allele[-CpG_in_primer_ampl]=="C")>1
            }
            for (t in 1:length(meth_allele)){
              if (meth_allele[t]=="C"){
                if(t %in% CpG_in_primer_ampl){
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info$count[i]/1000),
                      col="#ffed29",border=NA)
                } else {
                  if (cometh_i){
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info$count[i]/1000),
                      col="#46B1E1",border=NA)
                  } else {
                    rect(pos[t]-0.4,prev_y,
                        pos[t]+0.4,prev_y + (rb_info$count[i]/1000),
                        col="#CC5500",border=NA)
                  }
                }
              }
            }
          } else {
            for (t in 1:length(meth_allele)){
              if (meth_allele[t]=="C"){
                if(t %in% CpG_in_primer_ampl){
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info$count[i]/1000),
                      col="#ffed29",border=NA)
                } else {
                  rect(pos[t]-0.4,prev_y,
                      pos[t]+0.4,prev_y + (rb_info$count[i]/1000),
                      col="#CC5500",border=NA)
                }
              }
              
            }
            
          }
        }
        prev_y <- prev_y + (rb_info$count[i]/1000)
      }
    }
  }
  dev.off()
}
  




#################
#################
#################
###           ###
###           ###   
###   Stale   ###
###           ###
###           ###
#################
#################
#################



    # peak_height_i <- NULL
    # for (t in 1:ncol(dat_plot_rb_i)){
    #   if(sum(dat_plot_rb_i[,t]==0)==0){
    #     peak_height_i <- c(peak_height_i,rownames(dat_plot_rb_i)[nrow(dat_plot_rb_i)])
    #   } else if (which(dat_plot_rb_i[,t]==0)[1]==1){
    #     peak_height_i <- c(peak_height_i,0)
    #   } else {      
    #     peak_height_i <- c(peak_height_i,rownames(dat_plot_rb_i)[which(dat_plot_rb_i[,t]==0)[1]-1])
    #   }
    # }

    # area_i <- NULL
    # for (t in 1:ncol(dat_plot_rb_i)){
    #   if(sum(dat_plot_rb_i[,t]==0)==0){
    #     area_i <- c(area_i,rownames(dat_plot_rb_i)[nrow(dat_plot_rb_i)])
    #   } else if (which(dat_plot_rb_i[,t]==0)[1]==1){
    #     area_i <- c(area_i,0)
    #   } else {      
    #     area_i <- c(area_i,rownames(dat_plot_rb_i)[which(dat_plot_rb_i[,t]==0)[1]-1])
    #   }
    # }

    #print(peak_height_i)

    # # Normalize peak height vs max peak height (primer CpG incl)
    # meth_state[sample_j_ori,amplicon] <- max(as.numeric(peak_height_i))/sum(db_rb_anal[[amplicon]][[sample_j_ori]]$count)
    # Normalize peak height vs max peak height (primer CpG non incl)
    # meth_state[sample_j_ori,amplicon] <- max(as.numeric(peak_height_i))/as.numeric(rownames(dat_plot_rb_i)[rowSums(dat_plot_rb_i)!=0][sum(rowSums(dat_plot_rb_i)!=0)])
    
    # Normalize area vs max area
    # zone_oi <- dat_plot_rb_i[1:which(rownames(dat_plot_rb_i)==max(as.numeric(area_i))),]
    # zone_oi_incr <- as.numeric(c(rownames(zone_oi)[1],as.numeric(rownames(zone_oi)[2:nrow(zone_oi)])-as.numeric(rownames(zone_oi)[1:(nrow(zone_oi)-1)])))
    # meth_state2[sample_j_ori,amplicon] <- sum(rowSums(zone_oi)*zone_oi_incr)/sum(ncol(zone_oi)*zone_oi_incr)

    # new measures Tim
    # meth_state[sample_j_ori,amplicon] <- as.numeric(tail(rownames(dat_plot_rb_i[rowSums(dat_plot_rb_i)!=0,]),1))/as.numeric(tail(rownames(dat_plot_rb_i),1))
 
    # zone_oi <- dat_plot_rb_i[rowSums(dat_plot_rb_i)!=0,]
    # zone_oi_incr <- as.numeric(c(rownames(zone_oi)[1],as.numeric(rownames(zone_oi)[2:nrow(zone_oi)])-as.numeric(rownames(zone_oi)[1:(nrow(zone_oi)-1)])))
    # meth_state2[sample_j_ori,amplicon] <- sum(rowSums(zone_oi)*zone_oi_incr)/sum(ncol(zone_oi)*zone_oi_incr)


# ## Methrate plots
# for (amplicon in amplicons){
#   methrate <- NULL
#   print(sprintf("Plotting for %s",amplicon))
#   for (sample_j_ori in names(db_rb_anal[[amplicon]])){
#     dat_sample_j <- db_rb_anal[[amplicon]][[sample_j_ori]]
    
#     if(grepl("LY75",amplicon)){
#       betas_ampl <- betas_LY75[grepl(amplicon,rownames(betas_LY75)),]
#     } else if (grepl("GREM1",amplicon)){
#       betas_ampl <- betas_GREM1[grepl(amplicon,rownames(betas_GREM1)),]
#     } else if (grepl("NDRG4",amplicon)){
#       betas_ampl <- betas_NDRG4[grepl(amplicon,rownames(betas_NDRG4)),]
#     }
#     CpG_in_primer_ampl <- which(rownames(betas_ampl)%in%CpG_in_primer)

#     # Check rb dna meth counts vs biseq
#     #sample_annot[sample_annot$Sample_ID==sample_j_ori,]
#     #betas_ampl[,"P12"]
#     #sum(dat_sample_j$count[grepl("^C.*$",dat_sample_j$allele)])/sum(dat_sample_j$count)

#     for (i in 1:nrow(dat_sample_j)){
#       if (length(CpG_in_primer_ampl)>0){
#         methrate <-  c(methrate,rep(sum(grepl("C",unlist(strsplit(dat_sample_j$allele[i],""))[-CpG_in_primer_ampl]))/(nchar(dat_sample_j$allele[i])-length(CpG_in_primer_ampl)),dat_sample_j$count[i]))
#       } else {
#         methrate <-  c(methrate,rep(sum(grepl("C",unlist(strsplit(dat_sample_j$allele[i],""))))/nchar(dat_sample_j$allele[i]),dat_sample_j$count[i]))
#       }
#     }
#   }
#   if (!is.null(methrate)){
#     jpeg(paste0("waterfall_plots/methrate_",amplicon,".jpg"),width=512,height=512)
#     hist(methrate)
#     dev.off()
#   }
# }

# ## Meth_perc plots
# meth_state <- data.frame(matrix(rep(NA,length(amplicons)*nrow(sample_annot)),nrow=nrow(sample_annot)))
# colnames(meth_state) <- amplicons
# rownames(meth_state) <- sample_annot$Sample_ID
# meth_state2 <- meth_state
# cutoff_mono <- 0.5
# for (amplicon in amplicons){
#   print(sprintf("Plotting for %s",amplicon))
#   dat_methrate_ps <- NULL
#   for (sample_j_ori in names(db_rb_anal[[amplicon]])){
#     methrate <- NULL
#     dat_sample_j <- db_rb_anal[[amplicon]][[sample_j_ori]]
#     for (i in 1:nrow(dat_sample_j)){
#         methrate <-  c(methrate,rep(sum(grepl("C",unlist(strsplit(dat_sample_j$allele[i],""))))/nchar(dat_sample_j$allele[i]),dat_sample_j$count[i]))
#     }
#     # if (grepl("GREM1",amplicon)){
#     #       dat_methrate_ps <- rbind(dat_methrate_ps,c(sample_j_ori,
#     #                                            sum(methrate > cutoff_mono)/length(methrate),
#     #                                            sum(methrate < 0.25)/sum(methrate > 0.25 & methrate < 0.75)))
#     # } else {
#       dat_methrate_ps <- rbind(dat_methrate_ps,c(sample_j_ori,
#                                                 sum(methrate > cutoff_mono)/length(methrate),
#                                                 (sum(methrate > 0.75)+1)/(sum(methrate > 0.25 & methrate < 0.75)+1)))
#     # }
#   }
#   if (!is.null(dat_methrate_ps)){
#     jpeg(paste0("waterfall_plots/perc_meth_",amplicon,".jpg"),width=512,height=512)
#     hist(as.numeric(dat_methrate_ps[,2]),breaks=20)
#     dev.off()
#   }
#   meth_state[rownames(meth_state)%in%dat_methrate_ps[,1],amplicon] <- as.numeric(dat_methrate_ps[,2])
#   meth_state2[rownames(meth_state2)%in%dat_methrate_ps[,1],amplicon] <- as.numeric(dat_methrate_ps[,3])
#   #meth_state[rownames(meth_state)%in%dat_methrate_ps[,1],amplicon] <- c("non","mono","bi")[c(as.numeric(dat_methrate_ps[,2])>=0.1) + c(as.numeric(dat_methrate_ps[,2])>=0.6)+1]
# }

# amplicons_NDRG4 <- AIC_new[AIC_new[,2]>0,1]
# amplicons_NDRG4 <- amplicons[grepl("NDRG4",amplicons)]
# for (amplicon in amplicons_NDRG4){
#   CRC_annot_sort$meth_state <- NA
#   for (i in 1:nrow(CRC_annot_sort)){
#     CRC_annot_sort$meth_state[i] <- meth_state[sample_annot$Sample_ID[sample_annot$CRC==CRC_annot_sort$Sample_ID[i]],amplicon]
#   }
#   #CRC_annot_sort$meth_state <- c("yes","no")[c(CRC_annot_sort$meth_state>0.05) + 1]
#   formula_new <- formula(paste("Tissue ~ meth_state+Gender+Age+Patient_ID",sep=""))
#   cn_formula_new <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_new))),"+",fixed=TRUE)))
#   model_new <- glm(formula_new, 
#                   data = CRC_annot_sort, 
#                   family = binomial(link="logit"))
#   #print(summary(model_new))
#   AIC_new[AIC_new[,1]==amplicon,2] <- amplicon,AIC(model_new) - AIC(model_CF)
# }
# AIC_new

# AIC_new <- NULL
# amplicons_NDRG4 <- amplicons[grepl("NDRG4",amplicons)]
# for (amplicon in amplicons_NDRG4){
#   CRC_annot_sort$meth_state <- NA
#   for (i in 1:nrow(CRC_annot_sort)){
#     CRC_annot_sort$meth_state[i] <- meth_state[sample_annot$Sample_ID[sample_annot$CRC==CRC_annot_sort$Sample_ID[i]],amplicon]
#   }
#   #CRC_annot_sort$meth_state <- c("yes","no")[c(CRC_annot_sort$meth_state>0.05) + 1]
#   formula_new <- formula(paste("Tissue ~ meth_state+Gender+Age+Patient_ID",sep=""))
#   cn_formula_new <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_new))),"+",fixed=TRUE)))
#   model_new <- glm(formula_new, 
#                   data = CRC_annot_sort, 
#                   family = binomial(link="logit"))
#   #print(summary(model_new))
#   if (AIC(model_new) - AIC(model_CF)<0){
#     AIC_new <- rbind(AIC_new,c(amplicon,AIC(model_new) - AIC(model_CF)))
#   } else {
#     CRC_annot_sort$meth_state <- c("yes","no")[c(CRC_annot_sort$meth_state>0.05) + 1]
#     formula_new <- formula(paste("Tissue ~ meth_state+Gender+Age+Patient_ID",sep=""))
#     cn_formula_new <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_new))),"+",fixed=TRUE)))
#     model_new <- glm(formula_new, 
#                     data = CRC_annot_sort, 
#                     family = binomial(link="logit"))
#     AIC_new <- rbind(AIC_new,c(amplicon,AIC(model_new) - AIC(model_CF)))
#   }
# }
# amplicons_not_NDRG4 <- amplicons[!grepl("NDRG4",amplicons)]
# for (amplicon in amplicons_not_NDRG4){
#   CRC_annot_sort$meth_state <- NA
#   for (i in 1:nrow(CRC_annot_sort)){
#     CRC_annot_sort$meth_state[i] <- meth_state[sample_annot$Sample_ID[sample_annot$CRC==CRC_annot_sort$Sample_ID[i]],amplicon]
#   }
#   CRC_annot_sort$meth_state <- c("yes","no")[c(CRC_annot_sort$meth_state>0.05) + 1]
#   formula_new <- formula(paste("Tissue ~ meth_state+Gender+Age+Patient_ID",sep=""))
#   cn_formula_new <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_new))),"+",fixed=TRUE)))
#   model_new <- glm(formula_new, 
#                   data = CRC_annot_sort, 
#                   family = binomial(link="logit"))
#   #print(summary(model_new))
#   AIC_new <- rbind(AIC_new,c(amplicon,AIC(model_new) - AIC(model_CF)))
# }

#####################################################
# AIC analyses before implementing it as a function #
#####################################################

# #NDRG4
# set.seed(seed_perm)
# pval_list_NDRG4 <- list()
# for (k in 1:nrow(primer_annot_NDRG4)){
#   print(paste0("Primer: ",primer_annot_NDRG4$name[k]))
#   print(paste0("Location: chr ",primer_annot_NDRG4$chromosome[k]," from ",primer_annot_NDRG4$begin[k]," to ",primer_annot_NDRG4$end[k]))
#   mask_amplicon <- grepl(paste0(primer_annot_NDRG4$name[k],"$"),rownames(Mval_NDRG4_sort))
#   dat_amplicon <- t(Mval_NDRG4_sort[mask_amplicon,])
  
#   dim(dat_amplicon)
#   sum(CRC_annot_sort$Sample_ID==rownames(dat_amplicon))
  
#   dat_amplicon <- cbind(CRC_annot_sort[,colnames(CRC_annot_sort)%in%c("Tissue")],dat_amplicon)
#   colnames(dat_amplicon) <- c("Tissue", paste0(rep("CpG_",sum(mask_amplicon)),rownames(Mval_NDRG4_sort[mask_amplicon,])))
  
#   cn_names <- colnames(dat_amplicon)[grepl("CpG",colnames(dat_amplicon))]
#   cn_names <- gsub("CpG_","",cn_names)
  
#   print(paste0(length(cn_names)," CpGs"))
  
#   if (length(cn_names)==0){
#     pval_list_NDRG4[[k]] <- NULL
#     next()
#   }
  
#   dat_amplicon <- data.frame(dat_amplicon,stringsAsFactors = F)
#   dat_amplicon$Tissue <- as.factor(dat_amplicon$Tissue)
  
#   head(dat_amplicon)
  
#   AIC_amplicon <- NULL
#   for (i in 1:length(cn_names)){
#     dat_test_i <- dat_amplicon[,c("Tissue",paste0("CpG_",cn_names[i]))]
    
#     CpG_name_i <- colnames(dat_test_i)[grepl("CpG",colnames(dat_test_i))]
#     CpG_name_i <- gsub("^CpG_","",CpG_name_i)
    
#     formula_i <- formula(paste("Tissue ~ CpG_",CpG_name_i,sep=""))
#     cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
#     logreg_BS_i <- glm(formula_i, data = dat_test_i, family = binomial(link="logit"))
    
#     AIC_i <- AIC(logreg_BS_i)
#     AIC_amplicon <- rbind(AIC_amplicon,c(CpG_name_i,AIC_i))
#   }
#   AIC_amplicon <- data.frame(AIC_amplicon, stringsAsFactors = F)
#   colnames(AIC_amplicon) <- c("CpG","AIC")
#   AIC_amplicon$AIC <- as.numeric(AIC_amplicon$AIC)
    
#   combos <- ncol(combn(1:length(cn_names),2))
#   pval <- data.frame(CpG_1=rep("",combos),
#                      CpG_2=rep("",combos),
#                      AIC_1=rep("",combos),
#                      AIC_2=rep("",combos),
#                      delta_AIC=rep("",combos),
#                      pval=rep(1,combos),
#                      stringsAsFactors = F)
#   counter <- 1
#   for (i in 1:(length(cn_names)-1)){
#     for (j in (i+1):length(cn_names)){
#       print(paste0("CpG comparison: ",paste(cn_names[c(i,j)],collapse=" vs ")))
      
#       pval[counter,1:2] <- cn_names[c(i,j)]
      
#       dat_test_i <- dat_amplicon[,c("Tissue",paste0(rep("CpG_",2),cn_names[c(i,j)]))]
#       pval[counter,3:5] <- calculate_delta_AIC_diagnostic(dat_test_i)

#       distr_AIC <- mclapply(1:10000,bootstrap_AIC_diagnostic,mc.cores=40)
#       distr_AIC <- unlist(distr_AIC)

#       jpeg(paste0("AIC_plots/CpG_plots/CRC_CpG_",cn_names[i],"vs",cn_names[j],".jpg"))
#       hist(distr_AIC,breaks=50)
#       abline(v=pval[counter,5],col="red")
#       dev.off()
      
#       pval[counter,6] <- 2*min(sum(as.numeric(pval[counter,5])>distr_AIC)/length(distr_AIC),sum(as.numeric(pval[counter,5])<distr_AIC)/length(distr_AIC))
#       counter <- counter + 1
#     }
#   } 
#   pval$fdr <- p.adjust(pval$pval)
#   pval_list_NDRG4[[k]] <- pval 
  
#   if (sum(pval$fdr<fdr_treshold_toplot)>0){
#     pval_sign <- pval[pval$fdr<fdr_treshold_toplot,]
#     print(pval_sign)
#     if (fdr_treshold_toplot%%100/10==5){
#       sign_level <- (pval_sign$fdr < 0.05) + (pval_sign$fdr < 0.005) + (pval_sign$fdr < 0.0005) + (pval_sign$fdr < 0.0005)
#     } else {
#       sign_level <- (pval_sign$fdr < 0.1) + (pval_sign$fdr < 0.01) + (pval_sign$fdr < 0.001) + (pval_sign$fdr < 0.001)
#     }
    
#     pos_x <-as.numeric(gsub("^.+_","",gsub("_[NDRG|LY|GREM].+$","",AIC_amplicon$CpG)))
    
#     jpeg(paste0("AIC_plots/",primer_annot_NDRG4$name[k],".jpg"))
#     plot(pos_x,AIC_amplicon$AIC,type="l") # this line is just to get limits for plots in next command
#     add_space <- (par("usr")[4]-par("usr")[3])*0.05*(nrow(pval_sign)+1)
#     # print(add_space)
#     plot(pos_x,AIC_amplicon$AIC,type="l",main=primer_annot_NDRG4$name[k],
#          ylim=c((par("usr")[3]-add_space),par("usr")[4]),xlab="Pos",ylab="AIC",cex.lab=cex_AIC_plots, cex.axis=cex_AIC_plots)
#     increment_space <- add_space/(nrow(pval_sign)+1)
#     for (t in 1:nrow(pval_sign)){
#       segments(pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])],par("usr")[3]+increment_space*t,
#                pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])],par("usr")[3]+increment_space*t)
#       mid <- (pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])] +
#                 pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])])/2
#       text(mid,par("usr")[3]+increment_space*t-(0.33*increment_space),paste(rep("*",sign_level[t]),collapse=""))
#     }
#     points(pos_x,AIC_amplicon$AIC,pch=16)
#     n_bp_add <- 2
#     rect(primer_annot_list[["NDRG4"]]$begin[k]-n_bp_add,(par("usr")[3]-add_space),
#          primer_annot_list[["NDRG4"]]$begin[k]+primer_annot_list[["NDRG4"]]$len_fwd[k]+n_bp_add,par("usr")[4],
#          col=alpha("cadetblue",alpha=0.20))
#     rect(primer_annot_list[["NDRG4"]]$end[k]-primer_annot_list[["NDRG4"]]$len_rev[k]-n_bp_add,(par("usr")[3]-add_space),
#          primer_annot_list[["NDRG4"]]$end[k]+n_bp_add,par("usr")[4],
#          col=alpha("cadetblue",alpha=0.20))
#     dev.off()
#   } else {
#     print("No significant differences!")
#   }
# }
# names(pval_list_NDRG4) <- primer_annot_NDRG4$name
# saveRDS(pval_list_NDRG4,paste0("pval_list_NDRG4_seed",seed_perm,".rds"))

# #GREM1
# set.seed(seed_perm)
# pval_list_GREM1 <- list()
# for (k in 1:nrow(primer_annot_GREM1)){
#   print(paste0("Primer: ",primer_annot_GREM1$name[k]))
#   print(paste0("Location: chr ",primer_annot_GREM1$chromosome[k]," from ",primer_annot_GREM1$begin[k]," to ",primer_annot_GREM1$end[k]))
#   mask_amplicon <- grepl(paste0(primer_annot_GREM1$name[k],"$"),rownames(Mval_GREM1_sort))
#   dat_amplicon <- t(Mval_GREM1_sort[mask_amplicon,])
  
#   dim(dat_amplicon)
#   sum(ccRCC_annot_sort$Sample_ID==rownames(dat_amplicon))
  
#   dat_amplicon <- cbind(ccRCC_annot_sort[,colnames(ccRCC_annot_sort)%in%c("Followup_time","NK_death")],dat_amplicon)
#   colnames(dat_amplicon) <- c(colnames(ccRCC_annot_sort)[colnames(ccRCC_annot_sort)%in%c("Followup_time","NK_death")],
#                               paste0(rep("CpG_",sum(mask_amplicon)),rownames(Mval_GREM1_sort[mask_amplicon,])))
  
#   cn_names <- colnames(dat_amplicon)[grepl("CpG",colnames(dat_amplicon))]
#   cn_names <- gsub("CpG_","",cn_names)
  
#   print(paste0(length(cn_names)," CpGs"))
  
#   if (length(cn_names)==0){
#     pval_list_GREM1[[k]] <- "No CpGs"
#     next()
#   }
  
#   dat_amplicon <- data.frame(dat_amplicon,stringsAsFactors = F)
  
#   head(dat_amplicon)
  
#   AIC_amplicon <- NULL
#   for (i in 1:length(cn_names)){
#     dat_test_i <- dat_amplicon[,c("Followup_time","NK_death",paste0("CpG_",cn_names[i]))]
    
#     CpG_name_i <- colnames(dat_test_i)[grepl("CpG",colnames(dat_test_i))]
#     CpG_name_i <- gsub("^CpG_","",CpG_name_i)
    
#     formula_i <- formula(paste("Surv(Followup_time, NK_death) ~ CpG_",CpG_name_i,sep=""))
#     cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
#     cox_BS_i <- coxph(formula_i, data = dat_test_i)

#     AIC_i <- AIC(cox_BS_i)
#     AIC_amplicon <- rbind(AIC_amplicon,c(CpG_name_i,AIC_i))
#   }
#   AIC_amplicon <- data.frame(AIC_amplicon, stringsAsFactors = F)
#   colnames(AIC_amplicon) <- c("CpG","AIC")
#   AIC_amplicon$AIC <- as.numeric(AIC_amplicon$AIC)
  
  
#   combos <- ncol(combn(1:length(cn_names),2))
#   pval <- data.frame(CpG_1=rep("",combos),
#                      CpG_2=rep("",combos),
#                      AIC_1=rep("",combos),
#                      AIC_2=rep("",combos),
#                      delta_AIC=rep("",combos),
#                      pval=rep(1,combos),
#                      stringsAsFactors = F)
#   counter <- 1
#   for (i in 1:(length(cn_names)-1)){
#     print(sprintf("Comparing %s with %s CpGs",cn_names[i],length((i+1):length(cn_names))))
#     for (j in (i+1):length(cn_names)){
#       #print(paste0("CpG comparison: ",paste(cn_names[c(i,j)],collapse=" vs ")))
      
#       pval[counter,1:2] <- cn_names[c(i,j)]
      
#       dat_test_i <- dat_amplicon[,c("Followup_time","NK_death",paste0(rep("CpG_",2),cn_names[c(i,j)]))]
#       pval[counter,3:5] <- calculate_delta_AIC_prognostic(dat_test_i)
      
#       distr_AIC <- mclapply(1:10000,bootstrap_AIC_prognostic,mc.cores=40)
#       distr_AIC <- unlist(distr_AIC)

#       jpeg(paste0("AIC_plots/CpG_plots/ccRCC_CpG_",cn_names[i],"vs",cn_names[j],".jpg"))
#       hist(distr_AIC,breaks=50)
#       abline(v=pval[counter,5],col="red")
#       dev.off()
#       pval[counter,6] <- 2*min(sum(as.numeric(pval[counter,5])>distr_AIC)/length(distr_AIC),sum(as.numeric(pval[counter,5])<distr_AIC)/length(distr_AIC))
#       counter <- counter + 1
#     }
#   }
#   pval$fdr <- p.adjust(pval$pval)
#   pval_list_GREM1[[k]] <- pval 
  
#   if (sum(pval$fdr<fdr_treshold_toplot)>0){
#     pval_sign <- pval[pval$fdr<fdr_treshold_toplot,]
#     print(pval_sign)
#     if (fdr_treshold_toplot%%100/10==5){
#       sign_level <- (pval_sign$fdr < 0.05) + (pval_sign$fdr < 0.005) + (pval_sign$fdr < 0.0005) + (pval_sign$fdr < 0.0005)
#     } else {
#       sign_level <- (pval_sign$fdr < 0.1) + (pval_sign$fdr < 0.01) + (pval_sign$fdr < 0.001) + (pval_sign$fdr < 0.001)
#     }
    
#     pos_x <-as.numeric(gsub("^.+_","",gsub("_[NDRG|LY|GREM].+$","",AIC_amplicon$CpG)))
    
#     jpeg(paste0("AIC_plots/",primer_annot_GREM1$name[k],".jpg"))
#     plot(pos_x,AIC_amplicon$AIC,type="l") # this line is just to get limits for plots in next command
#     add_space <- (par("usr")[4]-par("usr")[3])*0.05*(nrow(pval_sign)+1)
#     # print(add_space)
#     plot(pos_x,AIC_amplicon$AIC,type="l",main=primer_annot_GREM1$name[k],
#          ylim=c((par("usr")[3]-add_space),par("usr")[4]),xlab="Pos",ylab="AIC",cex.lab=cex_AIC_plots,cex.axis=cex_AIC_plots)
#     increment_space <- add_space/(nrow(pval_sign)+1)
#     for (t in 1:nrow(pval_sign)){
#       segments(pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])],par("usr")[3]+increment_space*t,
#                pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])],par("usr")[3]+increment_space*t)
#       mid <- (pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])] +
#                 pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])])/2
#       text(mid,par("usr")[3]+increment_space*t-(0.33*increment_space),paste(rep("*",sign_level[t]),collapse=""))
#     }
#     points(pos_x,AIC_amplicon$AIC,pch=16)
#     dev.off()
#   } else {
#     print("No significant differences!")
#   }
# }
# names(pval_list_GREM1) <- primer_annot_GREM1$name
# saveRDS(pval_list_GREM1,paste0("pval_list_GREM1_seed",seed_perm,".rds"))

# #LY75
# colnames(melanoma_annot_sort)[colnames(melanoma_annot_sort)=="MK_death"] <- "NK_death"
# set.seed(seed_perm)
# pval_list_LY75 <- list()
# for (k in 1:nrow(primer_annot_LY75)){
#   print(paste0("Primer: ",primer_annot_LY75$name[k]))
#   print(paste0("Location: chr ",primer_annot_LY75$chromosome[k]," from ",primer_annot_LY75$begin[k]," to ",primer_annot_LY75$end[k]))
#   mask_amplicon <- grepl(paste0(primer_annot_LY75$name[k],"$"),rownames(Mval_LY75_sort))
#   dat_amplicon <- t(Mval_LY75_sort[mask_amplicon,])
  
#   dim(dat_amplicon)
#   sum(melanoma_annot_sort$Sample_ID==rownames(dat_amplicon))
  
#   dat_amplicon <- cbind(melanoma_annot_sort[,colnames(melanoma_annot_sort)%in%c("Followup_time","NK_death")],dat_amplicon)
#   colnames(dat_amplicon) <- c(colnames(melanoma_annot_sort)[colnames(melanoma_annot_sort)%in%c("Followup_time","NK_death")],
#                               paste0(rep("CpG_",sum(mask_amplicon)),rownames(Mval_LY75_sort[mask_amplicon,])))
  
#   cn_names <- colnames(dat_amplicon)[grepl("CpG",colnames(dat_amplicon))]
#   cn_names <- gsub("CpG_","",cn_names)
  
#   print(paste0(length(cn_names)," CpGs"))
  
#   if (length(cn_names)==0){
#     pval_list_LY75[[k]] <- "No CpGs"
#     next()
#   }
  
#   dat_amplicon <- data.frame(dat_amplicon,stringsAsFactors = F)

#   head(dat_amplicon)
  
#   AIC_amplicon <- NULL
#   for (i in 1:length(cn_names)){
#     dat_test_i <- dat_amplicon[,c("Followup_time","NK_death",paste0("CpG_",cn_names[i]))]
    
#     CpG_name_i <- colnames(dat_test_i)[grepl("CpG",colnames(dat_test_i))]
#     CpG_name_i <- gsub("^CpG_","",CpG_name_i)
    
#     formula_i <- formula(paste("Surv(Followup_time, NK_death) ~ CpG_",CpG_name_i,sep=""))
#     cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
#     cox_BS_i <- coxph(formula_i, data = dat_test_i)
    
#     AIC_i <- AIC(cox_BS_i)
#     AIC_amplicon <- rbind(AIC_amplicon,c(CpG_name_i,AIC_i))
#   }
#   AIC_amplicon <- data.frame(AIC_amplicon, stringsAsFactors = F)
#   colnames(AIC_amplicon) <- c("CpG","AIC")
#   AIC_amplicon$AIC <- as.numeric(AIC_amplicon$AIC)
  
  
#   combos <- ncol(combn(1:length(cn_names),2))
#   pval <- data.frame(CpG_1=rep("",combos),
#                      CpG_2=rep("",combos),
#                      AIC_1=rep("",combos),
#                      AIC_2=rep("",combos),
#                      delta_AIC=rep("",combos),
#                      pval=rep(1,combos),
#                      stringsAsFactors = F)
#   counter <- 1
#   for (i in 1:(length(cn_names)-1)){
#     for (j in (i+1):length(cn_names)){
#       print(paste0("CpG comparison: ",paste(cn_names[c(i,j)],collapse=" vs ")))

#       pval[counter,1:2] <- cn_names[c(i,j)]
      
#       dat_test_i <- dat_amplicon[,c("Followup_time","NK_death",paste0(rep("CpG_",2),cn_names[c(i,j)]))]
#       pval[counter,3:5] <- calculate_delta_AIC_prognostic(dat_test_i)
#       distr_AIC <- rep(0,1000)
      
#       distr_AIC <- mclapply(1:10000,bootstrap_AIC_prognostic,mc.cores=40)
#       distr_AIC <- unlist(distr_AIC)

#       jpeg(paste0("AIC_plots/CpG_plots/melanoma_CpG_",cn_names[i],"vs",cn_names[j],".jpg"))
#       hist(distr_AIC,breaks=50)
#       abline(v=pval[counter,5],col="red")
#       dev.off()
#       pval[counter,6] <- 2*min(sum(as.numeric(pval[counter,5])>distr_AIC)/length(distr_AIC),sum(as.numeric(pval[counter,5])<distr_AIC)/length(distr_AIC))
#       counter <- counter + 1
#     }
#   }
#   pval$fdr <- p.adjust(pval$pval)
#   pval_list_LY75[[k]] <- pval 
  
#   if (sum(pval$fdr<fdr_treshold_toplot)>0){
#     pval_sign <- pval[pval$fdr<fdr_treshold_toplot,]
#     print(pval_sign)
#     if (fdr_treshold_toplot%%100/10==5){
#       sign_level <- (pval_sign$fdr < 0.05) + (pval_sign$fdr < 0.005) + (pval_sign$fdr < 0.0005) + (pval_sign$fdr < 0.0005)
#     } else {
#       sign_level <- (pval_sign$fdr < 0.1) + (pval_sign$fdr < 0.01) + (pval_sign$fdr < 0.001) + (pval_sign$fdr < 0.001)
#     }
    
#     pos_x <-as.numeric(gsub("^.+_","",gsub("_[NDRG|LY|GREM].+$","",AIC_amplicon$CpG)))
    
#     jpeg(paste0("AIC_plots/",primer_annot_LY75$name[k],".jpg"))
#     plot(pos_x,AIC_amplicon$AIC,type="l") # this line is just to get limits for plots in next command
#     add_space <- (par("usr")[4]-par("usr")[3])*0.05*(nrow(pval_sign)+1)
#     # print(add_space)
#     plot(pos_x,AIC_amplicon$AIC,type="l",main=primer_annot_LY75$name[k],
#          ylim=c((par("usr")[3]-add_space),par("usr")[4]),xlab="Pos",ylab="AIC",cex.lab=cex_AIC_plots,cex.axis=cex_AIC_plots)
#     increment_space <- add_space/(nrow(pval_sign)+1)
#     for (t in 1:nrow(pval_sign)){
#       segments(pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])],par("usr")[3]+increment_space*t,
#                pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])],par("usr")[3]+increment_space*t)
#       mid <- (pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])] +
#                 pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])])/2
#       text(mid,par("usr")[3]+increment_space*t-(0.33*increment_space),paste(rep("*",sign_level[t]),collapse=""))
#     }
#     points(pos_x,AIC_amplicon$AIC,pch=16)
#     dev.off()
#   } else {
#     print("No significant differences!")
#   }
# }
# names(pval_list_LY75) <- primer_annot_LY75$name
# saveRDS(pval_list_LY75,paste0("pval_list_LY75_seed",seed_perm,".rds"))
