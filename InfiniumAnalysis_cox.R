#################
#################
##             ##
## Meth & CNV  ##
##             ##
#################
#################

# Get cancer type
cancertype <- commandArgs(TRUE)
print(cancertype)

if(identical(cancertype,character(0))){
  cancertype <- "KIRC"
  warning("Cancertype set as KIRC due to missing specification")
}

# Load libraries
library(parallel)
library(survival)
library(biomaRt)
library(survminer)
library(gridExtra)
library(ggplot2)

# Sample annotation
sample_annot <- read.table("gdc_sample_sheet.2021-02-23.tsv",header=T,sep="\t",stringsAsFactors=F)
sample_annot <- sample_annot[sample_annot$Project.ID==paste("TCGA-",cancertype,sep=""),]

control_annot <- sample_annot$Sample.ID[grep("Normal",sample_annot$Sample.Type)]
control_annot <- gsub("-","_",control_annot)
control_annot <- gsub("A$","",control_annot)

tumor_annot <- sample_annot$Sample.ID[grep("Tumor",sample_annot$Sample.Type)]
tumor_annot <- gsub("-","_",tumor_annot)
tumor_annot <- gsub("A$","",tumor_annot)
tumor_annot <- gsub("B$","",tumor_annot)
# Remove duplicate samples
tumor_annot <- unique(tumor_annot)

# Methylation annotation
meth450_annot_grch38 <- read.table("EPIC.anno.GRCh38.tsv",header=T,sep="\t",stringsAsFactors = F)
dim(meth450_annot_grch38)
meth450_annot_grch38 <- meth450_annot_grch38[sort(meth450_annot_grch38$probeID, index.return=T)$ix,]

# Methylation Data
dat_meth450 <- read.table(paste(cancertype,"_HumanMethylation450",sep=""),header=T,stringsAsFactors= F)
colnames(dat_meth450) <- gsub("\\.","\\_",colnames(dat_meth450))

dat_meth450_rn <- dat_meth450[,1]
dat_meth450 <- dat_meth450[,-1]

print(dim(dat_meth450))

dat_meth450 <- data.frame(dat_meth450)

dat_meth450_tumor <- dat_meth450[,colnames(dat_meth450)%in%tumor_annot]
dat_meth450_tumor <- dat_meth450_tumor[,sort(colnames(dat_meth450_tumor),index.return=T)$ix]
dat_meth450_tumor <- sapply(dat_meth450_tumor,as.numeric)
rownames(dat_meth450_tumor) <- dat_meth450_rn

dat_meth450_control <- dat_meth450[,colnames(dat_meth450)%in%control_annot]
dat_meth450_control <- dat_meth450_control[,sort(colnames(dat_meth450_control),index.return=T)$ix]
dat_meth450_control <- sapply(dat_meth450_control,as.numeric)
rownames(dat_meth450_control) <- dat_meth450_rn

print(dim(dat_meth450_tumor))
print(head(dat_meth450_tumor[,1:10]))

print(dim(dat_meth450_control))
print(head(dat_meth450_control[,1:10]))

# Only samples with expression data => DO NOT DO THIS
# samples_oi <- sort(header_expression[(header_expression%in%header_meth450) & (!(header_expression%in%control_annot))])
# dat_expression <- dat_expression[,colnames(dat_expression)%in%samples_oi]
# dat_meth450 <- dat_meth450[,colnames(dat_meth450)%in%samples_oi]
# 

# 
# if ((sum(colnames(dat_expression)==samples_oi)!=length(samples_oi)) |
#     (sum(colnames(dat_meth450)==samples_oi)!=length(samples_oi))){
#   stop("Problem with samples!")
# }

# which(!(c("age_at_initial_pathologic_diagnosis","days_to_birth","days_to_death",
#         "days_to_last_followup","gender","histological_type","pathologic_stage")%in%colnames(pat_info)))

pat_info <- read.table(paste(cancertype,"_clinicalMatrix",sep=""),header=T,row.names = 1,sep = "\t")
pat_info <- pat_info[,c("age_at_initial_pathologic_diagnosis","days_to_birth","days_to_death",
                       "days_to_last_followup","gender","histological_type", "pathologic_stage")]
rownames(pat_info) <- gsub("\\-","_",toupper(rownames(pat_info)))
colnames(pat_info) <- c("Age","DaysToBirth","DaysToDeath","DaysToLastFollowup","Gender",
                        "HistType","Stage")

pat_info <- data.frame(pat_info,stringsAsFactors = F)
pat_info$Age <- as.numeric(pat_info$Age)
pat_info$DaysToBirth <- as.numeric(pat_info$DaysToBirth)
pat_info$DaysToDeath <- as.numeric(pat_info$DaysToDeath)
pat_info$DaysToLastFollowup <- as.numeric(pat_info$DaysToLastFollowup)
pat_info$HistType <- as.factor(pat_info$HistType)
pat_info$Stage[!grepl("Stage",pat_info$Stage)] <-NA
pat_info$Stage <- as.factor(as.character(pat_info$Stage))

pat_info <- pat_info[sort(rownames(pat_info),index.return=T)$ix,]

pat_info_meth450 <- pat_info[rownames(pat_info)%in%colnames(dat_meth450_tumor),]

dim(pat_info_meth450)
dim(dat_meth450_tumor)

if (sum(rownames(pat_info_meth450)==colnames(dat_meth450_tumor))!=nrow(pat_info_meth450)){
  stop("Dimensions wrong!")
}

pat_info_meth450$event <- as.numeric(!is.na(pat_info_meth450$DaysToDeath))
pat_info_meth450$time <- rep("",nrow(pat_info_meth450))
pat_info_meth450$time[pat_info_meth450$event==1] <- pat_info_meth450$DaysToDeath[pat_info_meth450$event==1]
pat_info_meth450$time[pat_info_meth450$event==0] <- pat_info_meth450$DaysToLastFollowup[pat_info_meth450$event==0]
pat_info_meth450$time <- as.numeric(pat_info_meth450$time)

print(head(pat_info_meth450))

cte <- 0.01
dat_meth450_tumor_M <- log2((dat_meth450_tumor+cte)/(1-dat_meth450_tumor+cte))
head(dat_meth450_tumor_M[,1:6])
rownames(dat_meth450_tumor_M) <- rownames(dat_meth450_tumor)

pat_info_meth450$Survival <- Surv(pat_info_meth450$time,pat_info_meth450$event)

dat_meth450_tumor_filt <- dat_meth450_tumor[rowSums(is.na(dat_meth450_tumor_M))<ceiling(ncol(dat_meth450_tumor_M)),]
dat_meth450_tumor_M_filt <- dat_meth450_tumor_M[rowSums(is.na(dat_meth450_tumor_M))<ceiling(ncol(dat_meth450_tumor_M)),]

# Model Survival ifo Age, Gender and Stage
summary(pat_info_meth450)


# Basic barplot
df_plot <- data.frame(Gender=as.factor(unique(pat_info_meth450$Gender)),
                 Frequency=c(sum(pat_info_meth450$Gender==as.factor(unique(pat_info_meth450$Gender))[1]),
                             sum(pat_info_meth450$Gender==as.factor(unique(pat_info_meth450$Gender))[2])),
                 stringsAsFactors = F)
p <- ggplot(data=df_plot, aes(x=Gender,y=Frequency,fill=Gender)) +
  geom_bar(stat="identity") + 
  scale_fill_brewer(palette="Dark2") + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,250)) +
  coord_flip() + 
  geom_text(aes(label=Frequency),nudge_y=7,size=3.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.y = element_blank(),
        axis.line.x = element_line(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Gender distribution") 
ggsave(paste(cancertype,"_Gender_distr.jpg",sep = ""),p, height = 2.5, width = 9)

df_plot <- data.frame(Stage=as.factor(unique(pat_info_meth450$Stage)),
                      Frequency=c(sum((pat_info_meth450$Stage==as.factor(unique(pat_info_meth450$Stage))[1]) & !is.na(pat_info_meth450$Stage)),
                                  sum((pat_info_meth450$Stage==as.factor(unique(pat_info_meth450$Stage))[2]) & !is.na(pat_info_meth450$Stage)),
                                  sum((pat_info_meth450$Stage==as.factor(unique(pat_info_meth450$Stage))[3]) & !is.na(pat_info_meth450$Stage)),
                                  sum((pat_info_meth450$Stage==as.factor(unique(pat_info_meth450$Stage))[4]) & !is.na(pat_info_meth450$Stage)),
                                  sum(is.na(pat_info_meth450$Stage))),
                      stringsAsFactors = F)
df_plot$Stage <- factor(df_plot$Stage, levels = rev(levels(df_plot$Stage)))
p <- ggplot(data=df_plot, aes(x=Stage,y=Frequency,fill=Stage)) +
  geom_bar(stat="identity") + 
  scale_fill_brewer(palette="Dark2") + 
  scale_y_continuous(expand = c(0, 0), limits=c(0,200)) +
  coord_flip() + 
  geom_text(aes(label=Frequency),nudge_y=7,size=3.5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.y = element_blank(),
        axis.line.x = element_line(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm"),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Stage distribution") 
ggsave(paste(cancertype,"_Stage_distr.jpg",sep = ""),p, height = 4, width = 9)

p <- qplot(pat_info_meth450$Age, geom="histogram") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = "none",
        axis.line.x = element_line(), axis.line.y = element_line(),
        plot.margin = unit(c(1,1,1,1),"cm"),
        plot.title = element_text(hjust = 0.5)) + 
  labs(x ="Age", y="Frequency") + 
  ggtitle("Age distribution") 
ggsave(paste(cancertype,"_Age_distr.jpg",sep = ""),p, height = 6, width = 6)



# Model Survival ifo Age, Gender and Stage
dat_cox_0 <- data.frame(Surv = pat_info_meth450$Survival,
                      Age = pat_info_meth450$Age,
                      Gender = pat_info_meth450$Gender,
                      Stage = pat_info_meth450$Stage)

fit.coxph_0 <- coxph(Surv ~ Age + Gender + Stage, 
                     data = dat_cox_0)

print(summary(fit.coxph_0))


# Model Survival ifo methylation per probe (with Age, Gender and Stage as cofounders)
# for (i in 1:nrow(dat_meth450_tumor_M)){
myfun <- function(i){
  # print(i)
  if (i%%20000==0){print(i)}
  dat_cox <- data.frame(Surv = pat_info_meth450$Survival,
                        Mval = dat_meth450_tumor_M_filt[i,],
                        Age = pat_info_meth450$Age,
                        Gender = pat_info_meth450$Gender,
                        Stage = pat_info_meth450$Stage)
  
  fit_coxph <- coxph(Surv ~ Mval + Age + Gender + Stage, 
                     data = dat_cox)
  
  results <- c(rownames(dat_meth450_tumor_M_filt)[i],
               mean(dat_meth450_tumor_filt[i,],na.rm=T),sd(dat_meth450_tumor_filt[i,],na.rm=T),
               summary(fit_coxph)$coefficients[1,1],
               summary(fit_coxph)$conf.int[1,c(1,3,4)],
               summary(fit_coxph)$coefficients[1,5])
  return(results)
}
print(nrow(dat_meth450_tumor_M_filt))
results_cox_list <- mclapply(c(1:nrow(dat_meth450_tumor_M_filt)), myfun, mc.cores=20, mc.preschedule = T)
results_cox <- do.call(rbind,results_cox_list)
results_cox <- data.frame(results_cox,stringsAsFactors = F)
colnames(results_cox) <- c("ID","Mean_meth_tumor","SD_meth_tumor","Coef","HR","HR_lower95","HR_upper95","Pval")

results_cox$Mean_meth_tumor <- as.numeric(results_cox$Mean_meth_tumor)
results_cox$SD_meth_tumor <- as.numeric(results_cox$SD_meth_tumor)

results_cox$Coef <- as.numeric(results_cox$Coef)
results_cox$Abs_Coef <- abs(results_cox$Coef)
results_cox$HR <- as.numeric(results_cox$HR)
results_cox$HR_lower95 <- as.numeric(results_cox$HR_lower95)
results_cox$HR_upper95 <- as.numeric(results_cox$HR_upper95)
results_cox$Pval <- as.numeric(results_cox$Pval)

results_cox <- results_cox[sort(results_cox$Pval,index.return=T)$ix,]
results_cox$FDR <- p.adjust(results_cox$Pval)
head(results_cox)


meth450_annot_grch38_sorted <- meth450_annot_grch38
results_cox_sorted <- results_cox[sort(results_cox$ID,index.return=T)$ix,]
results_cox_sorted <- results_cox_sorted[results_cox_sorted$ID%in%meth450_annot_grch38_sorted$probeID,]
meth450_annot_grch38_sorted <- meth450_annot_grch38_sorted[meth450_annot_grch38_sorted$probeID%in%results_cox_sorted$ID,]
# Deficit in results_cox_sorted due to SNPs, i.e. rs id's

if(sum(meth450_annot_grch38_sorted$probeID!=results_cox_sorted$ID)!=0){
  stop("Problem with meta annotation")
}

dat_meth450_control_filt <- dat_meth450_control[rownames(dat_meth450_control)%in%results_cox_sorted$ID,]
dat_meth450_control_filt_sorted <- dat_meth450_control_filt[sort(rownames(dat_meth450_control_filt),index.return=T)$ix,]

results_cox_sorted$Mean_meth_control <- rowMeans(dat_meth450_control_filt_sorted,na.rm=T)
results_cox_sorted$SD_meth_control <- apply(dat_meth450_control_filt_sorted,1,sd,na.rm=T)

results_cox_sorted$Mean_meth_diff <- results_cox_sorted$Mean_meth_tumor - results_cox_sorted$Mean_meth_control

results_cox_sorted$Chr <- meth450_annot_grch38_sorted$chrm
results_cox_sorted$Pos <- meth450_annot_grch38_sorted$start
results_cox_sorted$Gene <- meth450_annot_grch38_sorted$GeneNames
results_cox_sorted$pos_TSS <- meth450_annot_grch38_sorted$PosTSS
results_cox_sorted$CGI <- meth450_annot_grch38_sorted$CGI
results_cox_sorted$pos_CGI <- meth450_annot_grch38_sorted$CGIPosition

results_cox_sorted <- results_cox_sorted[,c("ID","Chr","Pos","Gene",
                                            "Mean_meth_tumor","SD_meth_tumor","Mean_meth_control","SD_meth_control","Mean_meth_diff",
                                            "Abs_Coef","Coef","HR","Pval","FDR","HR_lower95","HR_upper95",
                                            "pos_TSS","CGI","pos_CGI")]

### Additional annotation
grch38 <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
annotation <- getBM(mart=grch38,attribute=c("ensembl_gene_id","external_gene_name",
                                            "transcription_start_site"))
dim(annotation)
annotation <- annotation[!duplicated(annotation$external_gene_name),]
dim(annotation)


# Expression Data
dat_expression <- read.table(paste("TCGA-",cancertype,".htseq_fpkm-uq.tsv",sep=""),
                             header=T,stringsAsFactors= F)


dat_expression_rn <- gsub("\\.[0-9]+$","",dat_expression[,1])
rownames(dat_expression) <- dat_expression_rn
dat_expression <- dat_expression[,-1]
colnames(dat_expression) <- gsub("A$","",colnames(dat_expression))
dim(dat_expression)
dat_expression <- dat_expression[,!grepl("B$",colnames(dat_expression))]
dim(dat_expression)
colnames(dat_expression) <- gsub("\\.","\\_",colnames(dat_expression))

dat_expression <- dat_expression[,colnames(dat_expression)%in%colnames(dat_meth450_tumor)]
dat_expression <- dat_expression[,sort(colnames(dat_expression),index.return=T)$ix]

dat_expression <- sapply(dat_expression,as.numeric)
rownames(dat_expression) <- dat_expression_rn

dim(dat_expression)
head(dat_expression[,1:10])

subset_has_expression <- colnames(dat_meth450_tumor)%in%colnames(dat_expression)

dat_meth450_tumor_filt_sorted <- dat_meth450_tumor_filt[sort(rownames(dat_meth450_tumor_filt),index.return=T)$ix,]
dat_meth450_tumor_filt_sorted <- dat_meth450_tumor_filt_sorted[rownames(dat_meth450_tumor_filt_sorted)%in%
                                                                 meth450_annot_grch38_sorted$probeID,]

myfun2 <- function(i){
  if(i%%10000==0){print(i)}
  if (results_cox_sorted$Gene[i]!=" "){
    genes_oi <- unique(unlist(strsplit(results_cox_sorted$Gene[i],";")))
    n_annot <- sum(annotation$external_gene_name%in%genes_oi)
    if (n_annot!=0){
      if (n_annot>1){
        
        correl <- NULL
        for (j in 1:length(genes_oi)){
          if (identical(annotation$ensembl_gene_id[annotation$external_gene_name==genes_oi[j]],character(0))){
            cor_exprs <- c(NA,NA)
          } else {
            # correlation gene expression
            if (annotation$ensembl_gene_id[annotation$external_gene_name==genes_oi[j]]%in%rownames(dat_expression)){
              correl <- rbind(correl,
                              c(cor.test(dat_meth450_tumor_filt_sorted[i,subset_has_expression],
                                         dat_expression[annotation$ensembl_gene_id[annotation$external_gene_name==genes_oi[j]][1],])$p.value,
                                cor.test(dat_meth450_tumor_filt_sorted[i,subset_has_expression],
                                         dat_expression[annotation$ensembl_gene_id[annotation$external_gene_name==genes_oi[j]][1],])$estimate))
            } else {
              correl <- rbind(correl,c(NA,NA))
            }
          }
        }
        cor_exprs <- c(paste(correl[,1],collapse=";"),
                       paste(correl[,2],collapse=";"))
      } else {
        
        if (identical(annotation$ensembl_gene_id[annotation$external_gene_name==genes_oi],character(0))){
          cor_exprs <- c(NA,NA)
        } else {
         # correlation gene expression
          if (annotation$ensembl_gene_id[annotation$external_gene_name==genes_oi]%in%rownames(dat_expression)){
            cor_exprs <- c(cor.test(dat_meth450_tumor_filt_sorted[i,subset_has_expression],
                                    dat_expression[annotation$ensembl_gene_id[annotation$external_gene_name==genes_oi][1],])$p.value,
                           cor.test(dat_meth450_tumor_filt_sorted[i,subset_has_expression],
                                    dat_expression[annotation$ensembl_gene_id[annotation$external_gene_name==genes_oi][1],])$estimate)
          } else {
            cor_exprs <- c(NA,NA)
          }
        }

        
      }
      
    } else {
      cor_exprs <- c(NA,NA)
    }
  } else {
    cor_exprs <- c(NA,NA)
  }
  
  return(cor_exprs)
}
print(nrow(results_cox_sorted))
# results <- mclapply(1:1000,myfun2,mc.cores=30,mc.preschedule = F)
results_list_cor <- mclapply(1:nrow(results_cox_sorted),myfun2,mc.cores=20,mc.preschedule = T)
results_cor <- do.call(rbind,results_list_cor)
results_cor <- data.frame(results_cor)
results_cor <- results_cor[,c(2,1)]
colnames(results_cor) <- c("cor","pval_cor")

results_cor$cor <- as.character(results_cor$cor)
results_cor$pval_cor <- as.character(results_cor$pval_cor)

dim(results_cox_sorted)
results_cox_sorted <- cbind(results_cox_sorted,results_cor)
dim(results_cox_sorted)

results_cox <- results_cox_sorted[sort(results_cox_sorted$Abs_Coef,
                                       index.return=T,decreasing = T)$ix,]
head(results_cox)


# iAI data 
dat_iAI <- read.table("/data/louisc/project_Tine/Data_new/Data_plot_allchr_iAI_sign_SNPs_allstages.txt",
                      header=T,stringsAsFactors =F,quote = "",sep="\t")

genes_iAI <- do.call(rbind,(strsplit(dat_iAI$gene,",")))
genes_iAI <- gsub("\\*","",genes_iAI)
myfun3 <- function(i){
  if(i%%10000==0){print(i)}
  if (results_cox$Gene[i]!=" "){
    genes_oi <- unique(unlist(strsplit(results_cox$Gene[i],";")))
    n_annot <- length(unique(gsub("\\*","",dat_iAI$gene[rowSums(matrix(genes_iAI%in%genes_oi,ncol=ncol(genes_iAI)))>0])))
    if (n_annot==0){
      mean_med_ai <- NA
      pval_med_ai <- NA
    } else if (n_annot>1){
      # print("more than one gene")
      # print(genes_iAI[rowSums(matrix(genes_iAI%in%genes_oi,ncol=ncol(genes_iAI)))>0,])
      # print(dat_iAI$Mean_AI[rowSums(matrix(genes_iAI%in%genes_oi,ncol=ncol(genes_iAI)))>0])
      # print(genes_oi)
      genes_temp <- unique(c(genes_iAI[rowSums(matrix(genes_iAI%in%genes_oi,ncol=ncol(genes_iAI)))>0,]))
      order_temp <- sort(genes_temp,index.return=T)$ix
      mean_med_ai <- paste(dat_iAI$Mean_AI[rowSums(matrix(genes_iAI%in%genes_oi,ncol=ncol(genes_iAI)))>0][order_temp][sort(genes_temp)%in%genes_oi],collapse=";")
      pval_med_ai <- paste(dat_iAI$AI_harm_pval[rowSums(matrix(genes_iAI%in%genes_oi,ncol=ncol(genes_iAI)))>0][order_temp][sort(genes_temp)%in%genes_oi],collapse=";")
    } else {
      mean_med_ai <- unique(dat_iAI$Mean_AI[rowSums(matrix(genes_iAI%in%genes_oi,ncol=ncol(genes_iAI)))>0])
      pval_med_ai <- unique(dat_iAI$AI_harm_pval[rowSums(matrix(genes_iAI%in%genes_oi,ncol=ncol(genes_iAI)))>0])
    }
  }
  obj_return <- c(mean_med_ai,pval_med_ai)
  return(obj_return)
}
print(nrow(results_cox))
# results <- mclapply(1:1000,myfun2,mc.cores=30,mc.preschedule = F)
results_list_iAI <- mclapply(1:nrow(results_cox),myfun3,mc.cores=30,mc.preschedule = T)
results_iAI <- do.call(rbind,results_list_iAI)
results_iAI <- data.frame(results_iAI)
colnames(results_iAI) <- c("mean_iAI","harm_pval_iAI")

dim(results_cox)
results_cox_final <- cbind(results_cox,results_iAI)
dim(results_cox_final)

sum(results_cox_final$FDR<0.05)

write.table(results_cox_final,file=paste(cancertype,"_results_cox.txt",sep=""),col.names = T,row.names = F,sep = "\t",quote = F)

################
# manual check #
################

CpG_id <- "cg13108194"

# genomic annotation
meth450_annot_grch38[meth450_annot_grch38$probeID==CpG_id,
                     colnames(meth450_annot_grch38)%in%c("probeID","chrm","start",
                                                         "end","adressA","GeneNames",
                                                         "PosTSS","CGI","CGIPosition")]

# meth statistics
mean(dat_meth450_control[rownames(dat_meth450_control)==CpG_id,],na.rm=T)
sd(dat_meth450_control[rownames(dat_meth450_control)==CpG_id,],na.rm=T)
mean(dat_meth450_tumor[rownames(dat_meth450_control)==CpG_id,],na.rm=T)
sd(dat_meth450_tumor[rownames(dat_meth450_control)==CpG_id,],na.rm=T)

# correlation
gene_id_test <- c("CARS2","ING1")
for (i in 1:length(gene_id_test)){
  print(annotation[annotation$external_gene_name==gene_id_test[i],])
}

mask_dat_meth450_tumor_hasexpression <- colnames(dat_meth450_tumor)%in%colnames(dat_expression)
sum(colnames(dat_expression)==colnames(dat_meth450_tumor)[mask_dat_meth450_tumor_hasexpression])

# dat_expression[rownames(dat_expression)=="ENSG00000149115",1:6]
# dat_meth450_tumor[,mask_dat_meth450_tumor_hasexpression][rownames(dat_meth450_control)==CpG_id,1:6]

ens_id_test <- c("ENSG00000134905","ENSG00000153487")
for (i in 1:length(ens_id_test)){
  if(ens_id_test[i]%in%rownames(dat_expression)){
    print(cor.test(dat_expression[rownames(dat_expression)==ens_id_test[i],],
                   dat_meth450_tumor[rownames(dat_meth450_control)==CpG_id,mask_dat_meth450_tumor_hasexpression]))
  } else {
    print(paste("No expression for:",ens_id_test[i]))
  }

}


# dat iAI 
for (i in 1:length(gene_id_test)){
  print(dat_iAI[grep(gene_id_test[i],dat_iAI$gene),1:10])
}

         
