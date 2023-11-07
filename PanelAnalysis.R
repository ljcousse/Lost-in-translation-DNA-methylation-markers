#####################
#####################
##                 ##
##  PanelAnalysis  ##
##                 ##
#####################
#####################

# Get cancer type
cancertype <- commandArgs(TRUE)
print(cancertype)

if(identical(args,character(0))){
  cancertype <- "KIRC"
}

# Library
# library("survminer")
library(xlsx)
library(parallel)
library(survival)
library(biomaRt)

results_cox <- read.table(file=paste(cancertype,"_results_cox.txt",sep=""),
                          header=T,sep = "\t",quote = "")
panel_probes <- read.xlsx(file="Suppl_Table5.xlsx",1,header=T)

print(head(results_cox))
print(panel_probes)

dim(panel_probes)
sum(panel_probes$Probe.ID%in%results_cox$ID)
which(results_cox$ID%in%panel_probes$Probe.ID)

results_cox_panel <- results_cox[results_cox$ID%in%panel_probes$Probe.ID,]
write.table(results_cox_panel,file=paste(cancertype,"_results_cox_panel.txt",sep=""),
            col.names=T,row.names = F,sep="\t",quote=F)

#################
#################
##             ##
## Meth & CNV  ##
##             ##
#################
#################


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



#################
#################
##             ##
## Meth & CNV  ##
##             ##
#################
#################

dat_cox <- data.frame(Surv = pat_info_meth450$Survival,
                      Age = pat_info_meth450$Age,
                      Gender = pat_info_meth450$Gender,
                      Stage = pat_info_meth450$Stage)

Mval_mat <- NULL
for (i in 1:nrow(panel_probes)){
  Mval_mat <- cbind(Mval_mat,dat_meth450_tumor_M_filt[rownames(dat_meth450_tumor_M_filt)==panel_probes$Probe.ID[i],])
}
colnames(Mval_mat) <- paste(rep("Mval_",nrow(panel_probes)),panel_probes$Probe.ID,sep="")

dat_cox_panel <- cbind(dat_cox,Mval_mat)

# Nog niet voor parallelisatie lopen anders werkt het niet meer
# formula <- as.formula(paste("Surv ~ ",paste(colnames(dat_cox)[grep("Mval",colnames(dat_cox))],collapse=" + ")," + Age + Gender + Stage"))
# fit.coxph <- coxph(formula, 
#                    data = dat_cox)
# summary(fit.coxph)

wh_in_panel <- which(rownames(dat_meth450_tumor_M_filt)%in%panel_probes$Probe.ID)
myfun <- function(i){
  # print(i)
  if (i%%20000==0){print(i)}
  if (i%in%wh_in_panel){
    return(NULL)
  }
  
  dat_cox_i <- cbind(dat_cox_panel,dat_meth450_tumor_M_filt[i,])
  colnames(dat_cox_i) <- c(colnames(dat_cox_panel),"Mval_new")
  
  # print(head(dat_cox_i))

  formula_i <- as.formula(paste("Surv ~ ",paste(rev(colnames(dat_cox_i)[grep("Mval",colnames(dat_cox_i))]),collapse=" + "),
                                " + Age + Gender + Stage",sep=""))
  # print(formula_i)
  
  fit_coxph_i <- coxph(formula_i,
                       data = dat_cox_i)
  # print(fit.coxph_i)
  # print(summary(fit_coxph_i))
  results <- c(rownames(dat_meth450_tumor_M_filt)[i],
               mean(dat_meth450_tumor_filt[i,],na.rm=T),sd(dat_meth450_tumor_filt[i,],na.rm=T),
               summary(fit_coxph_i)$coefficients[rownames(summary(fit_coxph_i)$coefficients)=="Mval_new",1],
               summary(fit_coxph_i)$conf.int[rownames(summary(fit_coxph_i)$conf.int)=="Mval_new",c(1,3,4)],
               summary(fit_coxph_i)$coefficients[rownames(summary(fit_coxph_i)$coefficients)=="Mval_new",5])
  
  # fit.coxph_i <- summary(coxph(formula_i,data = dat_cox_i))
  # 
  # results <- c(rownames(dat_meth450_M_filt)[i],
  #              mean(dat_meth450_filt[i,],na.rm=T),sd(dat_meth450_filt[i,],na.rm=T),
  #              fit.coxph_i$conf.int[rownames(fit.coxph_i$conf.int)=="Mval_new",c(1,3,4)],
  #              fit.coxph_i$coefficients[rownames(fit.coxph_i$coefficients)=="Mval_new",5])
  # print(results)
  return(results)
}
print(nrow(dat_meth450_tumor_M_filt))
# results_cox_panel_list <- mclapply(c(1:20), myfun, mc.cores=20, mc.preschedule = T)
results_cox_panel_list <- mclapply(c(1:nrow(dat_meth450_tumor_M_filt)), myfun, mc.cores=20, mc.preschedule = T)
# results_cox_panel <- results
results_cox_panel <- do.call(rbind,results_cox_panel_list)
results_cox_panel <- data.frame(results_cox_panel,stringsAsFactors = F)
colnames(results_cox_panel) <- c("ID","Mean_meth_tumor","SD_meth_tumor","Coef","HR",
                                 "HR_lower95","HR_upper95","Pval")
results_cox_panel$Coef <- as.numeric(results_cox_panel$Coef)
results_cox_panel$Abs_Coef <- abs(results_cox_panel$Coef)
results_cox_panel$HR <- as.numeric(results_cox_panel$HR)
results_cox_panel$HR_lower95 <- as.numeric(results_cox_panel$HR_lower95)
results_cox_panel$HR_upper95 <- as.numeric(results_cox_panel$HR_upper95)
results_cox_panel$Pval <- as.numeric(results_cox_panel$Pval)

results_cox_panel$Mean_meth_tumor <- as.numeric(results_cox_panel$Mean_meth_tumor)
results_cox_panel$SD_meth_tumor <- as.numeric(results_cox_panel$SD_meth_tumor)

results_cox_panel <- results_cox_panel[sort(results_cox_panel$Pval,index.return=T)$ix,]
results_cox_panel$FDR <- p.adjust(results_cox_panel$Pval)
head(results_cox_panel)

results_cox_panel_sorted <- results_cox_panel[sort(results_cox_panel$ID,index.return=T)$ix,]
results_cox_panel_sorted <- results_cox_panel_sorted[results_cox_panel_sorted$ID%in%meth450_annot_grch38$probeID,]

results_cox_sorted <- results_cox[sort(as.character(results_cox$ID),index.return=T)$ix,]

mask_not_in_panel <- results_cox_sorted$ID%in%results_cox_panel_sorted$ID
sum(results_cox_panel_sorted$ID==results_cox_sorted$ID[mask_not_in_panel])

results_cox_panel_sorted$Chr <- as.character(results_cox_sorted$Chr[mask_not_in_panel])
results_cox_panel_sorted$Pos <- results_cox_sorted$Pos[mask_not_in_panel]
results_cox_panel_sorted$Gene <- as.character(results_cox_sorted$Gene[mask_not_in_panel])

results_cox_panel_sorted$Mean_meth_control <- results_cox_sorted$Mean_meth_control[mask_not_in_panel]
results_cox_panel_sorted$SD_meth_control <- results_cox_sorted$SD_meth_control[mask_not_in_panel]
results_cox_panel_sorted$Mean_meth_diff <- results_cox_sorted$Mean_meth_diff[mask_not_in_panel]

results_cox_panel_sorted$pos_TSS <- results_cox_sorted$pos_TSS[mask_not_in_panel]
results_cox_panel_sorted$CGI <- results_cox_sorted$CGI[mask_not_in_panel]
results_cox_panel_sorted$pos_CGI <- results_cox_sorted$pos_CGI[mask_not_in_panel]
results_cox_panel_sorted$cor <- as.numeric(as.character(results_cox_sorted$cor[mask_not_in_panel]))
results_cox_panel_sorted$pval_cor <- as.numeric(as.character(results_cox_sorted$pval_cor[mask_not_in_panel]))
results_cox_panel_sorted$mean_iAI <- as.numeric(as.character(results_cox_sorted$mean_iAI[mask_not_in_panel]))
results_cox_panel_sorted$harm_pval_iAI <- as.numeric(as.character(results_cox_sorted$harm_pval_iAI[mask_not_in_panel]))

results_cox_panel_sorted <- results_cox_panel_sorted[,colnames(results_cox_sorted)]
results_cox_panel_final <- results_cox_panel_sorted[sort(results_cox_panel_sorted$Abs_Coef,
                                                         decreasing = T,index.return=T)$ix,]
write.table(results_cox_panel_final,file=paste(cancertype,"_results_cox_panel.txt",sep=""),
            col.names = T,row.names = F,sep = "\t",quote = F)


results_cox_panel_sign <- results_cox_panel_sorted[results_cox_panel_sorted$FDR<0.05,]
results_cox_panel_sign_sort <- results_cox_panel_sign[sort(results_cox_panel_sign$Abs_Coef,
                                                           decreasing = T,index.return=T)$ix,]
print(dim(results_cox_panel_sign))

results_cox_sign <- results_cox[results_cox$FDR<0.05,]
results_cox_sign_sort <- results_cox_sign[sort(results_cox_sign$Abs_Coef,
                                               decreasing = T,index.return=T)$ix,]
print(dim(results_cox_sign))

sum(results_cox_panel_sign_sort$ID%in%results_cox_sign_sort$ID)
CG_oi <- results_cox_panel_sign_sort$ID[results_cox_panel_sign_sort$ID%in%results_cox_sign_sort$ID]
CG_oi
results_cox_sign_sort[results_cox_sign_sort$ID%in%CG_oi,]

# 
# # Plotting
# "cg15211499"
# 
# for (i in 1:nrow(results_cox_panel_sign_sort)){
#   wh <- which(rownames(dat_meth450_M_filt)==results_cox_panel_sign_sort$ID[i])
# 
#   dat_cox_wh <- cbind(dat_cox,dat_meth450_M_filt[wh,])
#   dat_cox_wh <- cbind(dat_cox_wh,as.factor(dat_meth450_M_filt[wh,]>mean(dat_meth450_M_filt[wh,])))
#   colnames(dat_cox_wh) <- c(colnames(dat_cox),"Mval_new","Mval_new_fact")
# 
#   formula_i <- as.formula(paste("Surv ~ Mval_new_fact + ",paste(colnames(dat_cox_wh)[grep("Mval",colnames(dat_cox_i))],collapse=" + "),
#                                 " + Age + Gender + Stage",sep=""))
#   fit_plot <- survfit(formula_i, data = dat_cox_wh)
#   
#   df <- data.frame(val = dat_meth450_M_filt[wh,],stringsAsFactors = F)
#   tiff(paste("/data/louisc/Project_KimSmits/Plots/",results_cox_panel_sign_sort$ID[i],"_survival.tiff",sep=""),
#        units="in", height = 4, width=8, res=75)
#   plot1 <- ggsurvplot(fit_plot, data = dat_cox_wh)
#   plot2 <- ggplot(df,aes(x=val)) + geom_histogram() + geom_vline(xintercept=mean(dat_meth450_M_filt[wh,]),color="red") 
#   grid.arrange(plot1$plot, plot2, ncol=2)
#   dev.off()
# }
# 
# 
# cg_to_plot <- "cg15211499"
# for (i in 1:length(cg_to_plot)){
#   wh <- which(rownames(dat_meth450_M_filt)==cg_to_plot[i])
#   
#   dat_cox_wh <- cbind(dat_cox,dat_meth450_M_filt[wh,])
#   dat_cox_wh <- cbind(dat_cox_wh,as.factor(dat_meth450_M_filt[wh,]>mean(dat_meth450_M_filt[wh,])))
#   colnames(dat_cox_wh) <- c(colnames(dat_cox),"Mval_new","Mval_new_fact")
#   
#   formula_i <- as.formula(paste("Surv ~ Mval_new_fact + ",paste(colnames(dat_cox_wh)[grep("Mval",colnames(dat_cox_i))],collapse=" + "),
#                                 " + Age + Gender + Stage",sep=""))
#   fit_plot <- survfit(formula_i, data = dat_cox_wh)
#   
#   df <- data.frame(val = dat_meth450_M_filt[wh,],stringsAsFactors = F)
#   tiff(paste("/data/louisc/Project_KimSmits/Plots/",results_cox_panel_sign_sort$ID[i],"_survival.tiff",sep=""),
#        units="in", height = 4, width=8, res=75)
#   plot1 <- ggsurvplot(fit_plot, data = dat_cox_wh)
#   plot2 <- ggplot(df,aes(x=val)) + geom_histogram() + geom_vline(xintercept=mean(dat_meth450_M_filt[wh,]),color="red") 
#   grid.arrange(plot1$plot, plot2, ncol=2)
#   dev.off()
# }
