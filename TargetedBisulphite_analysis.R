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
library(survival)
library(scales)
library(parallel)
library(survival)
library(biomaRt)
library(survminer)
library(gridExtra)
library(ggplot2)

######################
# plotGenes function #
######################

plotPlots_AIC <- function(dat_gene,locs_gene,gene,MSP_amplicon,primer_annot_gene,gene_cor=NULL,
                          gene_obj,transcript_ids,chrom,object_box,elevate_amplicon=matrix(c(10^6,10^6),nrow=1),max_y=NULL){
  chrom_gene <- as.character(data.frame(locs_gene)[1,1])

  col_gene <- c("cadetblue","blue","navy","lightgrey")

  dat_gene <- dat_gene[sort(dat_gene$ID,index.return=T)$ix,]
  
  multiplier <- 2
  cex_axis <- 1*multiplier
  cex_lab <- 1*multiplier
  cex_main <- 1*multiplier
  cex_legend <- 1*multiplier
  cex_labels <- 0.8*multiplier
  title_line <- 4.5
  diff_units_transcript_id <- 5
  
  cex_gene <- 1*multiplier
  
  lwd_gene <- 1*multiplier
  
  if (is.null(max_y)){
    max_y <- max(c(dat_gene$delta_AIC,6))
  }

  # Plot delta AIC
  par(mar=c(6,8,2,40))
  plot(dat_gene$Pos,dat_gene$delta_AIC,type="n",
       ylab="",xlab="",
       xaxt="n",yaxt="n",cex.main=cex_main,
       cex=cex_gene,col=col_gene,
       ylim=c(min(c(dat_gene$delta_AIC,-6)),max_y),
       xlim=c(min(primer_annot_gene$begin),max(primer_annot_gene$end)))
  

  cols <- NULL
  s <- 0
  for (i in 1:nrow(primer_annot_gene)){
    if (sum(dat_gene$Amplicon==primer_annot_gene$name[i])>0){
      if (i == 1){
        cols <- c(cols,1)
        s <- 1
      } else {
        if (i%in%elevate_amplicon[,1]){
          cols <- c(cols,1+elevate_amplicon[elevate_amplicon[,1]==i,2])
          s <- i
        } else if ((primer_annot_gene$end[s]>primer_annot_gene$begin[i])){
          cols <- c(cols,1+i-s)
        } else {
          cols <- c(cols,1)
          s <- i
        }
      }
    } else {
        cols <- c(cols,4)
    }
  }
  
  level_i <- NULL
  s <- 0
  for (i in 1:nrow(primer_annot_gene)){
    if (i == 1){
      level_i <- c(level_i,1)
      s <- 1
    } else {
      if (i%in%elevate_amplicon[,1]){
        level_i <- c(level_i,1+elevate_amplicon[elevate_amplicon[,1]==i,2])
        s <- i
      } else if ((primer_annot_gene$end[s]>primer_annot_gene$begin[i])){
        level_i <- c(level_i,1+i-s)
      } else {
        level_i <- c(level_i,1)
        s <- i
      }
    }
  }
  
  ydelta_text <- 1.5*(abs(par("usr")[4])+abs(par("usr")[3]))/20
  level_0 <- max(c(0,dat_gene$delta_AIC)) + ydelta_text
  level_delta <- 1.5*(abs(par("usr")[4])+abs(par("usr")[3]))/20*2
  
  for (i in 1:nrow(primer_annot_gene)){
    level <- level_0 + (level_i[i]-1)*level_delta
    segments(primer_annot_gene$begin[i],level,primer_annot_gene$end[i],level,col=col_gene[cols[i]],lwd=2)
    text(primer_annot_gene$begin[i]+(primer_annot_gene$end[i]-primer_annot_gene$begin[i])/2,level+ydelta_text,
           primer_annot_gene$name[i],adj=0.5,offset=0,col=col_gene[cols[i]],cex=cex_lab)
  }
  
  for (i in 1:nrow(primer_annot_gene)){
    if (sum(dat_gene$Amplicon==primer_annot_gene$name[i])>0){
      dat_gene_plot_line_i <- dat_gene[dat_gene$Amplicon==primer_annot_gene$name[i],]
      lines(dat_gene_plot_line_i$Pos,dat_gene_plot_line_i$delta_AIC,lwd=2,col=col_gene[cols[i]])
      points(dat_gene_plot_line_i$Pos,dat_gene_plot_line_i$delta_AIC,pch=18,col=col_gene[cols[i]])
    }
  }
  
  abline(h=0,lty=1,col="grey",lwd=lwd_gene)
  abline(h=2*log(0.05),lty=2,col="grey",lwd=lwd_gene)
  mtext(paste("Chromosome ",chrom_gene,sep=""),1,cex=cex_lab,line=title_line)
  mtext(expression(paste(Delta,"AIC")),2,cex=cex_lab,line=title_line+1)
  mtext("RL",4,cex=cex_lab,line=title_line+3.5)
  axis(1,at=pretty(c(par("usr")[1],par("usr")[2])),labels=pretty(c(par("usr")[1],par("usr")[2])),cex.axis=cex_axis,padj=0.75)
  y_ax_t <- pretty(c(par("usr")[3],par("usr")[4]))
  y_ax_t <- y_ax_t[y_ax_t<1]
  axis(2,at=y_ax_t,labels=format(y_ax_t,nsmall=1),cex.axis=cex_axis,hadj=1,las=1)
  axis(4,at=c(0,-6,-12,-18),labels=sprintf("%0.4f",round(exp(c(0,-6,-12,-18)/2),digits=4)),cex.axis=cex_axis,las=1)
  rect(MSP_amplicon[1],par("usr")[3],MSP_amplicon[2],par("usr")[4],col=alpha("cadetblue",alpha=0.20),border=NA,lwd=lwd_gene)
  

  
  # Check for gene structures not in plotting range
  to_remove <- NULL
  for (i in 1:length(transcript_ids)){
    min_x <- min(c(gene_obj$exon_chrom_start[gene_obj$ensembl_transcript_id==transcript_ids[i]],
                   gene_obj$exon_chrom_end[gene_obj$ensembl_transcript_id==transcript_ids[i]]))
    max_x <- max(c(gene_obj$exon_chrom_start[gene_obj$ensembl_transcript_id==transcript_ids[i]],
                   gene_obj$exon_chrom_end[gene_obj$ensembl_transcript_id==transcript_ids[i]]))
    if ((min_x>par("usr")[2]) | (max_x<par("usr")[1])){
      to_remove <- c(to_remove,i)
    }
  }
  
  if (!is.null(to_remove)){
    transcript_ids <- transcript_ids[-to_remove] 
  }
  
  # Plot Gene structure
  par(mar=c(6,8,2,40))
  plot(dat_gene$Pos,dat_gene$delta_AIC,type="n",ylim=c(0,length(transcript_ids)+1),xlab="",
       ylab="",bty="n",yaxt="n",xaxt="n")
  mtext(paste("Chromosome ",chrom_gene,sep=""),1,cex=cex_lab,line=title_line)
  axis(1,at=pretty(c(par("usr")[1],par("usr")[2])),labels=pretty(c(par("usr")[1],par("usr")[2])),cex.axis=cex_axis,padj=0.75)
  
  
  for (i in 1:length(transcript_ids)){
    dat_transcript <- gene_obj[gene_obj$ensembl_transcript_id==transcript_ids[i],]
    dat_transcript <- dat_transcript[sort(as.numeric(dat_transcript$exon_chrom_start),index.return=T)$ix,]
    # print(dat_transcript)
    for (j in 1:nrow(dat_transcript)){
      rect(dat_transcript$exon_chrom_start[j],i+0.35,dat_transcript$exon_chrom_end[j],i-0.35,col="red")
      if (j>1){
        segments(dat_transcript$exon_chrom_end[j-1],i,dat_transcript$exon_chrom_start[j],i)
      }
    }
    if(sum((dat_transcript$exon_chrom_start-dat_transcript$transcription_start_site[1])>0)>((nrow(dat_transcript)/2)-0.05)){ #-0.05 for in the case you only have 2 exons
      # print(par("usr"))
      x_poly <- c(dat_transcript$transcription_start_site[1],dat_transcript$transcription_start_site[1],
                  dat_transcript$transcription_start_site[1]+((par("usr")[2]-par("usr")[1])/100))
      y_poly <- c(i-0.4,i+0.4,i)
      # print(x_poly)
      # print(y_poly)
      polygon(x_poly,y_poly,col="green")
      text(dat_transcript$transcription_start_site[1]-((par("usr")[2]-par("usr")[1])/100*diff_units_transcript_id),i,
           labels=dat_transcript$ensembl_transcript_id[1],cex=cex_labels)
    } else {
      # print(par("usr"))
      x_poly <- c(dat_transcript$transcription_start_site[1],dat_transcript$transcription_start_site[1],
                  dat_transcript$transcription_start_site[1]-((par("usr")[2]-par("usr")[1])/100))
      y_poly <- c(i-0.4,i+0.4,i)
      # print(x_poly)
      # print(y_poly)
      polygon(x_poly,y_poly,col="green")
      text(dat_transcript$transcription_start_site[1]+((par("usr")[2]-par("usr")[1])/100*diff_units_transcript_id),i,
           labels=dat_transcript$ensembl_transcript_id[1],cex=cex_labels)
    }
  }
  
  par(xpd=T)
  if (length(to_remove)==1){
    legend(par("usr")[2]+((par("usr")[2]-par("usr")[1])/100*2.5),par("usr")[4],
           c("TSS","Exon","Intron",paste(length(to_remove)," transcript was ommited (not ",sep=""),"in plotting range)"),
           pch=c(15,15,NA,NA,NA),lty=c(NA,NA,1,NA,NA),col=c("green","red","black",NA,NA),cex=cex_legend)
  } else {
    legend(par("usr")[2]+((par("usr")[2]-par("usr")[1])/100*2.5),par("usr")[4],
           c("TSS","Exon","Intron",paste(length(to_remove)," transcripts were ommited (not ",sep=""),"in plotting range)"),
           pch=c(15,15,NA,NA,NA),lty=c(NA,NA,1,NA,NA),col=c("green","red","black",NA,NA),cex=cex_legend)
  }
  par(xpd=F)
  
}

tiff(paste(output_dir,"LegendFirstPanel.jpg",sep=""),width=2048+512+256,height=1024*2/3)
par(mar=c(6,8,2,40))
plot(1,1,type="n")
multiplier <- 2
cex_legend <- 1*multiplier
col_gene <- c("cadetblue","blue","navy","lightgrey")
par(xpd=T)
legend(par("usr")[2]+((par("usr")[2]-par("usr")[1])/100*2.5),par("usr")[4],
       c(expression(paste(Delta,"AIC between base model and")),
         "   methylation model","","","Amplicons","","","Amplicons without sufficient data","",
         expression(paste(Delta,"AIC=0; Base model equally")),
         "probable as methylation model",
         expression(paste(Delta,"AIC=-6; Base only 0.05x as ")),
         "probable as methylation model"),
       pch=c(18,18,18,NA,NA,NA,NA,NA,NA,NA,NA,NA),lty=c(NA,NA,NA,NA,1,1,1,1,NA,1,NA,2,NA),
       lwd=c(NA,NA,NA,NA,3,3,3,3,NA,3,NA,3,NA),
       col=c(col_gene[1:3],NA,col_gene,NA,"grey",NA,"grey",NA),cex=cex_legend)
par(xpd=F)
dev.off()


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
sample_annot[19,2] <- "111N_2057T_M3"
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

# ccRCC
ccRCC_annot <- read.xlsx(paste(output_dir,"Klinische data location III CRC, ccRCC, melanoma- compleet 17022022.xlsx",sep=""),
                         2,startRow=1,colNames=T)
ccRCC_annot$Sample_ID <- paste(ccRCC_annot[,1],ccRCC_annot[,2],sep="")

sum(sample_annot$ccRCC%in%ccRCC_annot$Sample_ID)

# melanoma
melanoma_annot <- read.xlsx(paste(output_dir,"Klinische data location III CRC, ccRCC, melanoma- compleet 17022022_edited.xlsx",sep=""),
                         3,startRow=1,colNames=T)
melanoma_annot[6,2] <- "P32" 
melanoma_annot$Sample_ID <- melanoma_annot[,2]

melanoma_annot_2 <- read.xlsx(paste(output_dir,"Clinical follow-up data additional melanoma samples - Location-III.xlsx",sep=""),1,colNames=T)
melanoma_annot_2[,3] <- tolower(melanoma_annot_2[,3])
melanoma_annot_2 <- cbind(melanoma_annot_2[,1:8],rep(NA,nrow(melanoma_annot_2)),melanoma_annot_2[,9:ncol(melanoma_annot_2)],melanoma_annot_2[,2])
melanoma_annot_2$Date.last.contact[is.na(melanoma_annot_2$Date.last.contact)] <- "?"
colnames(melanoma_annot_2) <- colnames(melanoma_annot)

melanoma_annot <- rbind(melanoma_annot,melanoma_annot_2)

sum(sample_annot$melanoma%in%melanoma_annot$Sample_ID)

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
      
    
    sample_annot_amplicon_i <- sample_annot[sample_annot$Sample_ID%in%sample_amplicon_i,]
    
    # Check if annotation has right length
    if (length(files_amplicon_i)!=length(sample_annot_amplicon_i[,1])){
      stop("colData has wrong dimensions")
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

suff_cov <- 10
perc_suff_cov <- 0.5
mask_cov <- rowSums(totalReads(allrrbs)[,col_mask_CRC]>suff_cov)>perc_suff_cov*sum(col_mask_CRC) | 
  rowSums(totalReads(allrrbs)[,col_mask_ccRCC]>suff_cov)>perc_suff_cov*sum(col_mask_ccRCC) | 
  rowSums(totalReads(allrrbs)[,col_mask_melanoma]>suff_cov)>perc_suff_cov*sum(col_mask_melanoma)
sum(mask_cov)

missing_val <- 0
perc_missig_val <- 0.9
mask_missing <- rowSums(totalReads(allrrbs)[,col_mask_CRC]>missing_val)>perc_missig_val*sum(col_mask_CRC) | 
  rowSums(totalReads(allrrbs)[,col_mask_ccRCC]>missing_val)>perc_missig_val*sum(col_mask_ccRCC) | 
  rowSums(totalReads(allrrbs)[,col_mask_melanoma]>missing_val)>perc_missig_val*sum(col_mask_melanoma)
sum(mask_missing)

mask_allrrbs <- mask_missing & mask_cov
sum(mask_allrrbs)
rrbs <- allrrbs[mask_allrrbs]

amplicons_allrrbs <- rep(amplicons,amplicon_allrrbs_rows)
amplicons_rrbs <- rep(amplicons,amplicon_allrrbs_rows)[mask_allrrbs]

# Set working directory
setwd(output_dir)

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


#####################################
# Check data in regions of interest #
#####################################

locs <- rowRanges(rrbs)
locs$name <- amplicons_rrbs

locs_allrrbs <- rowRanges(allrrbs)
locs_allrrbs$name <- amplicons_allrrbs

table(locs$name)
table(locs_allrrbs$name)
names(table(locs_allrrbs$name))[!(names(table(locs_allrrbs$name))%in%names(table(locs$name)))]

table(data.frame(locs)[,1])
# CRC: NDRG4 (chr16:58462846-58513628)
head(data.frame(locs)[data.frame(locs)[,1]==16,])
# ccRCC: GREM1 (chr15:32718004-32745106)
head(data.frame(locs)[data.frame(locs)[,1]==15,])
# melanoma: LY75 (chr2:15980335-159904756)
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

CRC_annot <- CRC_annot[CRC_annot$Sample_ID!="178T",]

CRC_annot_sort <- CRC_annot[sort(CRC_annot$Sample_ID,index.return=T)$ix,]
CRC_annot_sort <- CRC_annot_sort[CRC_annot_sort$Sample_ID%in%colnames(betas_NDRG4),]
betas_NDRG4_sort <- betas_NDRG4[,sort(colnames(betas_NDRG4),index.return=T)$ix]
betas_NDRG4_sort <- betas_NDRG4_sort[,colnames(betas_NDRG4_sort)%in%CRC_annot_sort$Sample_ID]
sum(colnames(betas_NDRG4_sort)==CRC_annot_sort$Sample_ID)
write.table(CRC_annot_sort,file="CRC_annot_sort.txt",col.names = T,row.names = F,quote = F,sep="\t")

# Correction for age and sex
colnames(CRC_annot_sort) <- c("Sample","Tissue","Gender","Age","Meth","Sample_ID")
CRC_annot_sort$Gender <- as.factor(CRC_annot_sort$Gender)
CRC_annot_sort$Meth <- as.factor(CRC_annot_sort$Meth)
CRC_annot_sort$Tissue <- as.factor(CRC_annot_sort$Tissue)
print(CRC_annot_sort$Gender)
print(CRC_annot_sort$Age)
print(CRC_annot_sort$Tissue)
print(CRC_annot_sort$Meth)

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
formula_CF <- formula(paste("Tissue ~ Gender+Age",sep=""))
cn_formula_CF <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_CF))),"+",fixed=TRUE)))
model_CF <- glm(formula_CF, 
                data = CRC_annot_sort_BS, 
                family = binomial(link="logit"))
print(summary(model_CF))

dat_NDRG4_ROC_logit_perCpG$logit <- NA
dat_NDRG4_ROC_logit_perCpG$AIC <- NA
dat_NDRG4_ROC_logit_perCpG$AIC_CF <- AIC(model_CF)
dat_NDRG4_ROC_logit_perCpG$est_logit <- NA
dat_NDRG4_ROC_logit_perCpG$delta_AIC <- NA
for (i in 1:nrow(Mval_NDRG4_sort)){
  CRC_annot_sort_BS[paste("Meth_CpG_",i,sep="")] <- Mval_NDRG4_sort[i,]
  
  # logit 
  formula_i <- formula(paste("Tissue~Meth_CpG_",i,sep=""))
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
sum(colnames(betas_LY75_sort)==melanoma_annot_sort$Sample_ID)
write.table(melanoma_annot_sort,file="melanoma_annot_sort.txt",col.names = T,row.names = F,quote = F,sep="\t")


# Correctie in model voor leeftijd, geslacht en Breslow
colnames(melanoma_annot_sort) <- c("Study_Nr","Sample","Gender","Age","Location","Date_diag",
                                   "Breslow","Ulceration","T","Date_last_followup","Outcome","Meth","Sample_ID")
melanoma_annot_sort$Gender <- toupper(melanoma_annot_sort$Gender)
melanoma_annot_sort$Gender <- as.factor(melanoma_annot_sort$Gender)
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
formula_CF <- formula(paste("Surv(Followup_time, MK_death) ~ Gender+Age",sep=""))
cn_formula_CF <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_CF))),"+",fixed=TRUE)))
cox_CF <- coxph(formula_CF, data = melanoma_annot_sort_BS)
print(summary(cox_CF))

dat_LY75_ROC_logit_perCpG$logit <- NA
dat_LY75_ROC_logit_perCpG$AIC <- NA
dat_LY75_ROC_logit_perCpG$AIC_CF <- AIC(cox_CF)
dat_LY75_ROC_logit_perCpG$est_logit <- NA
dat_LY75_ROC_logit_perCpG$delta_AIC <- NA
for (i in 1:nrow(Mval_LY75_sort)){
  melanoma_annot_sort_BS[paste("Meth_CpG_",i,sep="")] <- Mval_LY75_sort[i,]
  
  # cox regression
  # formula_i <- formula(paste("Surv(Followup_time, MK_death) ~ Meth_CpG_",i,"+Gender+Age+Breslow",sep=""))
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


dat_GREM1_ROC_logit_perCpG$logit <- NA
dat_GREM1_ROC_logit_perCpG$AIC <- NA
dat_GREM1_ROC_logit_perCpG$AIC_CF <- AIC(cox_CF)
dat_GREM1_ROC_logit_perCpG$est_logit <- NA
dat_GREM1_ROC_logit_perCpG$delta_AIC <- NA
dat_GREM1_ROC_logit_perCpG$prob <- NA
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
  
  dat_GREM1_ROC_logit_perCpG$AIC[i] <- AIC(cox_BS_i)
  dat_GREM1_ROC_logit_perCpG$delta_AIC[i] <- dat_GREM1_ROC_logit_perCpG$AIC[i] - dat_GREM1_ROC_logit_perCpG$AIC_CF[i]
  if (sign(dat_GREM1_ROC_logit_perCpG$AIC[i] - dat_GREM1_ROC_logit_perCpG$AIC_CF[i])==1){
    dat_GREM1_ROC_logit_perCpG$prob[i] <- exp((dat_GREM1_ROC_logit_perCpG$AIC[i] - dat_GREM1_ROC_logit_perCpG$AIC_CF[i])/2)
  } else {
    dat_GREM1_ROC_logit_perCpG$prob[i] <- exp((dat_GREM1_ROC_logit_perCpG$AIC[i] - dat_GREM1_ROC_logit_perCpG$AIC_CF[i])/2)
  }
  
}              
write.table(dat_GREM1_ROC_logit_perCpG,"GREM1_AIC_table.txt",col.names = T,row.names = F,quote = F, sep="\t")


####################
# Overview objects #
####################

head(dat_NDRG4_ROC_logit_perCpG)
head(dat_GREM1_ROC_logit_perCpG)
head(dat_LY75_ROC_logit_perCpG)

################################
# Plot CRC (NDRG4): DIAGNOSTIC #
################################

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
              max_y=23)
dev.off()


####################################
# Plot melanoma (LY75): PROGNOSTIC #
####################################

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
              max_y=16)
dev.off()



##################################
# Plot ccRCC (GREM1): PROGNOSTIC #
##################################

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
          max_y=25)
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


##############################
# Intra amplicon differences #
##############################

if(!dir.exists("AIC_plots")){dir.create("AIC_plots")}
if(!dir.exists("AIC_plots/CpG_plots")){dir.create("AIC_plots/CpG_plots")}

# Functions
calculate_delta_AIC_diagnostic <- function(dat){
  CpG_names <- colnames(dat)[grepl("CpG",colnames(dat))]
  CpG_names <- gsub("^CpG_","",CpG_names)
  
  # first CpG
  formula_1 <- formula(paste("Tissue ~ CpG_",CpG_names[1],sep=""))
  cn_formula_1 <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_1))),"+",fixed=TRUE)))
  logreg_BS_1 <- glm(formula_1, data = dat, family = binomial(link="logit"))
  # print(summary(cox_BS_1))
  
  AIC_1 <- AIC(logreg_BS_1)
  
  # second CpG
  formula_2 <- formula(paste("Tissue ~ CpG_",CpG_names[2],sep=""))
  cn_formula_2 <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_2))),"+",fixed=TRUE)))
  logreg_BS_2 <- glm(formula_2, data = dat, family = binomial(link="logit"))
  # print(summary(cox_BS_2))
  
  AIC_2 <- AIC(logreg_BS_2)
  
  delta_AIC <- AIC_2 - AIC_1
  return(c(AIC_1,AIC_2,delta_AIC))
}

calculate_delta_AIC_prognostic <- function(dat){
  CpG_names <- colnames(dat)[grepl("CpG",colnames(dat))]
  CpG_names <- gsub("^CpG_","",CpG_names)
  
  # first CpG
  formula_1 <- formula(paste("Surv(Followup_time, NK_death) ~ CpG_",CpG_names[1],sep=""))
  cn_formula_1 <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_1))),"+",fixed=TRUE)))
  cox_BS_1 <- coxph(formula_1, data = dat) 
  # print(summary(cox_BS_1))
  
  AIC_1 <- AIC(cox_BS_1)
  
  # second CpG
  formula_2 <- formula(paste("Surv(Followup_time, NK_death) ~ CpG_",CpG_names[2],sep=""))
  cn_formula_2 <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_2))),"+",fixed=TRUE)))
  cox_BS_2 <- coxph(formula_2, data = dat) 
  # print(summary(cox_BS_2))
  
  AIC_2 <- AIC(cox_BS_2)
  
  delta_AIC <- AIC_2 - AIC_1
  return(c(AIC_1,AIC_2,delta_AIC))
}

bootstrap_AIC_diagnostic <- function(l){
  dat_test_shuffled_l <- dat_test_i
  dat_test_shuffled_l[,paste0("CpG_",cn_names[i])] <- sample(dat_test_shuffled_l[,paste0("CpG_",cn_names[i])])
  dat_test_shuffled_l[,paste0("CpG_",cn_names[j])] <- sample(dat_test_shuffled_l[,paste0("CpG_",cn_names[j])])
  return(calculate_delta_AIC_diagnostic(dat_test_shuffled_l))
}

bootstrap_AIC_prognostic <- function(l){
  dat_test_shuffled_l <- dat_test_i
  dat_test_shuffled_l[,paste0("CpG_",cn_names[i])] <- sample(dat_test_shuffled_l[,paste0("CpG_",cn_names[i])])
  dat_test_shuffled_l[,paste0("CpG_",cn_names[j])] <- sample(dat_test_shuffled_l[,paste0("CpG_",cn_names[j])])
  return(calculate_delta_AIC_prognostic(dat_test_shuffled_l))
}

# all comparisons
cex_AIC_plots <- 1.25
seed_perm <- 42

#NDRG4
set.seed(seed_perm)
pval_list_NDRG4 <- list()
for (k in 1:nrow(primer_annot_NDRG4)){
  print(paste0("Primer: ",primer_annot_NDRG4$name[k]))
  print(paste0("Location: chr ",primer_annot_NDRG4$chromosome[k]," from ",primer_annot_NDRG4$begin[k]," to ",primer_annot_NDRG4$end[k]))
  mask_amplicon <- grepl(paste0(primer_annot_NDRG4$name[k],"$"),rownames(Mval_NDRG4_sort))
  dat_amplicon <- t(Mval_NDRG4_sort[mask_amplicon,])
  
  dim(dat_amplicon)
  sum(CRC_annot_sort$Sample_ID==rownames(dat_amplicon))
  
  dat_amplicon <- cbind(CRC_annot_sort[,colnames(CRC_annot_sort)%in%c("Tissue")],dat_amplicon)
  colnames(dat_amplicon) <- c("Tissue", paste0(rep("CpG_",sum(mask_amplicon)),rownames(Mval_NDRG4_sort[mask_amplicon,])))
  
  cn_names <- colnames(dat_amplicon)[grepl("CpG",colnames(dat_amplicon))]
  cn_names <- gsub("CpG_","",cn_names)
  
  print(paste0(length(cn_names)," CpGs"))
  
  if (length(cn_names)==0){
    pval_list_NDRG4[[k]] <- NULL
    next()
  }
  
  dat_amplicon <- data.frame(dat_amplicon,stringsAsFactors = F)
  dat_amplicon$Tissue <- as.factor(dat_amplicon$Tissue)
  
  head(dat_amplicon)
  
  AIC_amplicon <- NULL
  for (i in 1:length(cn_names)){
    dat_test_i <- dat_amplicon[,c("Tissue",paste0("CpG_",cn_names[i]))]
    
    CpG_name_i <- colnames(dat_test_i)[grepl("CpG",colnames(dat_test_i))]
    CpG_name_i <- gsub("^CpG_","",CpG_name_i)
    
    formula_i <- formula(paste("Tissue ~ CpG_",CpG_name_i,sep=""))
    cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
    logreg_BS_i <- glm(formula_i, data = dat_test_i, family = binomial(link="logit"))
    
    AIC_i <- AIC(logreg_BS_i)
    AIC_amplicon <- rbind(AIC_amplicon,c(CpG_name_i,AIC_i))
  }
  AIC_amplicon <- data.frame(AIC_amplicon, stringsAsFactors = F)
  colnames(AIC_amplicon) <- c("CpG","AIC")
  AIC_amplicon$AIC <- as.numeric(AIC_amplicon$AIC)
  
  
  combos <- ncol(combn(1:length(cn_names),2))
  pval <- data.frame(CpG_1=rep("",combos),
                     CpG_2=rep("",combos),
                     AIC_1=rep("",combos),
                     AIC_2=rep("",combos),
                     delta_AIC=rep("",combos),
                     pval=rep(1,combos),
                     stringsAsFactors = F)
  counter <- 1
  for (i in 1:(length(cn_names)-1)){
    for (j in (i+1):length(cn_names)){
      print(paste0("CpG comparison: ",paste(cn_names[c(i,j)],collapse=" vs ")))
      
      pval[counter,1:2] <- cn_names[c(i,j)]
      
      dat_test_i <- dat_amplicon[,c("Tissue",paste0(rep("CpG_",2),cn_names[c(i,j)]))]
      pval[counter,3:5] <- calculate_delta_AIC_diagnostic(dat_test_i)

      distr_AIC <- mclapply(1:10000,bootstrap_AIC_diagnostic,mc.cores=40)
      distr_AIC <- unlist(distr_AIC)

      jpeg(paste0("AIC_plots/CpG_plots/CRC_CpG_",cn_names[i],"vs",cn_names[j],".jpg"))
      hist(distr_AIC,breaks=50)
      abline(v=pval[counter,5],col="red")
      dev.off()
      pval[counter,6] <- 2*min(sum(pval[counter,5]>distr_AIC)/length(distr_AIC),
                               sum(pval[counter,5]<distr_AIC)/length(distr_AIC))
      counter <- counter + 1
    }
  } 
  pval$fdr <- p.adjust(pval$pval)
  pval_list_NDRG4[[k]] <- pval 
  
  if (sum(pval$fdr<0.05)>0){
    pval_sign <- pval[pval$fdr<0.05,]
    print(pval_sign)
    sign_level <- (pval_sign$fdr < 0.05) + (pval_sign$fdr < 0.005) + (pval_sign$fdr < 0.0005) + (pval_sign$fdr < 0.00005)
    
    pos_x <-as.numeric(gsub("^.+_","",gsub("_[NDRG|LY|GREM].+$","",AIC_amplicon$CpG)))
    
    jpeg(paste0("AIC_plots/",primer_annot_NDRG4$name[k],".jpg"))
    plot(pos_x,AIC_amplicon$AIC,type="l") # this line is just to get limits for plots in next command
    add_space <- (par("usr")[4]-par("usr")[3])*0.05*(nrow(pval_sign)+1)
    # print(add_space)
    plot(pos_x,AIC_amplicon$AIC,type="l",main=primer_annot_NDRG4$name[k],
         ylim=c((par("usr")[3]-add_space),par("usr")[4]),xlab="Pos",ylab="AIC",cex.lab=cex_AIC_plots, cex.axis=cex_AIC_plots)
    increment_space <- add_space/(nrow(pval_sign)+1)
    for (t in 1:nrow(pval_sign)){
      segments(pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])],par("usr")[3]+increment_space*t,
               pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])],par("usr")[3]+increment_space*t)
      mid <- (pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])] +
                pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])])/2
      text(mid,par("usr")[3]+increment_space*t-(0.33*increment_space),paste(rep("*",sign_level[t]),collapse=""))
    }
    points(pos_x,AIC_amplicon$AIC,pch=16)
    dev.off()
  } else {
    print("No significant differences!")
  }
}
names(pval_list_NDRG4) <- primer_annot_NDRG4$name
saveRDS(pval_list_NDRG4,paste0("pval_list_NDRG4_seed",seed_perm,".rds"))

#GREM1
set.seed(seed_perm)
pval_list_GREM1 <- list()
for (k in 1:nrow(primer_annot_GREM1)){
  print(paste0("Primer: ",primer_annot_GREM1$name[k]))
  print(paste0("Location: chr ",primer_annot_GREM1$chromosome[k]," from ",primer_annot_GREM1$begin[k]," to ",primer_annot_GREM1$end[k]))
  mask_amplicon <- grepl(paste0(primer_annot_GREM1$name[k],"$"),rownames(Mval_GREM1_sort))
  dat_amplicon <- t(Mval_GREM1_sort[mask_amplicon,])
  
  dim(dat_amplicon)
  sum(ccRCC_annot_sort$Sample_ID==rownames(dat_amplicon))
  
  dat_amplicon <- cbind(ccRCC_annot_sort[,colnames(ccRCC_annot_sort)%in%c("Followup_time","NK_death")],dat_amplicon)
  colnames(dat_amplicon) <- c(colnames(ccRCC_annot_sort)[colnames(ccRCC_annot_sort)%in%c("Followup_time","NK_death")],
                              paste0(rep("CpG_",sum(mask_amplicon)),rownames(Mval_GREM1_sort[mask_amplicon,])))
  
  cn_names <- colnames(dat_amplicon)[grepl("CpG",colnames(dat_amplicon))]
  cn_names <- gsub("CpG_","",cn_names)
  
  print(paste0(length(cn_names)," CpGs"))
  
  if (length(cn_names)==0){
    pval_list_GREM1[[k]] <- "No CpGs"
    next()
  }
  
  dat_amplicon <- data.frame(dat_amplicon,stringsAsFactors = F)
  
  head(dat_amplicon)
  
  AIC_amplicon <- NULL
  for (i in 1:length(cn_names)){
    dat_test_i <- dat_amplicon[,c("Followup_time","NK_death",paste0("CpG_",cn_names[i]))]
    
    CpG_name_i <- colnames(dat_test_i)[grepl("CpG",colnames(dat_test_i))]
    CpG_name_i <- gsub("^CpG_","",CpG_name_i)
    
    formula_i <- formula(paste("Surv(Followup_time, NK_death) ~ CpG_",CpG_name_i,sep=""))
    cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
    cox_BS_i <- coxph(formula_i, data = dat_test_i)

    AIC_i <- AIC(cox_BS_i)
    AIC_amplicon <- rbind(AIC_amplicon,c(CpG_name_i,AIC_i))
  }
  AIC_amplicon <- data.frame(AIC_amplicon, stringsAsFactors = F)
  colnames(AIC_amplicon) <- c("CpG","AIC")
  AIC_amplicon$AIC <- as.numeric(AIC_amplicon$AIC)
  
  
  combos <- ncol(combn(1:length(cn_names),2))
  pval <- data.frame(CpG_1=rep("",combos),
                     CpG_2=rep("",combos),
                     AIC_1=rep("",combos),
                     AIC_2=rep("",combos),
                     delta_AIC=rep("",combos),
                     pval=rep(1,combos),
                     stringsAsFactors = F)
  counter <- 1
  for (i in 1:(length(cn_names)-1)){
    for (j in (i+1):length(cn_names)){
      print(paste0("CpG comparison: ",paste(cn_names[c(i,j)],collapse=" vs ")))
      
      pval[counter,1:2] <- cn_names[c(i,j)]
      
      dat_test_i <- dat_amplicon[,c("Followup_time","NK_death",paste0(rep("CpG_",2),cn_names[c(i,j)]))]
      pval[counter,3:5] <- calculate_delta_AIC_prognostic(dat_test_i)
      distr_AIC <- rep(0,1000)
      
      distr_AIC <- mclapply(1:10000,bootstrap_AIC_prognostic,mc.cores=40)
      distr_AIC <- unlist(distr_AIC)

      jpeg(paste0("AIC_plots/CpG_plots/ccRCC_CpG_",cn_names[i],"vs",cn_names[j],".jpg"))
      hist(distr_AIC,breaks=50)
      abline(v=pval[counter,5],col="red")
      dev.off()
      pval[counter,6] <- min(sum(pval[counter,5]>distr_AIC)/length(distr_AIC),sum(pval[counter,5]<distr_AIC)/length(distr_AIC))
      counter <- counter + 1
    }
  }
  pval$fdr <- p.adjust(pval$pval)
  pval_list_GREM1[[k]] <- pval 
  
  if (sum(pval$fdr<0.05)>0){
    pval_sign <- pval[pval$fdr<0.05,]
    print(pval_sign)
    sign_level <- (pval_sign$fdr < 0.05) + (pval_sign$fdr < 0.005) + (pval_sign$fdr < 0.0005) + (pval_sign$fdr < 0.00005) 
    
    pos_x <-as.numeric(gsub("^.+_","",gsub("_[NDRG|LY|GREM].+$","",AIC_amplicon$CpG)))
    
    jpeg(paste0("AIC_plots/",primer_annot_GREM1$name[k],".jpg"))
    plot(pos_x,AIC_amplicon$AIC,type="l") # this line is just to get limits for plots in next command
    add_space <- (par("usr")[4]-par("usr")[3])*0.05*(nrow(pval_sign)+1)
    # print(add_space)
    plot(pos_x,AIC_amplicon$AIC,type="l",main=primer_annot_GREM1$name[k],
         ylim=c((par("usr")[3]-add_space),par("usr")[4]),xlab="Pos",ylab="AIC",cex.lab=cex_AIC_plots,cex.axis=cex_AIC_plots)
    increment_space <- add_space/(nrow(pval_sign)+1)
    for (t in 1:nrow(pval_sign)){
      segments(pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])],par("usr")[3]+increment_space*t,
               pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])],par("usr")[3]+increment_space*t)
      mid <- (pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])] +
                pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])])/2
      text(mid,par("usr")[3]+increment_space*t-(0.33*increment_space),paste(rep("*",sign_level[t]),collapse=""))
    }
    points(as.numeric(gsub("^.+_","",AIC_amplicon$CpG)),AIC_amplicon$AIC,pch=16)
    dev.off()
  } else {
    print("No significant differences!")
  }
}
names(pval_list_GREM1) <- primer_annot_GREM1$name
saveRDS(pval_list_GREM1,paste0("pval_list_GREM1_seed",seed_perm,".rds"))

#LY75
colnames(melanoma_annot_sort)[colnames(melanoma_annot_sort)=="MK_death"] <- "NK_death"
set.seed(seed_perm)
pval_list_LY75 <- list()
for (k in 1:nrow(primer_annot_LY75)){
  print(paste0("Primer: ",primer_annot_LY75$name[k]))
  print(paste0("Location: chr ",primer_annot_LY75$chromosome[k]," from ",primer_annot_LY75$begin[k]," to ",primer_annot_LY75$end[k]))
  mask_amplicon <- grepl(paste0(primer_annot_LY75$name[k],"$"),rownames(Mval_LY75_sort))
  dat_amplicon <- t(Mval_LY75_sort[mask_amplicon,])
  
  dim(dat_amplicon)
  sum(melanoma_annot_sort$Sample_ID==rownames(dat_amplicon))
  
  dat_amplicon <- cbind(melanoma_annot_sort[,colnames(melanoma_annot_sort)%in%c("Followup_time","NK_death")],dat_amplicon)
  colnames(dat_amplicon) <- c(colnames(melanoma_annot_sort)[colnames(melanoma_annot_sort)%in%c("Followup_time","NK_death")],
                              paste0(rep("CpG_",sum(mask_amplicon)),rownames(Mval_LY75_sort[mask_amplicon,])))
  
  cn_names <- colnames(dat_amplicon)[grepl("CpG",colnames(dat_amplicon))]
  cn_names <- gsub("CpG_","",cn_names)
  
  print(paste0(length(cn_names)," CpGs"))
  
  if (length(cn_names)==0){
    pval_list_LY75[[k]] <- "No CpGs"
    next()
  }
  
  dat_amplicon <- data.frame(dat_amplicon,stringsAsFactors = F)

  head(dat_amplicon)
  
  AIC_amplicon <- NULL
  for (i in 1:length(cn_names)){
    dat_test_i <- dat_amplicon[,c("Followup_time","NK_death",paste0("CpG_",cn_names[i]))]
    
    CpG_name_i <- colnames(dat_test_i)[grepl("CpG",colnames(dat_test_i))]
    CpG_name_i <- gsub("^CpG_","",CpG_name_i)
    
    formula_i <- formula(paste("Surv(Followup_time, NK_death) ~ CpG_",CpG_name_i,sep=""))
    cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
    cox_BS_i <- coxph(formula_i, data = dat_test_i)
    
    AIC_i <- AIC(cox_BS_i)
    AIC_amplicon <- rbind(AIC_amplicon,c(CpG_name_i,AIC_i))
  }
  AIC_amplicon <- data.frame(AIC_amplicon, stringsAsFactors = F)
  colnames(AIC_amplicon) <- c("CpG","AIC")
  AIC_amplicon$AIC <- as.numeric(AIC_amplicon$AIC)
  
  
  combos <- ncol(combn(1:length(cn_names),2))
  pval <- data.frame(CpG_1=rep("",combos),
                     CpG_2=rep("",combos),
                     AIC_1=rep("",combos),
                     AIC_2=rep("",combos),
                     delta_AIC=rep("",combos),
                     pval=rep(1,combos),
                     stringsAsFactors = F)
  counter <- 1
  for (i in 1:(length(cn_names)-1)){
    for (j in (i+1):length(cn_names)){
      print(paste0("CpG comparison: ",paste(cn_names[c(i,j)],collapse=" vs ")))

      pval[counter,1:2] <- cn_names[c(i,j)]
      
      dat_test_i <- dat_amplicon[,c("Followup_time","NK_death",paste0(rep("CpG_",2),cn_names[c(i,j)]))]
      pval[counter,3:5] <- calculate_delta_AIC_prognostic(dat_test_i)
      distr_AIC <- rep(0,1000)
      
      distr_AIC <- mclapply(1:10000,bootstrap_AIC_prognostic,mc.cores=40)
      distr_AIC <- unlist(distr_AIC)

      jpeg(paste0("AIC_plots/CpG_plots/melanoma_CpG_",cn_names[i],"vs",cn_names[j],".jpg"))
      hist(distr_AIC,breaks=50)
      abline(v=pval[counter,5],col="red")
      dev.off()
      pval[counter,6] <- min(sum(pval[counter,5]>distr_AIC)/length(distr_AIC),
                             sum(pval[counter,5]<distr_AIC)/length(distr_AIC))
      counter <- counter + 1
    }
  }
  pval$fdr <- p.adjust(pval$pval)
  pval_list_LY75[[k]] <- pval 
  
  if (sum(pval$fdr<0.05)>0){
    pval_sign <- pval[pval$fdr<0.05,]
    print(pval_sign)
    sign_level <- (pval_sign$fdr < 0.05) + (pval_sign$fdr < 0.005) + (pval_sign$fdr < 0.0005) + (pval_sign$fdr < 0.00005) 
    
    pos_x <-as.numeric(gsub("^.+_","",gsub("_[NDRG|LY|GREM].+$","",AIC_amplicon$CpG)))
    
    jpeg(paste0("AIC_plots/",primer_annot_LY75$name[k],".jpg"))
    plot(pos_x,AIC_amplicon$AIC,type="l") # this line is just to get limits for plots in next command
    add_space <- (par("usr")[4]-par("usr")[3])*0.05*(nrow(pval_sign)+1)
    # print(add_space)
    plot(pos_x,AIC_amplicon$AIC,type="l",main=primer_annot_LY75$name[k],
         ylim=c((par("usr")[3]-add_space),par("usr")[4]),xlab="Pos",ylab="AIC",cex.lab=cex_AIC_plots,cex.axis=cex_AIC_plots)
    increment_space <- add_space/(nrow(pval_sign)+1)
    for (t in 1:nrow(pval_sign)){
      segments(pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])],par("usr")[3]+increment_space*t,
               pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])],par("usr")[3]+increment_space*t)
      mid <- (pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])] +
                pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])])/2
      text(mid,par("usr")[3]+increment_space*t-(0.33*increment_space),paste(rep("*",sign_level[t]),collapse=""))
    }
    points(as.numeric(gsub("^.+_","",AIC_amplicon$CpG)),AIC_amplicon$AIC,pch=16)
    dev.off()
  } else {
    print("No significant differences!")
  }
}
names(pval_list_LY75) <- primer_annot_LY75$name
saveRDS(pval_list_LY75,paste0("pval_list_LY75_seed",seed_perm,".rds"))




# Average delta AIC table

avAIC_table <- NULL
for (k in 1:nrow(primer_annot_GREM1)){
  mask_amplicon <- grepl(primer_annot_GREM1$name[k],rownames(Mval_GREM1_sort))
  avAIC_table <- rbind(avAIC_table,
                       c(primer_annot_GREM1$name[k],sum(mask_amplicon),
                         mean(dat_GREM1_ROC_logit_perCpG$delta_AIC[mask_amplicon])))
}
for (k in 1:nrow(primer_annot_LY75)){
  mask_amplicon <- grepl(primer_annot_LY75$name[k],rownames(Mval_LY75_sort))
  avAIC_table <- rbind(avAIC_table,
                       c(primer_annot_LY75$name[k],sum(mask_amplicon),
                         mean(dat_LY75_ROC_logit_perCpG$delta_AIC[mask_amplicon])))
}
for (k in 1:nrow(primer_annot_NDRG4)){
  mask_amplicon <- grepl(primer_annot_NDRG4$name[k],rownames(Mval_NDRG4_sort))
  avAIC_table <- rbind(avAIC_table,
                       c(primer_annot_NDRG4$name[k],sum(mask_amplicon),
                         mean(dat_NDRG4_ROC_logit_perCpG$delta_AIC[mask_amplicon])))
}
colnames(avAIC_table) <- c("AmpliconID","Am_CpG","avAIC")
avAIC_table
write.table(avAIC_table,"avAIC_table.txt",col.names = T,row.names = F,quote = F,sep="\t")
 

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

res_dupl_CpG <- NULL
#
dupl_CpG_GREM1 <- gsub("_[NDRG|GREM|LY].+$","",rownames(Mval_GREM1_sort))
dupl_CpG_GREM1 <- unique(dupl_CpG_GREM1[duplicated(dupl_CpG_GREM1)])
for (i in 1:length(dupl_CpG_GREM1)){
  dat_dupl_CpG_GREM1 <- Mval_GREM1_sort[grep(dupl_CpG_GREM1[i],rownames(Mval_GREM1_sort)),]
  beta_dupl_CpG_GREM1 <- betas_GREM1_sort[grep(dupl_CpG_GREM1[i],rownames(betas_GREM1_sort)),]
  
  res_dupl_CpG <- rbind(res_dupl_CpG,
                              c("GREM1",dupl_CpG_GREM1[i],
                                gsub(paste0(dupl_CpG_GREM1[i],"_"),"",rownames(dat_dupl_CpG_GREM1)),
                                t.test(dat_dupl_CpG_GREM1[1,],dat_dupl_CpG_GREM1[2,])$p.val,
                                abs(mean(beta_dupl_CpG_GREM1[1,],na.rm=T)-mean(beta_dupl_CpG_GREM1[2,],na.rm=T))))
}
dupl_CpG_LY75 <- gsub("_[NDRG|GREM|LY].+$","",rownames(Mval_LY75_sort))
dupl_CpG_LY75 <- unique(dupl_CpG_LY75[duplicated(dupl_CpG_LY75)])
for (i in 1:length(dupl_CpG_LY75)){
  dat_dupl_CpG_LY75 <- Mval_LY75_sort[grep(dupl_CpG_LY75[i],rownames(Mval_LY75_sort)),]
  beta_dupl_CpG_LY75 <- betas_LY75_sort[grep(dupl_CpG_LY75[i],rownames(betas_LY75_sort)),]
  
  res_dupl_CpG <- rbind(res_dupl_CpG,
                        c("LY75",dupl_CpG_LY75[i],
                          gsub(paste0(dupl_CpG_LY75[i],"_"),"",rownames(dat_dupl_CpG_LY75)),
                          t.test(dat_dupl_CpG_LY75[1,],dat_dupl_CpG_LY75[2,])$p.val,
                          abs(mean(beta_dupl_CpG_LY75[1,],na.rm=T)-mean(beta_dupl_CpG_LY75[2,],na.rm=T))))
}
dupl_CpG_NDRG4 <- gsub("_[NDRG|GREM|LY].+$","",rownames(Mval_NDRG4_sort))
dupl_CpG_NDRG4 <- unique(dupl_CpG_NDRG4[duplicated(dupl_CpG_NDRG4)])
for (i in 1:length(dupl_CpG_NDRG4)){
  dat_dupl_CpG_NDRG4 <- Mval_NDRG4_sort[grep(dupl_CpG_NDRG4[i],rownames(Mval_NDRG4_sort)),]
  beta_dupl_CpG_NDRG4 <- betas_NDRG4_sort[grep(dupl_CpG_NDRG4[i],rownames(betas_NDRG4_sort)),]
  
  res_dupl_CpG <- rbind(res_dupl_CpG,
                        c("NDRG4",dupl_CpG_NDRG4[i],
                          gsub(paste0(dupl_CpG_NDRG4[i],"_"),"",rownames(dat_dupl_CpG_NDRG4)),
                          t.test(dat_dupl_CpG_NDRG4[1,],dat_dupl_CpG_NDRG4[2,])$p.val,
                          abs(mean(beta_dupl_CpG_NDRG4[1,],na.rm=T)-mean(beta_dupl_CpG_NDRG4[2,],na.rm=T))))
}
res_dupl_CpG <- data.frame(res_dupl_CpG)
colnames(res_dupl_CpG) <- c("Gene","CpG","Amplicon1","Amplicon2","pval","diff_meth")
res_dupl_CpG$pval <- as.numeric(res_dupl_CpG$pval)
res_dupl_CpG$diff_meth <- as.numeric(res_dupl_CpG$diff_meth)
write.table(res_dupl_CpG,file="res_dupl_CpG.txt",col.names = T,row.names = F,quote = F,sep="\t")

