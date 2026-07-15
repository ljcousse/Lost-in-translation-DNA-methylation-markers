library(Rsamtools)
library(openxlsx)
library(parallel)

#########################
# Set working directory #
#########################

setwd("/data/homes/louisc/Project_KimSmits/mapping_bis")

#################
# Get bam files #
#################

files_list <- list.files()
files_list <- files_list[grepl("^run.*.bam$",files_list)]

######################
# Load amplicon data #
######################

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


print("LY75")
GR_LY75_pre <- data.frame(cbind(gsub("chr","",primer_annot_LY75$chromosome),
                                primer_annot_LY75$begin,
                                primer_annot_LY75$end,primer_annot_LY75$name))
colnames(GR_LY75_pre) <- c("chr","start","end","names") 
GR_LY75 <- GRanges(GR_LY75_pre)

print("GREM1")
GR_GREM1_pre <- data.frame(cbind(gsub("chr","",primer_annot_GREM1$chromosome),
                                 primer_annot_GREM1$begin,
                                 primer_annot_GREM1$end,primer_annot_GREM1$name))
colnames(GR_GREM1_pre) <- c("chr","start","end","names") 
GR_GREM1 <- GRanges(GR_GREM1_pre)

print("NDRG4")
GR_NDRG4_pre <- data.frame(cbind(gsub("chr","",primer_annot_NDRG4$chromosome),
                                 primer_annot_NDRG4$begin,
                                 primer_annot_NDRG4$end,primer_annot_NDRG4$name))
colnames(GR_NDRG4_pre) <- c("chr","start","end","names") 
GR_NDRG4 <- GRanges(GR_NDRG4_pre)

GR_amplicons <- c(GR_LY75,GR_GREM1,GR_NDRG4)
length(GR_amplicons)

GR_amplicons <- GR_amplicons[sort(data.frame(GR_amplicons)$width,index.return=T)$ix]

####################################
# split bamfiles for each amplicon #
####################################

in_amplicons <- NULL
length_limit_perc <- 120
myfun <- function(j){
# for (j in 1:length(files_list)){
  file_loc <- files_list[j]
  print(paste0(file_loc," (",j,"/",length(files_list),")"))
  
  bamFile <- BamFile(file_loc)
  aln <- scanBam(bamFile,asMates=TRUE)
  
  init_ix <- sort(aln[[1]]$qname,index.return=T)$ix
  
  if(sum(aln[[1]]$isize[init_ix][(c(1:(length(aln[[1]]$rname)/2))*2-1)]==(-aln[[1]]$isize[init_ix][(c(1:(length(aln[[1]]$rname)/2))*2)]))!=(length(aln[[1]]$qname)/2)){
    print("Insert sizes don't match")
  }
  
  if(sum(aln[[1]]$qname[init_ix][c(1:(length(aln[[1]]$qname)/2))*2]==aln[[1]]$qname[init_ix][c(1:(length(aln[[1]]$qname)/2))*2-1])!=(length(aln[[1]]$qname)/2)){
    print("Pairs are not following one another")
  }
  
  GR_aln_map <- GRanges(data.frame(seqnames = as.character(aln[[1]]$rname),
                                   start = aln[[1]]$pos,
                                   end = aln[[1]]$pos+abs(aln[[1]]$isize)))
  
  # dat_test <- cbind(aln[[1]]$qname[init_ix],data.frame(GR_aln_map)[init_ix,],aln[[1]]$isize)
  # dat_test[dat_test[,3]==160760991,]
  
  if(!dir.exists("/data/homes/louisc/Project_KimSmits/tabs_amplicon")){dir.create("/data/homes/louisc/Project_KimSmits/tabs_amplicon")}
  jpeg(paste0("/data/homes/louisc/Project_KimSmits/tabs_amplicon/",gsub("_S[0-9]+.+.bam$","",file_loc),"_readsize_distr.jpg"))
  hist(data.frame(GR_aln_map)$width,xlab="Read size (bp)",main=gsub("_S[0-9]+.+.bam$","",file_loc))
  dev.off()
  
  bp_overlap <- matrix(rep(0,length(GR_amplicons)*length(GR_aln_map)),ncol=length(GR_amplicons))
  for (i in 1:length(GR_amplicons)){
    mask_i <- (as.character(data.frame(GR_aln_map)$seqnames)==as.character(data.frame(GR_amplicons)$seqnames[i])) &
      ((data.frame(GR_aln_map)$start<=data.frame(GR_amplicons)$start[i] & data.frame(GR_aln_map)$end>=data.frame(GR_amplicons)$end[i]) |
         (data.frame(GR_aln_map)$start<=data.frame(GR_amplicons)$start[i] & data.frame(GR_aln_map)$end>=data.frame(GR_amplicons)$start[i]) |
         (data.frame(GR_aln_map)$start<=data.frame(GR_amplicons)$end[i] & data.frame(GR_aln_map)$end>=data.frame(GR_amplicons)$end[i]) |
         (data.frame(GR_aln_map)$start>=data.frame(GR_amplicons)$start[i] & data.frame(GR_aln_map)$end<=data.frame(GR_amplicons)$end[i])) & 
      (data.frame(GR_aln_map)$width<data.frame(GR_amplicons)$width[i]*length_limit_perc/100)
   
    bp_overlap[mask_i,i] <- apply(cbind(data.frame(GR_aln_map)$end[mask_i],rep(data.frame(GR_amplicons)$end[i],sum(mask_i))),1,min) -
      apply(cbind(data.frame(GR_aln_map)$start[mask_i],rep(data.frame(GR_amplicons)$start[i],sum(mask_i))),1,max)
  }

  dim(bp_overlap)
  bp_overlap <- bp_overlap[init_ix,]
  bp_overlap_pairs <- bp_overlap[c(1:(nrow(bp_overlap)/2))*2-1,]
  
  bp_overlap_pairs <- t(t(bp_overlap_pairs)/(data.frame(GR_amplicons)$width-1))
  # print(sum(bp_overlap_pairs>1))
  # print(sum(bp_overlap_pairs==1))
  # print(sum(c(rowSums(bp_overlap_pairs)==1)>1))
  
  which_amplicon_pairs <- rep(0,nrow(bp_overlap_pairs))
  mask_no_amplicon_pairs <- rowSums(bp_overlap_pairs==0)==ncol(bp_overlap_pairs)
  which_amplicon_pairs[!mask_no_amplicon_pairs] <- apply(bp_overlap_pairs,1,which.max)[!mask_no_amplicon_pairs]
  # head(which_amplicon,40)

  length(which_amplicon_pairs)
  which_amplicon <- rep(which_amplicon_pairs,each=2)
  head(which_amplicon,50)
  which_amplicon <- which_amplicon[sort(init_ix,index.return=T)$ix]
  head(which_amplicon,50)
  length(which_amplicon)
  
  GR_aln_map$names <- c("NAN",GR_amplicons$names)[which_amplicon+1]
  GR_aln_map
  
  tab_amplicon <- NULL
  for (i in 1:length(unique(GR_aln_map$names))){
    # print(unique(GR_aln_map$names)[i])
    mask_amplicon_i <- GR_aln_map$names==unique(GR_aln_map$names)[i]
    # print(sum(mask_amplicon_i))
    
    tab_amplicon <- rbind(tab_amplicon,c(unique(GR_aln_map$names)[i],sum(mask_amplicon_i)))
    
    qname_i <- aln[[1]]$qname[mask_amplicon_i]
    
    filter <- FilterRules(list(inAmplicon = function(x) x$qname %in% qname_i))
    
    filterBam(file_loc, paste0(unique(GR_aln_map$names)[i],"_",file_loc), filter=filter)
    
    if(sum(table(aln[[1]]$qname[mask_amplicon_i])==2)!=length(unique(aln[[1]]$qname[mask_amplicon_i]))){
      print(paste0("Not all reads are paired for ",unique(GR_aln_map$names)[i]))
    }
  }
  colnames(tab_amplicon) <- c("Amplicon","Mapped_reads")
  print(tab_amplicon)
  write.table(tab_amplicon,paste0("/data/homes/louisc/Project_KimSmits/tabs_amplicon/",gsub("_S[0-9]+.+.bam$","",file_loc),"_tab_amplicon.txt"),
              col.names = T, row.names = F, quote = F, sep="\t")
  
  print(paste0(as.numeric(tab_amplicon[tab_amplicon[,1]=="NAN",2])/sum(as.numeric(tab_amplicon[,2]))*100,
               "% of reads were result from aspecific hybridisation."))
  
  in_amplicons <- rbind(in_amplicons,c(file_loc,as.numeric(tab_amplicon[tab_amplicon[,1]=="NAN",2])/sum(as.numeric(tab_amplicon[,2]))*100))
  
  unlink(file_loc)
  unlink(paste0(file_loc,".bai"))
}
mclapply(1:length(files_list),myfun,mc.preschedule=T,mc.cores=20)
