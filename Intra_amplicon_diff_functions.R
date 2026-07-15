#############
# Functions #
#############

# calculate_delta_AIC_diagnostic
###############

calculate_delta_AIC_diagnostic <- function(dat){
  CpG_names <- colnames(dat)[grepl("CpG",colnames(dat))]
  CpG_names <- gsub("^CpG_","",CpG_names)
  
  # first CpG
  formula_1 <- formula(paste("Tissue ~ CpG_",CpG_names[1],"+Patient_ID",sep=""))
  cn_formula_1 <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_1))),"+",fixed=TRUE)))
  logreg_BS_1 <- glm(formula_1, data = dat, family = binomial(link="logit"))
  # print(summary(cox_BS_1))
  
  AIC_1 <- AIC(logreg_BS_1)
  
  # second CpG
  formula_2 <- formula(paste("Tissue ~ CpG_",CpG_names[2],"+Patient_ID",sep=""))
  cn_formula_2 <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_2))),"+",fixed=TRUE)))
  logreg_BS_2 <- glm(formula_2, data = dat, family = binomial(link="logit"))
  # print(summary(cox_BS_2))
  
  AIC_2 <- AIC(logreg_BS_2)
  
  delta_AIC <- AIC_2 - AIC_1
  return(c(AIC_1,AIC_2,delta_AIC))
}

# calculate_delta_AIC_prognostic
###############

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

# bootstrap_AIC_diagnostic
###############

bootstrap_AIC_diagnostic <- function(i,j){
  i <- iter_i
  j <- iter_j
  dat_test_shuffled_l <- dat_test_i
  dat_test_shuffled_l[,paste0("CpG_",cn_names[i])] <- sample(dat_test_shuffled_l[,paste0("CpG_",cn_names[i])])
  dat_test_shuffled_l[,paste0("CpG_",cn_names[j])] <- sample(dat_test_shuffled_l[,paste0("CpG_",cn_names[j])])
  return(calculate_delta_AIC_diagnostic(dat_test_shuffled_l)[3])
}

# bootstrap_AIC_prognostic 
###############

bootstrap_AIC_prognostic <- function(l){
  i <- iter_i
  j <- iter_j
  dat_test_shuffled_l <- dat_test_i
  dat_test_shuffled_l[,paste0("CpG_",cn_names[i])] <- sample(dat_test_shuffled_l[,paste0("CpG_",cn_names[i])])
  dat_test_shuffled_l[,paste0("CpG_",cn_names[j])] <- sample(dat_test_shuffled_l[,paste0("CpG_",cn_names[j])])
  return(calculate_delta_AIC_prognostic(dat_test_shuffled_l)[3])
}

# perform AIC analysis
###############


perform_AIC_analysis <- function(primer_annot,Mval,sample_annot,gene,seed_perm=42){
  set.seed(seed_perm)
  pval_list <- list()
  counter_list <- 1
  for (k in 1:nrow(primer_annot)){
    print(paste0("Primer: ",primer_annot$name[k]))
    print(paste0("Location: chr ",primer_annot$chromosome[k]," from ",primer_annot$begin[k]," to ",primer_annot$end[k]))
    mask_amplicon <- grepl(paste0(primer_annot$name[k],"$"),rownames(Mval))
    dat_amplicon <- t(Mval[mask_amplicon,])

    if (gene=="NDRG4"){
      cn_oi <- c("Tissue","Patient_ID")
      part_formula <- "Tissue"
    } else {
      cn_oi <- c("NK_death","Followup_time")
      part_formula <- "Surv(Followup_time, NK_death)"
    }
    
    dim(dat_amplicon)
    sum(sample_annot$Sample_ID==rownames(dat_amplicon))
    
    dat_amplicon <- cbind(sample_annot[,colnames(sample_annot)%in%cn_oi],dat_amplicon)
    colnames(dat_amplicon) <- c(colnames(sample_annot[,colnames(sample_annot)%in%cn_oi]), paste0(rep("CpG_",sum(mask_amplicon)),rownames(Mval[mask_amplicon,])))
    
    cn_names <<- colnames(dat_amplicon)[grepl("CpG",colnames(dat_amplicon))]
    cn_names <<- gsub("CpG_","",cn_names)
    
    print(paste0(length(cn_names)," CpGs"))
    
    if (length(cn_names)==0){
      next()
    }
    
    dat_amplicon <- data.frame(dat_amplicon,stringsAsFactors = F)
    if(gene=="NDRG4"){
      dat_amplicon$Tissue <- as.factor(dat_amplicon$Tissue)
    }
    # print(head(dat_amplicon))
    # print(str(dat_amplicon))
    
    AIC_amplicon <- NULL
    for (i in 1:length(cn_names)){
      dat_test_i <- dat_amplicon[,c(cn_oi,paste0("CpG_",cn_names[i]))]
      
      CpG_name_i <- colnames(dat_test_i)[grepl("CpG",colnames(dat_test_i))]
      CpG_name_i <- gsub("^CpG_","",CpG_name_i)
      
      if (gene=="NDRG4"){
        formula_i <- formula(paste(part_formula," ~ CpG_",CpG_name_i,"+Patient_ID",sep=""))
      } else {
        formula_i <- formula(paste(part_formula," ~ CpG_",CpG_name_i,sep=""))
      }
      cn_formula_i <-  gsub(" ","",unlist(strsplit(gsub(".*~","",Reduce(paste, deparse(formula_i))),"+",fixed=TRUE)))
      if (gene=="NDRG4"){
        model_BS_i <- glm(formula_i, data = dat_test_i, family = binomial(link="logit"))
      } else {
        model_BS_i <- coxph(formula_i, data = dat_test_i)
      }
      
      AIC_i <- AIC(model_BS_i)
      AIC_amplicon <- rbind(AIC_amplicon,c(CpG_name_i,AIC_i))
    }
    AIC_amplicon <- data.frame(AIC_amplicon, stringsAsFactors = F)
    colnames(AIC_amplicon) <- c("CpG","AIC")
    AIC_amplicon$AIC <- as.numeric(AIC_amplicon$AIC)

    print(cn_names)
      
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
      print(sprintf("Comparing %s with %s CpGs",cn_names[i],length((i+1):length(cn_names))))
      for (j in (i+1):length(cn_names)){
        #print(paste0("CpG comparison: ",paste(cn_names[c(i,j)],collapse=" vs "))) 

        pval[counter,1:2] <- cn_names[c(i,j)]
        dat_test_i <<- dat_amplicon[,c(cn_oi,paste0(rep("CpG_",2),cn_names[c(i,j)]))]
        dat_test_i <- dat_amplicon[,c(cn_oi,paste0(rep("CpG_",2),cn_names[c(i,j)]))]

        #print(head(dat_test_i))

        iter_i <<- i
        iter_j <<- j

        if (gene=="NDRG4"){
          pval[counter,3:5] <- calculate_delta_AIC_diagnostic(dat_test_i)
          distr_AIC <- mclapply(1:10000,bootstrap_AIC_diagnostic,mc.cores=40)
        } else {
          pval[counter,3:5] <- calculate_delta_AIC_prognostic(dat_test_i)
          distr_AIC <- mclapply(1:10000,bootstrap_AIC_prognostic,mc.cores=40)
        }
        # print(pval[counter,3:5])
        # print(distr_AIC[[1]])
        distr_AIC <- unlist(distr_AIC)
        
        jpeg(paste0("AIC_plots/CpG_plots/CRC_CpG_",cn_names[i],"vs",cn_names[j],".jpg"))
        hist(distr_AIC,breaks=50)
        abline(v=pval[counter,5],col="red")
        dev.off()
        
        pval[counter,6] <- 2*min(sum(as.numeric(pval[counter,5])>distr_AIC)/length(distr_AIC),sum(as.numeric(pval[counter,5])<distr_AIC)/length(distr_AIC))
        counter <- counter + 1
      }
    } 
    pval$fdr <- p.adjust(pval$pval)
    pval_list[[counter_list]] <- pval 
    names(pval_list)[counter_list] <- primer_annot$name[k]
    counter_list <- counter_list + 1
    
    if (sum(pval$fdr<fdr_treshold_toplot)>0){
      pval_sign <- pval[pval$fdr<fdr_treshold_toplot,]
      print(pval_sign)
      if (fdr_treshold_toplot%%100/10==5){
        sign_level <<- (pval_sign$fdr < 0.05) + (pval_sign$fdr < 0.005) + (pval_sign$fdr < 0.0005) + (pval_sign$fdr < 0.0005)
      } else {
        sign_level <<- (pval_sign$fdr < 0.1) + (pval_sign$fdr < 0.01) + (pval_sign$fdr < 0.001) + (pval_sign$fdr < 0.001)
      }
      
      intra_amplicon_plot(primer_annot$name[k],AIC_amplicon,pval_sign,gene,k)

    } else {
      print("No significant differences!")
    }
  }
  return(pval_list)
}

intra_amplicon_plot <- function(name,AIC_amplicon,pval_sign,gene,k,dir="AIC_plots/"){
  pos_x <-as.numeric(gsub("^.+_","",gsub("_[NDRG|LY|GREM].+$","",AIC_amplicon$CpG)))

  jpeg(paste0(dir,name,".jpg"))

  plot(pos_x,AIC_amplicon$AIC,type="l") # this line is just to get limits for plots in next command
  add_space <- (par("usr")[4]-par("usr")[3])*0.05*(nrow(pval_sign)+1)
  # print(add_space)

  plot(pos_x,AIC_amplicon$AIC,type="l",main=name,
       ylim=c((par("usr")[3]-add_space),par("usr")[4]),
       xlim=c(primer_annot_list[[gene]]$begin[k]-n_bp_add,primer_annot_list[[gene]]$end[k]+n_bp_add),
       xlab="Pos",ylab="AIC",cex.lab=cex_AIC_plots, cex.axis=cex_AIC_plots)

  increment_space <- add_space/(nrow(pval_sign)+1)
  for (t in 1:nrow(pval_sign)){
    segments(pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])],par("usr")[3]+increment_space*t,
            pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])],par("usr")[3]+increment_space*t)
    mid <- (pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_1[t])] +
              pos_x[which(AIC_amplicon$CpG==pval_sign$CpG_2[t])])/2
    text(mid,par("usr")[3]+increment_space*t-(0.33*increment_space),paste(rep("*",sign_level[t]),collapse=""))
  }

  points(pos_x,AIC_amplicon$AIC,pch=16)

  rect(primer_annot_list[[gene]]$begin[k]-n_bp_add, (par("usr")[3]-add_space),
       primer_annot_list[[gene]]$begin[k]+primer_annot_list[[gene]]$len_fwd[k]+n_bp_add, par("usr")[4],
       col=alpha("cadetblue",alpha=0.20))

  rect(primer_annot_list[[gene]]$end[k]-primer_annot_list[[gene]]$len_rev[k]-n_bp_add, (par("usr")[3]-add_space),
       primer_annot_list[[gene]]$end[k]+n_bp_add, par("usr")[4],
       col=alpha("cadetblue",alpha=0.20))

  dev.off()
}

reduced_AIC_analysis <- function(gene,pval_list,pval_cutoff,dat_AIC){
  pval_list_red <- pval_list
  for (k in 1:nrow(primer_annot_list[[gene]])){
    if (is.null(pval_list_red[[primer_annot_list[[gene]]$name[k]]])){
       print(sprintf("%s skipped",primer_annot_list[[gene]]$name[k]))
       next()
    }
    dat_correl_k <- dat_correl[grepl(paste0(primer_annot_list[[gene]]$name[k],"$"),rownames(dat_correl)),]

    pval_list_red[[primer_annot_list[[gene]]$name[k]]] <- pval_list_red[[primer_annot_list[[gene]]$name[k]]][!((pval_list_red[[primer_annot_list[[gene]]$name[k]]]$CpG_1 %in% dat_correl_k$ID[(dat_correl_k$dist_CpG==0) & (!is.na(dat_correl_k$dist_CpG))]) | 
                                                            (pval_list_red[[primer_annot_list[[gene]]$name[k]]]$CpG_2 %in% dat_correl_k$ID[(dat_correl_k$dist_CpG==0) & (!is.na(dat_correl_k$dist_CpG))])),]
    # print(dim(pval_list_NDRG4_red[[k]]))
    pval <- pval_list_red[[primer_annot_list[[gene]]$name[k]]]

    if (sum(pval$pval<pval_cutoff)>0){
      pval_sign <- pval[pval$pval<pval_cutoff,]
      print(sprintf("%s significant p values for %s",nrow(pval_sign),primer_annot_list[[gene]]$name[k]))
      print(pval_sign[,c(1,2,5,6)])
      if (0.01%%100/10==5){
        sign_level <<- (pval_sign$pval < 0.05) + (pval_sign$pval < 0.005) + (pval_sign$pval < 0.0005) + (pval_sign$pval < 0.0005)
      } else {
        sign_level <<- (pval_sign$pval < 0.1) + (pval_sign$pval < 0.01) + (pval_sign$pval < 0.001) + (pval_sign$pval < 0.001)
      }
      
      AIC_amplicon <- dat_AIC[dat_AIC$ID%in%rownames(dat_correl_k),c("ID","AIC")]
      colnames(AIC_amplicon) <- c("CpG","AIC")
      intra_amplicon_plot(primer_annot_list[[gene]]$name[k],AIC_amplicon,pval_sign,gene,k,dir="AIC_plots_red/")
    } else {
      print(sprintf("No significant p values for %s",primer_annot_list[[gene]]$name[k]))
    }
  }
}