
plotPlots_AIC <- function(dat_gene,locs_gene,gene,MSP_amplicon,primer_annot_gene,gene_cor=NULL,
                          gene_obj,transcript_ids,chrom,object_box,primer_mask,elevate_amplicon=matrix(c(10^6,10^6),nrow=1),max_y=NULL){
  chrom_gene <- as.character(data.frame(locs_gene)[1,1])

  col_gene <- c("cadetblue","blue","navy","lightgrey")

  index_sort <- sort(dat_gene$ID,index.return=T)$ix
  primer_mask <- primer_mask[index_sort]
  dat_gene <- dat_gene[index_sort,]
  
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

  pch_plot <- c(18,4)[primer_mask+1]
  print(pch_plot)
  
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
      points(dat_gene_plot_line_i$Pos,dat_gene_plot_line_i$delta_AIC,pch=pch_plot[dat_gene$Amplicon==primer_annot_gene$name[i]],col=col_gene[cols[i]])
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
  
  n_transcript_ids <- length(transcript_ids)
  print(n_transcript_ids)
  for (i in 1:n_transcript_ids){
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
      if(n_transcript_ids>10){
        text(dat_transcript$transcription_start_site[1]-((par("usr")[2]-par("usr")[1])/100*diff_units_transcript_id),i,
            labels=dat_transcript$ensembl_transcript_id[1],cex=cex_labels*0.6)
      } else {
        text(dat_transcript$transcription_start_site[1]-((par("usr")[2]-par("usr")[1])/100*diff_units_transcript_id),i,
            labels=dat_transcript$ensembl_transcript_id[1],cex=cex_labels)
      }
    } else {
      # print(par("usr"))
      x_poly <- c(dat_transcript$transcription_start_site[1],dat_transcript$transcription_start_site[1],
                  dat_transcript$transcription_start_site[1]-((par("usr")[2]-par("usr")[1])/100))
      y_poly <- c(i-0.4,i+0.4,i)
      # print(x_poly)
      # print(y_poly)
      polygon(x_poly,y_poly,col="green")
      if(n_transcript_ids>10){
          text(dat_transcript$transcription_start_site[1]+((par("usr")[2]-par("usr")[1])/100*diff_units_transcript_id*0.2),i,
            labels=dat_transcript$ensembl_transcript_id[1],cex=cex_labels*0.35)
      } else {
        text(dat_transcript$transcription_start_site[1]+((par("usr")[2]-par("usr")[1])/100*diff_units_transcript_id),i,
            labels=dat_transcript$ensembl_transcript_id[1],cex=cex_labels)
      }
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
