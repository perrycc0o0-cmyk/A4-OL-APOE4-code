#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# Restored lightweight STEP4: regenerates the key tables/plots needed by later redraw scripts.
options(stringsAsFactors=FALSE); options(timeout=300)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(showtext);library(sysfonts)})
source(file.path("R", "load_config.R"))
ARIAL=A4OL_FONT_FILE; if(!file.exists(ARIAL))stop("Arial not found"); font_add("Arial",ARIAL); showtext_auto(); showtext_opts(dpi=300)
ROOT=file.path(A4OL_FIGURE_ROOT,"SCENIC_A4OL_FORMAL_V9"); OUT=file.path(ROOT,"SCENIC_EGR1_downstream_all"); PLOT=file.path(OUT,"plots"); TAB=file.path(OUT,"tables"); dir.create(PLOT,recursive=TRUE,showWarnings=FALSE); dir.create(TAB,recursive=TRUE,showWarnings=FALSE)
auc=as.matrix(readRDS(file.path(ROOT,"human/rds/human_auc_matrix.rds"))); cell=as.data.frame(readRDS(file.path(ROOT,"human/rds/human_cellInfo.rds"))); if("cell"%in%colnames(cell))rownames(cell)=cell$cell
stats=read.csv(file.path(ROOT,"human/tables/human_EGR1_regulon_AUC_stats.csv"),check.names=FALSE)
common=intersect(colnames(auc),rownames(cell)); regs=intersect(stats$regulon,rownames(auc))
human_egr1_auc_df=bind_rows(lapply(regs,function(reg)data.frame(species="Human",regulon=reg,cell=common,A4_status=as.character(cell[common,"A4_status"]),regulon_AUC=as.numeric(auc[reg,common]))))
write.csv(human_egr1_auc_df,file.path(TAB,"Human_EGR1_regulon_AUC_per_cell.csv"),row.names=FALSE)
UCELL=file.path(A4OL_FIGURE_ROOT,"UCell_A4OL","rds","human_olig_res04_with_A4_UCell.rds")
if(file.exists(UCELL) && requireNamespace("Seurat",quietly=TRUE)){
  u=readRDS(UCELL); meta=u@meta.data; cols=grep("UCell|A4",colnames(meta),value=TRUE); pref=cols[grepl("A4.*UCell|UCell.*A4",cols,ignore.case=TRUE)]; if(length(pref)==0)pref=cols[grepl("UCell",cols,ignore.case=TRUE)]; if(length(pref)>0){score_col=pref[1]; score=data.frame(cell=rownames(meta),A4_UCell_score=as.numeric(meta[[score_col]])); corr=left_join(human_egr1_auc_df,score,by="cell")%>%filter(is.finite(A4_UCell_score),is.finite(regulon_AUC)); write.csv(corr,file.path(TAB,"Human_EGR1_regulon_AUC_vs_A4_UCell_score_per_cell.csv"),row.names=FALSE); p=ggplot(corr,aes(A4_UCell_score,regulon_AUC,color=A4_status))+geom_point(size=.45,alpha=.35)+geom_smooth(method="lm",se=TRUE,linewidth=.55,color="black")+facet_wrap(~regulon,scales="free_y",ncol=2)+scale_color_manual(values=c("Other OLs"="#71A682","A4-OLs"="#D19246"))+theme_classic(base_family="Arial")+theme(legend.position="top",plot.title=element_text(face="bold",hjust=.5))+labs(title="Human EGR1 regulon activity correlates with A4 UCell score",x=paste0(score_col," score"),y="EGR1 regulon AUC"); ggsave(file.path(PLOT,"Human_EGR1_regulon_AUC_vs_A4_UCell_score.png"),p,width=8,height=5.2,dpi=300,bg="white")}
}
cat("DONE lightweight STEP4. For full targets/GO/network, run REDRAW_EGR1_NETWORKS_2x2_COMBINED_ARIAL.R and GO scripts.\n")
