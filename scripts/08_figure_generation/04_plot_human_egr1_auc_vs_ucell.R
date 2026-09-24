#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE); options(timeout=300)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(showtext);library(sysfonts)})
source(file.path("R", "load_config.R"))
ARIAL_TTF=A4OL_FONT_FILE
if(!file.exists(ARIAL_TTF)) stop("Arial font file not found: ",ARIAL_TTF)
font_add("arial",ARIAL_TTF); showtext_auto(); showtext_opts(dpi=300)
ROOT=file.path(A4OL_FIGURE_ROOT,"SCENIC_A4OL_FORMAL_V9")
OUT=file.path(ROOT,"SCENIC_EGR1_downstream_all"); PLOT=file.path(OUT,"plots"); TAB=file.path(OUT,"tables")
dir.create(PLOT,recursive=TRUE,showWarnings=FALSE); dir.create(TAB,recursive=TRUE,showWarnings=FALSE)
input=file.path(TAB,"Human_EGR1_regulon_AUC_vs_A4_UCell_score_per_cell.csv")
if(!file.exists(input)) stop("Missing table: ",input,"\nRun STEP4 first.")
df=read.csv(input,check.names=FALSE)
need=c("regulon","A4_UCell_score","regulon_AUC","A4_status"); if(length(setdiff(need,colnames(df)))>0) stop("Missing columns: ",paste(setdiff(need,colnames(df)),collapse=", "))
df=df%>%mutate(regulon_label=case_when(grepl("^EGR1 \\(",regulon)~"EGR1 core regulon",grepl("^EGR1_extended \\(",regulon)~"EGR1 extended regulon",TRUE~regulon),
             A4_status=factor(A4_status,levels=c("A4-OLs","Other OLs")))%>%
  filter(regulon_label%in%c("EGR1 core regulon","EGR1 extended regulon"),is.finite(A4_UCell_score),is.finite(regulon_AUC),!is.na(A4_status))
write.csv(df,file.path(TAB,"Human_EGR1_regulon_AUC_vs_A4_UCell_score_per_cell_renamed.csv"),row.names=FALSE)
calc=function(d){p=tryCatch(cor.test(d$A4_UCell_score,d$regulon_AUC,method="spearman",exact=FALSE)$p.value,error=function(e)NA_real_);data.frame(n_cells=nrow(d),rho=cor(d$A4_UCell_score,d$regulon_AUC,method="spearman",use="complete.obs"),p=p)}
stat=bind_rows(calc(filter(df,regulon_label=="EGR1 core regulon"))%>%mutate(regulon_label="EGR1 core regulon"),
               calc(filter(df,regulon_label=="EGR1 extended regulon"))%>%mutate(regulon_label="EGR1 extended regulon"))%>%mutate(FDR=p.adjust(p,"BH"))
write.csv(stat,file.path(TAB,"Human_EGR1_core_extended_AUC_vs_A4_UCell_score_split_correlation_stats.csv"),row.names=FALSE)
cols=c("A4-OLs"="#D19246","Other OLs"="#71A682")
theme0=theme_classic(base_family="arial",base_size=13)+theme(plot.title=element_text(size=18,face="bold",hjust=.5),axis.title=element_text(size=16),axis.text=element_text(size=13,color="black"),legend.position="top",legend.title=element_blank(),legend.text=element_text(size=13),panel.border=element_rect(color="black",fill=NA,linewidth=.7),panel.grid=element_blank())
plot_one=function(reg,title,y,out){
 d=filter(df,regulon_label==reg); s=filter(stat,regulon_label==reg)
 lab=paste0("Spearman rho = ",sprintf("%.3f",s$rho),", FDR = ",format(s$FDR,scientific=TRUE,digits=2))
 p=ggplot(d,aes(A4_UCell_score,regulon_AUC,color=A4_status))+geom_point(size=.65,alpha=.35)+geom_smooth(aes(group=1),method="lm",se=TRUE,linewidth=.65,color="black",fill="grey75")+
   scale_color_manual(values=cols,breaks=c("A4-OLs","Other OLs"))+theme0+labs(title=title,x="A4_OL_UCell score",y=y)+annotate("text",x=Inf,y=Inf,label=lab,hjust=1.05,vjust=1.4,family="arial",size=4)
 ggsave(file.path(PLOT,paste0(out,".png")),p,width=6.6,height=5.4,dpi=300,bg="white"); ggsave(file.path(PLOT,paste0(out,".pdf")),p,width=6.6,height=5.4,bg="white")
}
plot_one("EGR1 core regulon","Human EGR1 core regulon activity correlates with A4 UCell score","EGR1 core regulon AUC","Human_EGR1_core_regulon_AUC_vs_A4_UCell_score")
plot_one("EGR1 extended regulon","Human EGR1 extended regulon activity correlates with A4 UCell score","EGR1 extended regulon AUC","Human_EGR1_extended_regulon_AUC_vs_A4_UCell_score")
