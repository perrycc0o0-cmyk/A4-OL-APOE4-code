#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE); options(timeout=300)
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(showtext);library(sysfonts)})
source(file.path("R", "load_config.R"))
ARIAL_TTF=A4OL_FONT_FILE
if(!file.exists(ARIAL_TTF)) stop("Arial font file not found: ",ARIAL_TTF)
font_add("arial",ARIAL_TTF); showtext_auto(); showtext_opts(dpi=300)
ROOT=file.path(A4OL_FIGURE_ROOT,"SCENIC_A4OL_FORMAL_V9")
OUTDIR=file.path(ROOT,"SCENIC_EGR1_downstream_all","plots_redraw"); TABLEDIR=file.path(ROOT,"SCENIC_EGR1_downstream_all","tables_redraw")
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE); dir.create(TABLEDIR,recursive=TRUE,showWarnings=FALSE)
COL_OTHER="#71A682"; COL_A4="#D19246"; COL_LINE="grey65"
detect=function(df,cand){h=cand[cand%in%colnames(df)]; if(length(h)==0) NA_character_ else h[1]}
get_regs=function(auc,species){regs=rownames(auc); if(species=="human"){core=grep("^EGR1 \\(",regs,value=TRUE)[1]; ext=grep("^EGR1_extended \\(",regs,value=TRUE)[1]}else{core=grep("^Egr1 \\(",regs,value=TRUE)[1]; ext=grep("^Egr1_extended \\(",regs,value=TRUE)[1]}; if(is.na(core)||is.na(ext))stop("Cannot detect regulons"); list(core=core,ext=ext,label_core=core,label_ext=ext)}
build=function(auc_file,cell_file,species){auc=readRDS(auc_file); cell=as.data.frame(readRDS(cell_file)); cell$cell=rownames(cell); common=intersect(colnames(auc),cell$cell); auc=auc[,common,drop=FALSE]; cell=filter(cell,cell%in%common)
 sample_col=detect(cell,c("sample","orig.ident","sample_id","donor","Sample","sampleID")); a4_col=detect(cell,c("A4_status","a4_status","group","Group","A4_group"))
 if(is.na(sample_col)||is.na(a4_col))stop("Cannot detect sample/A4 columns")
 r=get_regs(auc,species)
 df=bind_rows(data.frame(cell=colnames(auc),regulon_AUC=as.numeric(auc[r$core,]),regulon=r$label_core),data.frame(cell=colnames(auc),regulon_AUC=as.numeric(auc[r$ext,]),regulon=r$label_ext))%>%
   left_join(cell,by="cell")%>%transmute(cell=cell,sample=as.character(.data[[sample_col]]),A4_status=as.character(.data[[a4_col]]),regulon=regulon,regulon_AUC=regulon_AUC)%>%
   filter(A4_status%in%c("Other OLs","A4-OLs"),!is.na(sample),sample!="",is.finite(regulon_AUC))%>%
   group_by(sample,A4_status,regulon)%>%summarise(mean_auc=mean(regulon_AUC,na.rm=TRUE),n_cells=n(),.groups="drop")
 keep=df%>%group_by(sample,regulon)%>%summarise(n_groups=n_distinct(A4_status),.groups="drop")%>%filter(n_groups==2)
 df2=df%>%inner_join(select(keep,sample,regulon),by=c("sample","regulon"))%>%mutate(A4_status=factor(A4_status,levels=c("Other OLs","A4-OLs")),regulon=factor(regulon,levels=c(r$label_core,r$label_ext)))
 write.csv(df2,file.path(TABLEDIR,paste0("SampleLevel_",species,"_EGR1_regulon_AUC_paired_connected_only_plot_data.csv")),row.names=FALSE); df2}
plotp=function(df,species,out){p=ggplot(df,aes(A4_status,mean_auc,group=sample))+geom_line(color=COL_LINE,linewidth=.6,alpha=.75)+geom_point(aes(fill=A4_status),shape=21,size=3.8,stroke=.7,color="black")+facet_wrap(~regulon,nrow=1,scales="free_y")+scale_fill_manual(values=c("Other OLs"=COL_OTHER,"A4-OLs"=COL_A4))+labs(title=paste0(species," sample-level EGR1/Egr1 regulon activity"),x=NULL,y="Mean regulon AUC per sample")+theme_classic(base_family="arial")+theme(plot.title=element_text(size=24,face="plain",hjust=.5),axis.text.x=element_text(size=17,angle=25,hjust=1,color="black"),axis.text.y=element_text(size=17,color="black"),axis.title.y=element_text(size=21),strip.background=element_rect(fill="grey92",color="grey65",linewidth=.8),strip.text=element_text(size=18,face="plain"),legend.position="none",panel.border=element_rect(color="black",fill=NA,linewidth=.8),panel.grid=element_blank(),axis.line=element_line(color="black",linewidth=.8))
 ggsave(file.path(OUTDIR,paste0(out,".png")),p,width=11.5,height=6.8,dpi=300,bg="white"); ggsave(file.path(OUTDIR,paste0(out,".pdf")),p,width=11.5,height=6.8,bg="white")}
h=build(file.path(ROOT,"human","rds","human_auc_matrix.rds"),file.path(ROOT,"human","rds","human_cellInfo.rds"),"human"); plotp(h,"Human","SampleLevel_Human_EGR1_regulon_AUC_paired_connected_only")
m=build(file.path(ROOT,"mouse","rds","mouse_auc_matrix.rds"),file.path(ROOT,"mouse","rds","mouse_cellInfo.rds"),"mouse"); plotp(m,"Mouse","SampleLevel_Mouse_Egr1_regulon_AUC_paired_connected_only")
