#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(patchwork)
  library(showtext)
  library(Seurat)
})

source(file.path("R", "load_config.R"))


############################################################
# 0. Font
############################################################

font_add(
  "Arial",
  A4OL_FONT_FILE
)

showtext_auto()


############################################################
# 1. Path
############################################################


base_dir <- file.path(A4OL_FIGURE_ROOT, "olig-c")


deg_dir <-
file.path(
  base_dir,
  "ALL4_OLIG_SUBCLUSTERS_DEG_GO_KEGG_TOP50_UMAP_20260802"
)


out_dir <-
file.path(
  base_dir,
  "OLIGO4_FINAL_UPDATE"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



############################################################
# 2. Marker dotplot
############################################################


seu_file <-
file.path(
  base_dir,
  "seurat_oligodendrocytes_subset_annotated_A4_OLs.rds"
)

obj <- readRDS(seu_file)


Idents(obj) <- "oligo_annotated_cluster"



marker_list <- list(

"A4-OLs" =
c(
"Uqcr11",
"Ndufc1",
"Rpl39"
),

"Oligo1" =
c(
"Enpp6",
"Cyp27a1",
"Cldn14"
),

"Oligo2" =
c(
"Dpf3",
"Zfp521",
"Grin2c"
),

"Oligo3" =
c(
"Myo16",
"Kcnmb2",
"Kirrel"
)

)



dot_genes <- unlist(marker_list)


dot_data <- DotPlot(
  obj,
  features = dot_genes,
  group.by = "oligo_annotated_cluster"
)$data



# 调整顺序

dot_data$id <- factor(
  dot_data$id,
  levels=c(
    "A4-OLs",
    "Oligo1",
    "Oligo2",
    "Oligo3"
  )
)


dot_data$features.plot <-
factor(
dot_data$features.plot,
levels = dot_genes
)



p_dot <- ggplot(
dot_data,
aes(
x=features.plot,
y=id
)
)+

geom_point(
aes(
size=pct.exp,
color=avg.exp.scaled
)
)+

scale_color_gradient(
low="#FFFFFF",
high="#D73027"
)+

scale_size_continuous(
range=c(2,12)
)+

theme_classic(base_family="Arial")+

theme(

axis.text.x =
element_text(
angle=45,
hjust=1,
size=28
),

axis.text.y =
element_text(
size=30
),

axis.title =
element_blank(),

legend.title =
element_text(
size=25
),

legend.text =
element_text(
size=22
),

plot.title =
element_text(
size=32,
face="bold"
)

)+

labs(
title=
"A4-OLs and Oligo1-3 marker expression"
)



ggsave(
file.path(
out_dir,
"OLIGO4_marker_dotplot.png"
),
p_dot,
width=12,
height=8,
dpi=300
)





############################################################
# 3. Enrichment plot function
############################################################


draw_enrichment <- function(
group_name,
up_terms,
down_terms
){


up_go <-
read.csv(
file.path(
deg_dir,
paste0(
"up_",
group_name,
"_vs_Other_GO_BP_Upregulated.csv"
)
)
)


up_kegg <-
read.csv(
file.path(
deg_dir,
paste0(
"up_",
group_name,
"_vs_Other_KEGG_Upregulated.csv"
)
)
)


down_go <-
read.csv(
file.path(
deg_dir,
paste0(
"down_",
group_name,
"_vs_Other_GO_BP_Downregulated.csv"
)
)
)


down_kegg <-
read.csv(
file.path(
deg_dir,
paste0(
"down_",
group_name,
"_vs_Other_KEGG_Downregulated.csv"
)
)
)



pick <- function(df,terms,db,direction){

df %>%
filter(
Description %in% terms
)%>%
mutate(
Database=db,
Direction=direction
)

}



df <- bind_rows(

pick(
up_go,
up_terms,
"GO",
"Up"
),

pick(
up_kegg,
up_terms,
"KEGG",
"Up"
),

pick(
down_go,
down_terms,
"GO",
"Down"
),

pick(
down_kegg,
down_terms,
"KEGG",
"Down"
)

)



df <-
df %>%
mutate(

log10p=-log10(p.adjust),

x=ifelse(
Direction=="Up",
log10p,
-log10p
)

)



df$Description <-
factor(
df$Description,
levels=
rev(
unique(
df$Description[
order(df$x)
]
)
)
)



p <- ggplot(
df,
aes(
x=x,
y=Description
)
)+


geom_vline(
xintercept=0,
linetype="dashed",
color="grey60"
)+


geom_bar(
aes(fill=Direction),
stat="identity",
width=.65
)+


geom_point(
aes(
size=Count,
color=p.adjust,
shape=Database
)
)+


scale_fill_manual(
values=c(
"Up"="#E64B35",
"Down"="#4DBBD5"
)
)+


scale_shape_manual(
values=c(
GO=16,
KEGG=18
)
)+


scale_color_gradient(
low="#2E8B57",
high="#FFD700"
)+


scale_size_continuous(
range=c(3,10)
)+


labs(
x=expression(-log[10](padj)),
y=NULL
)+


theme_minimal(base_family="Arial")+

theme(

axis.text.y=
element_text(
size=28
),

axis.text.x=
element_text(
size=25
),

legend.text=
element_text(
size=22
),

legend.title=
element_text(
size=25
),

panel.grid.major.y=
element_blank()

)



ggsave(
file.path(
out_dir,
paste0(
group_name,
"_selected_enrichment.png"
)
),
p,
width=14,
height=10,
dpi=300
)


}



############################################################
# 4. pathways
############################################################


pathways <- list(


A4.OLs=list(

up=c(
"Oxidative phosphorylation",
"Aerobic respiration",
"Aerobic electron transport chain",
"ATP synthesis coupled electron transport",
"Mitochondrial ATP synthesis coupled electron transport",
"Cytoplasmic translation",
"Ribosome"
),

down=c(
"Chromatin remodeling",
"DNA repair",
"Heterochromatin formation",
"Cell-cell junction"
)

),



Oligo1=list(

up=c(
"Myelination",
"Myelin assembly",
"Ensheathment of neurons",
"Actin cytoskeleton organization",
"Cell-cell junction organization"
),

down=c(
"Oxidative phosphorylation",
"Ribosome",
"Regulation of trans-synaptic signaling"
)

),



Oligo2=list(

up=c(
"Chromatin remodeling",
"Heterochromatin formation",
"DNA replication",
"Histone demethylase activity",
"Calcium signaling pathway"
),

down=c(
"Ribosome",
"Oxidative phosphorylation",
"Myelin sheath"
)

),



Oligo3=list(

up=c(
"Synapse assembly",
"Regulation of membrane potential",
"GABAergic synapse",
"Calcium signaling pathway",
"Neuroactive ligand-receptor interaction"
),

down=c(
"Ribosome",
"Oxidative phosphorylation",
"Myelin sheath"
)

)

)




############################################################
# 5. Run
############################################################


for(g in names(pathways)){

draw_enrichment(
g,
pathways[[g]]$up,
pathways[[g]]$down
)

}


cat(
"ALL DONE\nOutput:",
out_dir,
"\n"
)
