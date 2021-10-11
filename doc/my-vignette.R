## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## -----------------------------------------------------------------------------
## Load libraries 
library(screp)
library(GO.db)

## Load data set for tutorial and get markers
data(pbmc)
markers<-pbmc@misc$markers

## -----------------------------------------------------------------------------

### First example providing only the marker table
enrich_results<-enrich_test(markers_df=markers)

### Second example generating and providing a clust_list
### Generate a list of gene vectors from the marker table

sig.marks<-markers[markers$p_val_adj<0.1,]

marks<-foreach(i=1:length(unique(sig.marks$cluster))) %do% {
  sig.marks[sig.marks$cluster==unique(sig.marks$cluster)[i],]$gene
}

names(marks)<-paste("cluster",unique(sig.marks$cluster),sep="_")

## check size of gene lists associated with clusters

foreach(i=1:length(marks),.combine='c') %do% {length(marks[[i]])}

## example of running the function with the cluster gene list
## enrich_results<-enrich_test(clust_list=marks)

## -----------------------------------------------------------------------------
## accessing the enrich_results object
head(enrich_results$Enriched_lists$Cluster_0$GOBP$data)[,2:5]

head(enrich_results$Enriched_lists$Cluster_3$REACTOME$data)[,2:5]

head(enrich_results$Enriched_df$Cluster_0$Pathways)[,2:5]

head(enrich_results$Enriched_df$Cluster_2$GOBP)[,2:5]

## ---- warning=FALSE-----------------------------------------------------------
## Make the vector for the GO terms variable
## Load Go.db package from Bioconductor

x <- Term(GOTERM)
goterms<-names(x)
names(goterms)<-x

## Run the GO_visualization function 
go_vis_res1<-GO_visualization(enrich_results$Enriched_df,markers_df=markers,goterms=goterms,numcats=10)

## ----fig.asp=1, fig.width=8, out.width="100%"---------------------------------
go_vis_res1$plot


## ---- warning=FALSE-----------------------------------------------------------
## Take the results of the GO_visualization function and choose the categories to plot.
newcats<-names(sort(table(go_vis_res1$GO_sem$parentTerm)))[1:10]
newcats

go_vis_res2<-GO_viz_choose(go_vis_res1,chosen_cats=newcats,markers_df=markers)

## ----fig.asp=1, fig.width=8, out.width="100%"---------------------------------
go_vis_res2

## -----------------------------------------------------------------------------
## read in and prep the "Complete List of Pathways" file after downloading
## RP<-read.delim("ReactomePathways.txt",as.is=T,header=F)
## for this tutorial we will use the "Complete List of Pathways" file loaded as part of the package
data(RP)
RP<-RPprep(RP,"Homo sapiens")

## check for NAs before preceeding
which(is.na(RP))

## read in the "Pathways hierarchy relationship" file after downloading
## RPR<-read.delim("ReactomePathwaysRelation.txt",as.is=T,header=F)
## for this tutorial we will use the "Pathways hierarchy relationship"  file loaded as part of the package
data(RPR)

## -----------------------------------------------------------------------------
## this function requires the RP and RPR objects generated previously. We will first modify the RPR object by removing the Immune System nodes from the adjacency table. To do so you must obtain the pathway ID from the Reactome database and remove rows where its in the first column
mRPR<-RPR[-(grep("R-HSA-168256",RPR[,1])),]
## Now we can use the reactome prep function
RP_ready<-reactome_prep(enrich_results$Enriched_df,RP=RP,RP_adj=mRPR)
head(RP_ready)

## -----------------------------------------------------------------------------
react_vis1<-reactome_visualization(RP_ready,mRPR)
head(react_vis1$Plotting_dataframe)

## ----fig.asp=0.75, fig.width=10, out.width="100%"-----------------------------
react_vis1$plot_object

## -----------------------------------------------------------------------------
## now we will generate the plot with the Pathways in alphabetical order
z<-sort(unique(RP_ready$rootName))
z
react_vis2<-reactome_visualization(RP_ready,mRPR,path_order=z)

## ----fig.asp=0.75, fig.width=10, out.width="100%"-----------------------------
react_vis2$plot_object

