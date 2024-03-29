#' enrich_test
#'
#' A function to streamline enrichment tests using of cluster marker genes in single cell seq after analysis with Seurat.
#' @param markers_df is a data frame generated by the FindAllMarkers function in the Seurat package.
#' @param p_val is the FDR cutoff to generate a list of cluster associated genes from the markers df.
#' @param cats_list is a list of annotation lists created by msigdb_gsets function including BIOCARTA, KEGG, REACTOME, GO:BP, GO:MF, and GO:CC in that order. Function generates list if not provided using species_x which is human by default
#' @param genome_genes should be the number of genes in your data set.
#' @param clust_list can be a user generated list of genes associated with clusters. Must be a list of vectors.
#' @param species_x species names, default is "Homo sapiens", refer to ?msigdbr::msigdbr for available species
#' @import 'Seurat'
#' @import 'hypeR'
#' @import 'foreach'
#' @importFrom foreach %do%
#' @return returns a list of length 2. The first object is a list of of lists of the object dervied from the enrichment tests,
#'           and the second is a data frame of the data used as input for the visulaization functions.
#' @export
#' @examples
#' \dontrun{
#' GO_and_Reactome_results<-enrich_test(seurat_markers)
#' }

enrich_test<-function(markers_df=NULL,p_val=0.1,cats_list=NULL,species_x="Homo sapiens",genome_genes=23459,clust_list=NULL) {
  i<-NULL
  if((is.null(markers_df) && is.null(clust_list))==TRUE) {stop("Function requires markers provided to either markers_df or clust_list")}
  if(is.null(clust_list)) {
    markers<-markers_df[markers_df$p_val_adj<p_val,]
    clust_list<-foreach(i=1:length(levels(markers$cluster)))%do% {markers[markers$cluster==levels(markers$cluster)[i],]$gene}
    names(clust_list)<-levels(markers$cluster)
    clust_list<-clust_list[lengths(clust_list)>0]
    clustlength<-length(clust_list)
  }
  if(!(is.null(markers_df))) {
    markers<-markers_df[markers_df$p_val_adj<p_val,]
    clustlength<-length(levels(markers$cluster))
  }
  if(!(is.null(clust_list))) {
    clustlength<-length(clust_list)
  }

  if(is.null(cats_list)){
    BIOCARTA <- msigdb_gsets(species=species_x, category="C2", subcategory="CP:BIOCARTA")
    KEGG     <- msigdb_gsets(species=species_x, category="C2", subcategory="CP:KEGG")
    REACTOME <- msigdb_gsets(species=species_x, category="C2", subcategory="CP:REACTOME")
    GOBP <- msigdb_gsets(species=species_x, category="C5", subcategory="GO:BP")
    GOMF <- msigdb_gsets(species=species_x, category="C5", subcategory="MF")
    GOCC <- msigdb_gsets(species=species_x, category="C5", subcategory="CC")}
  else {
    BIOCARTA <- cats_list$BIOCARTA
    KEGG     <- cats_list$KEGG
    REACTOME <- cats_list$REACTOME
    GOBP <- cats_list$GOBP
    GOMF <- cats_list$GOMF
    GOCC <- cats_list$GOCC
  }

  cluster_enriched_list<-foreach(i=1:length(clust_list)) %do% {
    Path_BIOCARTA<-hypeR(clust_list[[i]], BIOCARTA, background=genome_genes, fdr=0.01)
    Path_Reactome<-hypeR(clust_list[[i]], REACTOME, background=genome_genes, fdr=0.01)
    Path_KEGG<-hypeR(clust_list[[i]], KEGG, background=genome_genes, fdr=0.01)
    GOBioProc<-hypeR(clust_list[[i]], GOBP, background=genome_genes, fdr=0.01)
    GOMoleFunc<-hypeR(clust_list[[i]], GOMF, background=genome_genes, fdr=0.01)
    GOCellComp<-hypeR(clust_list[[i]], GOCC, background=genome_genes, fdr=0.01)
    x<-list(Path_BIOCARTA,Path_Reactome,Path_KEGG,GOBioProc,GOMoleFunc,GOCellComp)
    names(x)<-c("BIOCARTA","REACTOME","KEGG","GOBP","GOMF","GOCC")
    x
  }
  names(cluster_enriched_list)<-names(clust_list)

  cluster_enriched_df<-foreach(i=1:length(clust_list)) %do% {
    Path1<-cluster_enriched_list[[i]][[1]]$data
    Path2<-cluster_enriched_list[[i]][[2]]$data
    Path3<-cluster_enriched_list[[i]][[3]]$data
    Pathways<-rbind(Path1,Path2,Path3)
    BP<-cluster_enriched_list[[i]][[4]]$data
    MF<-cluster_enriched_list[[i]][[5]]$data
    CC<-cluster_enriched_list[[i]][[6]]$data
    x<-list(Pathways,BP,MF,CC)
    names(x)<-c("Pathways","GOBP","GOMF","GOCC")
    x
  }
  
  names(cluster_enriched_df)<-names(clust_list)
  fin<-list(cluster_enriched_list,cluster_enriched_df)
  names(fin)<-c("Enriched_lists","Enriched_df")
  return(fin)
}

