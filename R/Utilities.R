
#'gg_color_hue
#'
#' internal function that acts as a ggplot color replicator helper function
#' @param n number of colors needed
#' @return returns a vector of colors for use in plotting
#' @export
#' @examples gg_color_hue(8)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' RPprep
#'
#' function to prepare the Reactome Pathway file downloaded from the ReactomeDB for downstream use.
#' @param RP file downloaded from the reactome database or the provided annotation file. If using a
#'   a file from the reactome database some NAs will likely result and need to be addressed before continuing.  A corrected RP file for Humans is included in the package.
#' @param species_x genus of target species. Default is "Homo" (human)
#' @returns returns a data frame suitable for use in reactome_prep function
#' @export
#' @examples RP<-RPprep(RP,"Homo sapiens")


RPprep<-function(RP,species_x="Homo"){
  RP<-RP[grep(species_x,RP[,3]),]
  rownames(RP)<-RP$V1
  RP$caps<-toupper(RP$V2)
  RP$caps<-gsub("\\(","",RP$caps)
  RP$caps<-gsub("\\)","",RP$caps)
  RP$caps<-gsub(" - "," ",RP$caps)
  RP$caps<-gsub("-"," ",RP$caps)
  RP$caps<-gsub(", "," ",RP$caps)
  RP$caps<-gsub(","," ",RP$caps)
  RP$caps<-gsub("\\+","",RP$caps)
  RP$caps<-gsub(": "," ",RP$caps)
  RP$caps<-gsub(":"," ",RP$caps)
  RP$caps<-gsub("; "," ",RP$caps)
  RP$caps<-gsub(";"," ",RP$caps)
  RP$caps<-gsub("\\/"," ",RP$caps)
  RP$caps<-gsub("\\.","",RP$caps)
  RP$caps<-gsub("\\'","",RP$caps)
  RP$caps<-gsub("& ","",RP$caps)
  RP$caps<-gsub("  "," ",RP$caps)
  RP$caps<-trimws(RP$caps)
  return(RP)
}

#' reactome_prep
#'
#' function to make a data to use with visualization function
#' @param cluster_enriched_df data frame produced by enrich_test function
#' @param RP product of RPprep function
#' @param RP_adj Pathway relationship file from Reactome. To decrease time required should limit to organism of interest
#'            by selecting for organism, e.g. RPH<-RP[grep("HSA",RPR[,1]),]. Can be modified by user
#'            to use some 2nd level terms (e.g. RP_adj<-RPH[-(grep("R-HSA-8953897",RPH[,1])),]  ).
#' @import 'foreach'
#' @importFrom foreach %do%
#' @returns returns a data frame for use in the reactome_visualization function
#' @export
#' @examples RP_ready<-reactome_prep(enrich_results$Enriched_df,RP=RP,RP_adj=mRPR)

reactome_prep<-function(cluster_enriched_df,RP,RP_adj){


  clx<-foreach(i=1:length(cluster_enriched_df),.combine='rbind') %do% {
    if(dim(cluster_enriched_df[[i]]$Pathways)==0} {} else{
    x<-cluster_enriched_df[[i]]$Pathways
    x$Cluster<-i-1
    x[grep("REACTOME",x$label),]
      }
  }
  x<-strsplit(clx$label,"REACTOME_")
  x<-unlist(x)[seq(2,length(unlist(x)),2)]
  x<-gsub("_", " ", x)
  clx$label1<-x
  clx$PathID<-foreach(i=1:length(clx$label1),.combine="c") %dopar% {
    x<-RP[which(RP$caps==clx$label1[i]),1]
    if(length(x)==0) {x<-NA} else{x}
  }

  x<-foreach(i=1:length(clx$PathID),.combine='c') %do% {
    y<-RP_adj[which(RP_adj[,2]==clx$PathID[i]),1]
    if(length(y)==0) {y<-clx$PathID[i]} else {y[1]}
  }

  repeat{
    z1<-length(unique(x))
    x<-foreach(i=1:length(x),.combine='c') %do% {
      y<-RP_adj[which(RP_adj[,2]==x[i]),1]
      if(length(y)==0) {y<-x[i]} else {y[1]}
    }
    if(length(unique(x))==z1) {break}
  }

  clx$rootPath<-x
  clx$rootName<-RP[x,]$V2
  return(clx)
}

#' getGoTerm
#' Get the description of a GO term, function is from the rrvgo package
#'
#' @param x GO terms
#' @return the Term slot in GO.db::GOTERM[[x]]
getGoTerm <- function(x) {
  sapply(x, function(x) tryCatch(GO.db::GOTERM[[x]]@Term, error=function(e) NA))
}

#' prepGOTERMS
#'
#'
#' function to generate a GO vector from the GO.db
#' @param x is the vector from Term(GOTERM) using the GO.db package
#' @import 'GO.db'
#' @import 'LSAfun'
#' @import 'stringr'
#' @export
 prepGOTERMS<-function(x){
  goterms<-names(x)
  x<-breakdown(x)
  x<-str_squish(x)
  names(goterms)<-x
  return(goterms)
 }




