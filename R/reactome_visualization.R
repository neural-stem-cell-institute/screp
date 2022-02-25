
#'reactome_visualization
#'
#'function to generate plot using the dataframe from the reactome_prep function
#' @param RP.df data frame produced by reactome_prep function
#' @param RP_adj same data frame as used in reactome_prep function
#' @param path_order a vector of parent pathway names from the rootName column in the RP.df. Can be used to
#'         plot the parent pathway names in desired order along the x-axis. Default is NULL which gives a random order.
#' @param RP data frame of Path IDS from Reactome database. IF left null uses the RP file incuded as data in the screp package.
#' @import 'foreach'
#' @import 'ggplot2'
#' @import 'reshape2'
#' @importFrom foreach %do%
#' @importFrom stats sd
#' @returns returns a list with the first item in the list being a data frame used in creating a
#'          plotting object and the second item being the ggplot2 object generated
#' @export
#' @examples
#' \dontrun{
#' react_vis2<-reactome_visualization(RP_ready,mRPR,path_order=z)
#' }

reactome_visualization<-function(RP.df,RP_adj,RP=NULL, path_order=NULL){
  b<-NULL
  
  if(is.null(RP)) {print("Function needs product of RPprep function")}

  x1<-table(RP.df$rootPath)
  x<-foreach(i=1:length(names(x1))) %do% {
    y<-RP_adj[which(RP_adj[,1]==names(x1[i])),2]
    if(length(y)==0) {y<-names(x1)[i]} else {y}
  }

  x<-foreach(b=1:length(x)) %do% {
    x2<-x[[b]]
    x3<-x2
    repeat{
      x3<-foreach(i=1:length(x3),.combine='c') %do% {
        y<-RP_adj[which(RP_adj[,1]==x3[i]),2]
      }
      if(length(x3)==0) {break} else {x2<-c(x2,x3)}
    }
    x2
  }

  RP_numbers<-foreach(i=1:length(x),.combine='c') %do% {length(x[[i]])}
  names(RP_numbers)<-RP[names(x1),]$V2



  path_mat<-foreach(i=0:max(RP.df$Cluster),.combine='cbind') %do% {
    y<-table(RP.df[RP.df$Cluster==i,]$rootName)
    z<-setdiff(names(RP_numbers),names(y))
    z1<-rep(0,length(z))
    names(z1)<-z
    y<-c(y,z1)
    y<-y[names(RP_numbers)]

  }
  colnames(path_mat)<-0:max(RP.df$Cluster)
  if(!(is.null(path_order))) {path_mat<-path_mat[path_order,]}

  RP_numbers<-RP_numbers[rownames(path_mat)]

  z1<-apply(path_mat,2,sum)
  z2<-apply(path_mat,2,function(x) x/sum(x)*100)
  z3<-apply(path_mat,2,function(x) (x-mean(x))/sd(x))

  y<-melt(path_mat)
  colnames(y)<-c("Path","Cluster","P.N")
  y1<-melt(z2)
  y$percent<-y1$value
  y1<-path_mat/RP_numbers*100
  y1<-melt(y1)
  y$Tree.per<-y1$value
  y<-y[-(which(y$`P.N`==0)),]
  y1<-y$Cluster
  for( i in 0:max(RP.df$Cluster)) {y1[which(y1==i)]<-rev(gg_color_hue(max(RP.df$Cluster)+1))[i+1]}
  y$colors<-y1
  p<-ggplot(y, aes(y=factor(Cluster), x=Path, size=percent, alpha=Tree.per,color=colors))
  p<-p + geom_point() + geom_text(aes(label=P.N,vjust=-1), alpha=1.0, size=4)
  p<-p + scale_alpha_continuous(range=c(0.3, 1)) + scale_size_area(max_size = 6)
  p<-p + theme_bw() + theme(axis.line = element_blank(),
                            axis.title = element_blank(),
                            panel.border = element_blank(),
                            panel.grid.major.y=element_line(color="gray20"),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            axis.text.x= element_text(angle=90,vjust=1,hjust=1),aspect.ratio=1/2 )

  p<-p + ggtitle("Reactome Pathways") + theme(plot.title = element_text(hjust = 0.5))
  p<-p + guides(color="none") + guides(alpha=guide_legend(override.aes=list(size=6,color=y$colors[1])))
  p<-p + theme(plot.title = element_text(size=24,margin = margin(10, 0, 10, 0)),
               legend.key=element_rect(fill='snow2'))
  p

  fin<-list(y,p)
  names(fin)<-c("Plotting_dataframe","plot_object")
  return(fin)
}
