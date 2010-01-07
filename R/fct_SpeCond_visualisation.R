plotNormalMixture <- function(probeset_name,prefix.file,expression_values_p,param_p,class_col,M_s_p,null_p){
  
  plot_names=c("exp_profile_","hist_")

  p_names=sapply(1:length(plot_names), function(i) paste(prefix.file,"_",plot_names[i],probeset_name,c(".png",".pdf"),sep=""))
  link_html=vector()
  
  png(p_names[1,1])
  plot(expression_values_p,ylab="expression value",xlab= "condition",col=class_col,pch=20)
  points(x=which(M_s_p==-1),y=expression_values_p[M_s_p==-1],col="red",pch=6,cex=2)
  points(x=which(M_s_p==1),y=expression_values_p[M_s_p==1],col="red",pch=2,cex=2)

  link_html[1]=hwriteImage(p_names[1,1])
  dev.off()
  
  G=ncol(param_p)

  png(p_names[1,2])
  x <- seq( min(expression_values_p-1), max(expression_values_p+1), length = 500 )
  hist( expression_values_p, breaks = 15, prob = TRUE, main="",xlab="expression value" )
  points(expression_values_p, rep(0, length(expression_values_p)), pch = 1)

  f1 <- function(x)  param_p[1,1]*dnorm( x,param_p[2,1],param_p[3,1])
  f2 <- function(x)  param_p[1,2]*dnorm( x,param_p[2,2],param_p[3,2])
  f3 <- function(x)  param_p[1,3]*dnorm( x,param_p[2,3],param_p[3,3])
  f4 <- function(x)  param_p[1,4]*dnorm( x,param_p[2,4],param_p[3,4])
  f5 <- function(x)  param_p[1,5]*dnorm( x,param_p[2,5],param_p[3,5])
  
  if(G==2){
    lines( x, f1(x), col = "blue" )
    lines( x, f2(x), col = "green" )
    legend("topright",legend=c("Normal 1","Normal 2","null distribution"),col=c("blue","green","red"),lty=1)
  }
  
  if(G==3){
    lines( x, f1(x), col = "blue" )
    lines( x, f2(x), col = "green" )
    lines( x, f3(x), col = "orange" )
  legend("topright",legend=c("Normal 1","Normal 2","Normal 3","null distribution"),col=c("blue","green","orange","red"),lty=1)
  }

   if(G==4){
     lines( x, f1(x), col = "blue" )
     lines( x, f2(x), col = "green" )
     lines( x, f3(x), col = "orange" )
     lines( x, f4(x), col = "purple" )
     legend("topright",legend=c("Normal 1","Normal 2","Normal 3","Normal 4"),col=c("blue","green","orange","purple"),lty=1)
   }
  
  if(G==5){    
     lines( x, f1(x), col = "blue" )
     lines( x, f2(x), col = "green" )
     lines( x, f3(x), col = "orange" )
     lines( x, f4(x), col = "purple" )
     lines( x, f5(x), col = "cyan" )
     
     legend("topright",legend=c("Normal 1","Normal 2","Normal 3","Normal 4", "Normal 5"),col=c("blue","green","orange","purple", "cyan"),lty=1)
   }
  
  ##  Combination of the normal distributions
  n=which(null_p==1)
  W_sum=sum(param_p[1,n])
  W=param_p[1,n]/W_sum

  f_null<- function(x)  apply((sapply(1:length(n), function(i) W[i]*dnorm(x,param_p[2,n[i]],param_p[3,n[i]]))),1,sum)
  lines( x, f_null(x), col = "red")
  link_html[2]=hwriteImage(p_names[1,2])
  dev.off()  
  return(link_html)
}


getExpressionpatternLegend<- function(G){

  if(G==1){
    cases=matrix('&nbsp;', nr=2, nc=1)
    bgcolor=c('#000DFE','#FE0202') ##bleu,red
    patchwork=hwrite(cases, bgcolor=bgcolor,col.width=rep('18px',1),table.style='text-align:center; border-collapse:collapse')
    legend=c("Normal distribution 1","Specific expresison")
    h_l=hwrite(legend,dim=c(2,1),border=0)
    h_legend=hwrite(c(patchwork,h_l),border=0)
  }

  if(G==2){
    cases=matrix('&nbsp;', nr=2, nc=1)
    bgcolor=c('#000DFE','#04C610') ##bleu,green
    patchwork=hwrite(cases, bgcolor=bgcolor,col.width=rep('18px',1),table.style='text-align:center; border-collapse:collapse')
    legend=c("Normal distribution 1","Normal distribution 2")
    h_l=hwrite(legend,dim=c(2,1),border=0)
    cases2=matrix('&nbsp;', nr=2, nc=1)
    bgcolor2=c('#FE0202',"") ## red
    patchwork2=hwrite(cases2, bgcolor=bgcolor2,col.width=rep('18px',1),table.style='text-align:center; border-collapse:collapse')
    legend2=c("Specific expression","-")
    h_l2=hwrite(legend2,dim=c(2,1),border=0)
    h_legend=hwrite(c(patchwork,h_l,patchwork2,h_l2),border=0)
  }

  if(G==3){
    cases=matrix('&nbsp;', nr=2, nc=1)
    bgcolor=c('#000DFE','#04C610') ##bleu,green
    patchwork=hwrite(cases, bgcolor=bgcolor,col.width=rep('18px',1),table.style='text-align:center; border-collapse:collapse')
    legend=c("Normal distribution 1","Normal distribution 2")
    h_l=hwrite(legend,dim=c(2,1),border=0)
    cases2=matrix('&nbsp;', nr=2, nc=1)
    bgcolor2=c('#FFA413','#FE0202') ##orange,red
    patchwork2=hwrite(cases2, bgcolor=bgcolor2,col.width=rep('18px',1),table.style='text-align:center; border-collapse:collapse')
    legend2=c("Normal distribution 3","Specific expression")
    h_l2=hwrite(legend2,dim=c(2,1),border=0)
    h_legend=hwrite(c(patchwork,h_l,patchwork2,h_l2),border=0)
  }

  if(G>3){
    h_legend=hwrite("Legend! (G>3 !!!!!!!)")
  }
  return(h_legend)
}

createSingleGeneHtmlPage <- function(index.html.link,prefix.file,p_id,gene.name,expression_values_p,param.detection=NULL,fit_p,selective_p,param_p,M_class_col_p,M_s_p,L_pv_s_p,null_p,min_lk_values_p,rsd_p,n_link){

  p=openPage(paste(prefix.file,"_",gene.name,".html",sep=""))
  
  if(!is.null(index.html.link)){
    
    if(p_id==1){
      hwrite(c('Previous','Up','Next'),p,link=c(index.html.link,index.html.link,n_link[p_id+1]),border=0)
  }
    else{
      if(p_id==length(n_link)){
        hwrite(c('Previous','Up','Next'),p,link=c(n_link[p_id-1],index.html.link,index.html.link),border=0)
    }
      else{
        hwrite(c('Previous','Up','Next'),p,link=c(n_link[p_id-1],index.html.link,n_link[p_id+1]),border=0)
      }
    }
  }
    
  hwrite(c(gene.name,selective_p),p,heading=2,dim=c(1,2),border=0,style=matrix(c('color:#000DFE',NA),nr=1,nc=2),cellspacing=10)
  
  if(!is.null(param.detection)){
      hwrite('Parameter set used for the detection of the specific condition(s):',p,heading=2)
      hwrite(param.detection,p,row.style=list('text-align:center, font-weight:bold'),table.style='text-align:center; border-collapse:collapse',col.width=c('75px','75px','75px','75px','75px','75px','75px','75px'))
    }

  nb_G_initial=fit_p$G_initial
  nb_G=fit_p$G
  t1=hwrite("Number of distribution detected= ")
  hwrite(c(t1,nb_G_initial,nb_G),p, border=0)

  link_html=plotNormalMixture(gene.name,prefix.file,expression_values_p,param_p,M_class_col_p,M_s_p,null_p)
  h_legend_expression=hwrite(paste("Figure 1: Expression values of ", gene.name," across conditions",sep=""),style='font-size:10pt')
  h_legend_hist=hwrite("Figure 2: Density of the expression values, normal(s) and the null distribution",style='font-size:10pt')
  
  param=rbind(round(param_p,5),as.character(null_p))
  rownames(param)=c(rownames(param_p),"null distribution")
  rowname_color=c('#637FF7','#00E50B','#FFAC33','#A533FF','#FFFFFF')
  
  h_L_param=hwrite(param,table.style='text-align:center; border-collapse:collapse',row.style=list('text-align:center'),col.width=c('120px',rep('150px',ncol(param))),row.bgcolor=list(c(NA,rowname_color[1:nb_G]),'null distribution'=c(NA,rep('#EF1D1D',nb_G))))
  h_legend=getExpressionpatternLegend(nb_G)
  
  if(!is.null(min_lk_values_p) && !is.null(rsd_p)){
      h_min_lk_values=hwrite(c("min loglikelihood",min_lk_values_p),dim=c(1,(1+length(min_lk_values_p))),table.style='text-align:center; border-collapse:collapse',col.width=c('120px',rep('150px',length(min_lk_values_p))))
      h_rsd=hwrite(c("sd ratio",rsd_p),dim=c(1,(1+length(rsd_p))),table.style='text-align:center; border-collapse:collapse',col.width=c('120px',rep('150px',length(rsd_p))))
      hwrite(c(link_html[1],h_legend_expression,"",h_legend,"","",link_html[2],h_legend_hist,"",h_L_param,h_min_lk_values,h_rsd),dim=c(6,2),p,border=0)
    }
  else{
    if(!is.null(rsd_p)){
      h_rsd=hwrite(c("sd ratio",rsd_p),dim=c(1,(1+length(rsd_p))),table.style='text-align:center; border-collapse:collapse',col.width=c('120px',rep('150px',length(rsd_p))))
      hwrite(c(link_html[1],h_legend_expression,"",h_legend,"",link_html[2],h_legend_hist,"",h_L_param,h_rsd),dim=c(5,2),p,border=0)
    }
    else{
      hwrite(c(link_html[1],h_legend_expression,"",h_legend,link_html[2],h_legend_hist,"",h_L_param),dim=c(4,2),p,border=0)
    }
  }
  nb_normal=sum(null_p)
  hwrite(paste("The number of normal distribution used in the null model is: ",nb_normal,sep=""),p)
  
  if(!length(L_pv_s_p)==0){
    hwrite("Significant p_values",p,heading=2)
    result_pv_s=hwrite(as.matrix(round(L_pv_s_p,5)),table.style='text-align:center; border-collapse:collapse',col.width=c('200px','150px'))
    hwrite(c(result_pv_s,""),p,border=0)
  }
  
  hwrite("",p,br=TRUE)
  hwrite("",p,br=TRUE)
  rversion = sessionInfo()$R.version$version.string
  hwrite(paste("Page created by the ", hwrite("SpeCond",link="http://www.bioconductor.org/packages/release/bioc/")," package using hwriter under ",rversion,collapse=""), p, style='font-size:8pt')
  closePage(p)
}


getProfileHeatmap <- function(M_profile,file_name,colors=NULL){

  file_n=paste(file_name,c("png","pdf"),sep=".")
  print(file_n)
  png(file_n[1],width=800,height=800)
  heatmap(M_profile,scale='none',margin=c(8,8),col=colors)
  dev.off()
  pdf(file_n[2])
  heatmap(M_profile,scale='none',margin=c(8,8),col=colors)
  dev.off()
}

getHeatmap <- function(M_values,file_name,probe_set=NULL,colors=NULL,colRowSide=NULL,colrange=NULL){
  
  file_n=paste(file_name,c("png","pdf"),sep=".")
  print(file_name)
  png(file_n[1],width=800,height=800)
  if(is.null(probe_set)){
    if(is.null(colRowSide)){
      heatmap(M_values,scale='none',margin=c(8,8),col=colors)
    }
    else{
      heatmap(M_values,scale='none',margin=c(8,8),col=colors,RowSideColors=colRowSide)
      colorbar.plot(0,1,c(0:length(colrange)),col=colrange)
    }
  }
  else{
    heatmap(M_values[rownames(M_values) %in% probe_set,],scale='none',margin=c(8,8),col=colors)
  }
  dev.off()

  pdf(file_n[2])
  if(is.null(probe_set)){
    if(is.null(colRowSide)){
      heatmap(M_values,scale='none',margin=c(8,8),col=colors)
    }
    else{
      heatmap(M_values,scale='none',margin=c(8,8),col=colors,RowSideColors=colRowSide)
      colorbar.plot(0,1,c(0:length(colrange)),col=colrange)
    }
  }
  else{
    heatmap(M_values[rownames(M_values) %in% probe_set,],scale='none',margin=c(8,8),col=colors)
  }
  dev.off()

}
## return(M_s,M_s_sum_row,M_s_sum_column,L_pv_s,selective)
## getListSelectiveResult()
## return(M_selective,M_selective_sum_row,M_selective_sum_column)
## getProfile()
## return(M_selective_profile,M_selective_profile_unique,M_selective_profile_table)

writeSpeCondResult <- function(L.specific.result, file.name.profile="specific_profile.txt", file.specific.gene="list_specific_gene.txt",file.name.unique.profile="specific_unique_profile.txt"){

  print(paste("write file:",file.name.profile))
  print(paste("write file:",file.specific.gene))
  write.table(L.specific.result$M.specific,file=file.name.profile,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
  write.table(rownames(L.specific.result$M.specific),file=file.specific.gene,sep=",",row.names=FALSE,col.names=FALSE,quote=FALSE)

  L.specific.result.profile=getProfile(L.specific.result$M.specific)
  ##   the rows in M.specific.profile.unique and M.specific.profile.table should correspond
  ##   return(M.specific.profile,M.specific.profile.unique,M.specific.profile.table)
  print(paste("write file:",file.name.unique.profile))
  M.unique.profile=cbind(L.specific.result.profile$M.specific.profile.unique,L.specific.result.profile$M.specific.profile.table)
  write.table(M.unique.profile,file=file.name.unique.profile,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
}
## if full.list.gene=TRUE, write at the end of the row the probe set names which have the profile described on the row
writeUniqueProfileSpecificResult<- function(L.specific.result, file.name.unique.profile="specific.unique_profile.txt", full.list.gene=FALSE){
  
  L.specific.result.profile=getProfile(L.specific.result$M.specific)
  print(paste("write file:",file.name.unique.profile))
  if(full.list.gene){
    List.gene=lapply(1:nrow(L.specific.result.profile[[3]]),function(i) which(L.specific.result.profile[[1]][,1]==L.specific.result.profile[[3]][i,1]))
    names(List.gene)=L.specific.result.profile[[3]][,1]
    gene.names=sapply(1:length(List.gene),function(i) paste(names(List.gene[[i]]),collapse=","))
    M.unique.profile=cbind(L.specific.result.profile[[2]],L.specific.result.profile[[3]],gene.names)
  }
  else{
    M.unique.profile=cbind(L.specific.result.profile[[2]],L.specific.result.profile[[3]])
  }
  write.table(M.unique.profile,file=file.name.unique.profile,row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

}

## ## Write a file containing the list of probe_sets, if they have been detected as condition-specific or not (S/N), for how many conditions in total, how many condition as up-regulated, how many condition as down-regulated, in which condition for up-regulated and down-regulated.

writeGeneResult <- function(L.specific.result,file.name.result.gene="gene_summary_result.txt",gene.names=NULL){
  
  M.result.gene=getSelectiveResultTable(L.specific.result$M.specific.all,L.specific.result$M.specific.sum.row,gene.names)
  print(paste("write file:",file.name.result.gene))
  write.table(M.result.gene,file=file.name.result.gene,row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
}

## Create the full experiment result html file
getFullHtmlSpeCondResult<- function(SpeCondResult=NULL, L.specific.result=NULL, param.detection=NULL, page.name="SpeCond_result", page.title="Condition-specific analysis results", prefix.file=NULL, outdir="General_Result", sort.condition="all", gene.page.info=NULL, heatmap.profile=TRUE, heatmap.expression=FALSE, heatmap.unique.profile=FALSE, expressionMatrix=NULL){

  nb_table=1
  nb_figure=1

    if(is.null(prefix.file)){
      if(!is.null(SpeCondResult)){
        prefix.file=SpeCondResult$prefix.file
      }
      else{
        stop("Error you need to enter a prefix.file value or to implement the SpeCondResult attribute (sp_list class object from the SpeCond function)")
      }
    }
  
  olddir=getwd()
  on.exit(setwd(olddir))
  outdir=paste(prefix.file,outdir,sep="_")
  p=openPage(paste(prefix.file,"_",page.name,".html",sep=""))
  message(sprintf("The full result html page '%s' will be created in the current directory '%s'",paste(prefix.file,"_",page.name,".html",sep=""),getwd()))

  createOutdir(outdir)

  if(!is.null(SpeCondResult) && !is(SpeCondResult,"sp_list")){
    stop("Error, the SpeCondResult attribute must be of type sp_list (i.e output of the SpeCond function)")
    
  }
  if(is.null(L.specific.result)) {
    if(!is.null(SpeCondResult) && is.null(L.specific.result)){
      L.specific.result=SpeCondResult$specificResult$L.specific.result
    }
    else{
      print("Error, at least the SpeCondResult or the L.specific.result must not be NULL")
      return()
    }
  }

  nb_condition=ncol(L.specific.result$M.specific)
  hwrite(page.title,p,heading=1,center=TRUE)

  if(is.null(param.detection)){
    if(!is.null(SpeCondResult)){
      param.detection=SpeCondResult$param.detection
    }
  }
  if(!is.null(param.detection)){
    hwrite('Parameter set used to detect the specific condition(s)',p,heading=2)
    hwrite(param.detection,p,table.style='text-align:center; border-collapse:collapse',col.width=c('75px','75px','75px','75px','75px','75px','75px','75px'),row.bgcolor='#9af')
    hwrite("",p,br=TRUE)
    hwrite(paste("Table ",nb_table,": The parameter set used for the analysis",sep=""),p,br=TRUE)
    nb_table=nb_table+1
  }
  else{
    print("Advice: Implement the param.detection attribute or used the SpeCondResult attribute")
    print(paste("to present the parameter set used on the ",prefix.file,"_",page.name,".html"," page.",sep=""))
  }
  
  if(!is.null(L.specific.result$M.specific)){
    M.specific.sum.row=apply(abs(L.specific.result$M.specific),1,sum)

    tab=t(as.matrix(table(M.specific.sum.row)))
    specific_table_nb_tissues=matrix(c(colnames(tab),as.vector(tab)),nrow=2,ncol=ncol(tab),byrow=TRUE)
    rownames(specific_table_nb_tissues)=c("# conditions","# genes")
    specific_table_nb_tissues=cbind(specific_table_nb_tissues,c("-","-"),c("sum",sum(tab)))
    
    hwrite("Specific genes - conditions",p,heading=2)
    hwrite(specific_table_nb_tissues,p,table.style='text-align:center; border-collapse:collapse', row.bgcolor=list("# conditions"='#9af',NA),col.width=c('100px',rep('50px',ncol(specific_table_nb_tissues))))
    hwrite("",p,br=TRUE)
    hwrite(paste("Table ",nb_table,": Number of genes detected as specific in the different number of conditions",sep=""),p,br=TRUE)
    nb_table=nb_table+1

    print("The following files are created in the directory:")
    print(getwd())
    barplot_name=paste(prefix.file,"_barplot_nb_tissue_nb_genes.png",sep="")
    png(barplot_name)
    print(barplot_name)
    barplot(tab,col="purple",xlab="# conditions",ylab="# genes")
    dev.off()

##   plot the number of condition-specific genes by conditions
    nb_specific_in_source_name=paste(prefix.file,"_nb_specific_gene_in_condition.png",sep="")
    png(nb_specific_in_source_name)
    print(nb_specific_in_source_name)
    par(mar=c(14.1,4.1,4.1,2.1),mgp=c(3,1,0))
    
    if(sort.condition=="all"){
      M_ord=order(L.specific.result$M.specific.sum.column[3,])
    }
    else{
      if(sort.condition=="positive"){
        M_ord=order(L.specific.result$M.specific.sum.column[1,])
      }
      else{
        if(sort.condition=="negative"){
          M_ord=order(L.specific.result$M.specific.sum.column[2,])
        }
        else{
          M_ord=c(1:ncol(L.specific.result$M_specific_sum_colunm))
        }
      }
    }
    plot(L.specific.result$M.specific.sum.column[3,M_ord],ylim=c(min(0,min(L.specific.result$M.specific.sum.column)-10),max(L.specific.result$M.specific.sum.column)+10),pch=8,axes=FALSE,xlab="",ylab="# specific genes")
    axis(1,1:ncol(L.specific.result$M.specific.sum.column),colnames(L.specific.result$M.specific.sum.column)[M_ord],las=3, cex=0.3)
    axis(2,seq(0,max(L.specific.result$M.specific.sum.column),by=200))
    points(L.specific.result$M.specific.sum.column[1,M_ord],pch=3,col="blue")
    points(L.specific.result$M.specific.sum.column[2,M_ord],pch=45,col="green")
    legend("topleft",legend=c("# specific genes ","# up specific genes","# down specific genes"),col=c("black","blue","green"),pch=c(8,3,45))
    mtext("condition",1,10)
    dev.off()
    
    hwriteImage(paste(outdir,"/",barplot_name,sep=""),p,br=TRUE)
    hwrite(paste("Figure ",nb_figure,": Barplot representing the number of specific genes (y axis) in function of the number of their specific condition (x axis)",sep=""),p,br=TRUE)
    nb_figure=nb_figure+1
    hwriteImage(paste(outdir,"/",nb_specific_in_source_name,sep=""),p,br=TRUE)
    hwrite(paste("Figure ",nb_figure,": The number of specific genes detected in each of the conditions",sep=""),p,br=TRUE)
    nb_figure=nb_figure+1

    ##  Heapmap visualisation:
    ##   library('RColorBrewer')
    colHeatmap=colorRampPalette(brewer.pal(10,"RdBu"))(256)
    hwrite("Heatmap representation",p,heading=2)
    probe_set_names=rownames(L.specific.result$M.specific)
    if(heatmap.expression){
      hwrite("Expression values",p,heading=3)
      getHeatmap(expressionMatrix,paste(prefix.file,'_expression_values_heatmap',sep=""),probe_set_names,colHeatmap)
      hwriteImage(paste(outdir,"/",prefix.file,'_expression_values_heatmap.png',sep=""),p,center=TRUE)
      hwrite(paste("Figure ",nb_figure,": Heatmap representing the genes (y axis) in each condition (x axis), the color correspond to the expresion level (red to blue, i.e low to hight values)", sep=""),p,br=TRUE)
      nb_figure=nb_figure+1

    }
    if(heatmap.profile){
      hwrite("Profiles; 0,1 or -1 values",p,heading=3)
      getProfileHeatmap(L.specific.result$M.specific,paste(prefix.file,'_profile_heatmap',sep=""),colHeatmap)
      hwriteImage(paste(outdir,"/",prefix.file,'_profile_heatmap.png',sep=""),p,center=TRUE)
      hwrite(paste("Figure ",nb_figure,": Profile heatmap representing the gene specificity (y axis) over the conditions (x axis), the colors correspond to the specific level; red, white, blue: down (-1), non (0), up (1) specific or just non (0) and specific (1) (red and blue respectively) if only two colors are present", sep=""),p,br=TRUE)
      nb_figure=nb_figure+1
    }
    
    ##   More elaborated heatmap for profiles:
    ##   library('fields')
    if(heatmap.unique.profile){
      L_profile_unique=getProfile(L.specific.result$M.specific)
      col_p=colorRampPalette(brewer.pal(10,"RdBu"))(max(as.numeric(L_profile_unique[[3]][,2])))
      col_p_heatmap=col_p[as.numeric(L_profile_unique[[3]][,2])]
      print(length(col_p_heatmap))
      getHeatmap(L_profile_unique[[2]],paste(prefix.file,'_unique_profile_heatmap',sep=""),colors=colHeatmap,colRowSide=col_p_heatmap,colrange=col_p)
      hwrite("Unique profiles",p,heading=3)
      hwriteImage(paste(outdir,"/",prefix.file,"_unique_profile_heatmap.png",sep=""),p,center=TRUE)
      hwrite(paste("Figure ",nb_figure,": Profile heatmap representing the group of genes having the same profile (specific pattern) (y axis) over the condition (x axis), the colors correspond to the specific level; red, white, blue: down (-1), non (0), up (1) specific or just non (0) and specific (1) (red and blue respectively) if only two colors are present. The color bar on the left indicate the number of genes present in each group (row)", sep=""),p,br=TRUE)
      nb_figure=nb_figure+1
    }
    
    hwrite("Specific genes table result",p,heading=2)
    M_result_probeset=getSelectiveResultTable(L.specific.result$M.specific,M.specific.sum.row)
    linkNames=rep("NA",nrow(M_result_probeset))
    M_result_probeset=cbind(M_result_probeset,linkNames)
    colnames(M_result_probeset)=c("Gene","Type","# condition specific","# condition up","# condition down","Condition up","Condition down","Gene_page")
    pNAName=paste(prefix.file,"_NA_page.html",sep="")
    pNA=openPage(pNAName)
    hwrite("Sorry, the page required has not been previously generated",pNA,p,heading=2,style='color:#FE0202',br=TRUE)
      hwrite("You can generate this page(s) with the \"getGeneHtmlPage()\" function. Then use the output of the function as the \"gene.page.info\" argument of the \"getFullHtmlSpeCondResult()\" function (function used to generate the previous page).",pNA,br=TRUE,style="color:#0000ff")
        hwrite("",pNA,br=TRUE)
        hwrite("Back to full result page",pNA,link=paste(olddir,"/",prefix.file,"_",page.name,".html",sep=""), table.style='font-weight:bold',br=TRUE)
        hwrite("",pNA,br=TRUE)
    hwrite("",pNA,br=TRUE)
    rversion = sessionInfo()$R.version$version.string
    hwrite(paste("Page created by the ", hwrite("SpeCond",link="http://www.bioconductor.org/packages/release/bioc/")," package using hwriter under ",rversion,collapse=""), pNA, style='font-size:8pt')

        closePage(pNA)
    
    if(!is.null(gene.page.info)){
      ##       Some gene result pages have been created and will be linked from the gene result table
      geneLink=as.vector(gene.page.info$geneLink[match(rownames(M_result_probeset),gene.page.info$geneLink),2])
      geneLinkNA=geneLink
      geneLink=paste(gene.page.info$genePageDirectory,"/",geneLink,sep="")
      M_result_probeset[!is.na(geneLinkNA),"Gene_page"] <- "gene_page"
      geneLink[is.na(geneLinkNA)] <- paste(outdir,"/",pNAName,sep="")
      geneLinkTable=geneLink
      
      h_link=hwrite(M_result_probeset[,"Gene_page"],link=geneLinkTable,target='_blank',table=FALSE)
      h_table_link=hwrite(c("Gene page",h_link),dim=c(nrow(M_result_probeset)+1,1),row.bgcolor="#ffffaa",br=TRUE)
      h_table_result=hwrite(M_result_probeset[,c(1:7)],br=TRUE,center=TRUE, row.names=FALSE, row.bgcolor='#ffffaa',table.style='border-collapse: collapse; text-align:center;',row.style=list('text-align:center'))
      hwrite(c(h_table_result,h_table_link),p,border=0)
      hwrite(paste("Table ",nb_table,": Specific gene result: list of all the specific genes detected, the number of condition in which they are specific (all, up and down) and the respective condition names", sep=""),p,br=TRUE)
      nb_table=nb_table+1
      hwrite("",p,br=TRUE)
    }    
    else{
      h_link=hwrite(M_result_probeset[,"Gene_page"],link=rep(paste(outdir,"/",pNAName,sep=""),nrow(M_result_probeset)),target='_blank',table=FALSE)
      h_table_link=hwrite(c("Gene page*",h_link),dim=c(nrow(M_result_probeset)+1,1),row.bgcolor="#ffffaa",table.style='border-collapse: collapse; text-align:center;',br=TRUE)
      h_table_result=hwrite(M_result_probeset[,c(1:7)],center=TRUE, row.names=FALSE, row.bgcolor='#ffffaa',table.style='border-collapse: collapse; text-align:center;',row.style=list('text-align:center'),col.style=list('text-align:center'),br=TRUE)
      hwrite(c(h_table_result,h_table_link),p,border=0)
      hwrite(paste("Table ",nb_table,": Specific gene result: list of all the specific genes detected, the number of condition in which they are specific (all, up and down) and the respective condition names", sep=""),p,br=TRUE)
      nb_table=nb_table+1
      hwrite("",p,br=TRUE)
      hwrite("*, ",p)
      hwrite("NA ???: ",p,style='color:#FF8A00')
      hwrite("You can generate the gene result page(s) of the gene present in the above table with the \"getGeneHtmlPage()\" function. Using its result object as the  \"gene.page.info\" argument of the \"getFullHtmlSpeCondResult()\" function (function generating this page), you will have a direct link to the gene pages",p,br=TRUE,style='color:#0000ff')
        hwrite("",p,br=TRUE)
    }
    hwrite(paste("The corresponding file is : ",outdir,"/",prefix.file,"_result_specific_probeset.txt",sep=""),p,br=TRUE)
    
    ##   Write this gene result table in the file result_specific_probeset.txt
    write.table(M_result_probeset,file=paste(prefix.file,"_result_specific_probeset.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
    print(paste(prefix.file,"_result_specific_probeset.txt",sep=""))
  }
  else{
    hwrite("No probe set have been detected as specific in this experiment",p,heading=3,style='color:#FE0202')
    }
  hwrite("",p,br=TRUE)
  hwrite("",p,br=TRUE)
  rversion = sessionInfo()$R.version$version.string
  hwrite(paste("Page created by the ", hwrite("SpeCond",link="http://www.bioconductor.org/packages/release/bioc/")," package using hwriter under ",rversion,collapse=""), p, table.style='font-size:8pt')

  closePage(p)
  setwd(olddir)
}

