getPValueMean <- function(expression_values,param_p,null_p){
  
  G=ncol(param_p)

  if(G==1){
    pv=sapply(1:length(expression_values), function(i) min(pnorm(expression_values[i],param_p[2,1],param_p[3,1]),(1-pnorm(expression_values[i],param_p[2,1],param_p[3,1]))))
    null_mean=param_p[2,1]
  }
  
  if(G==2){
    if(sum(null_p)==2){
      pv=sapply(1:length(expression_values), function(i) min((param_p[1,1]*pnorm(expression_values[i],param_p[2,1],param_p[3,1]) + param_p[1,2] * pnorm(expression_values[i],param_p[2,2],param_p[3,2])), (1-(param_p[1,1]*pnorm(expression_values[i],param_p[2,1],param_p[3,1])+ param_p[1,2] * pnorm(expression_values[i],param_p[2,2],param_p[3,2])))))
      null_mean =(param_p[1,1]* param_p[2,1]) + (param_p[1,2] *param_p[2,2])
    }
    else{
      ##       Use the normal distribution which have the largest proportion and well separated of the other one
      n=which(null_p==1)
      pv=sapply(1:length(expression_values), function(i) min(pnorm(expression_values[i],param_p[2,n],param_p[3,n]), (1-pnorm(expression_values[i],param_p[2,n],param_p[3,n]))))
      null_mean=param_p[2,n]
    }
  }
  if(G==3){
    if(sum(null_p)==3){
      ##       Combination of the three normal distributions
      pv=sapply(1:length(expression_values), function(i) min((param_p[1,1]*pnorm(expression_values[i],param_p[2,1],param_p[3,1]) + param_p[1,2] * pnorm(expression_values[i],param_p[2,2],param_p[3,2])+ param_p[1,3]*pnorm(expression_values[i],param_p[2,3],param_p[3,3])),(1- (param_p[1,1]*pnorm(expression_values[i],param_p[2,1],param_p[3,1]) + param_p[1,2] * pnorm(expression_values[i],param_p[2,2],param_p[3,2])+ param_p[1,3]*pnorm(expression_values[i],param_p[2,3],param_p[3,3])))))
      null_mean=(param_p[1,1]* param_p[2,1]) + (param_p[1,2] *param_p[2,2]) + (param_p[1,3] *param_p[2,3])
    }
    else{
      if(sum(null_p)==2){
        ##       Combination of the two normal which have the largest proportion
        n=which(null_p==1)
        W1=param_p[1,n[1]]/(param_p[1,n[1]]+param_p[1,n[2]])
        W2=param_p[1,n[2]]/(param_p[1,n[1]]+param_p[1,n[2]])

        pv=sapply(1:length(expression_values), function(i) min((W1*pnorm(expression_values[i],param_p[2,n[1]],param_p[3,n[1]]) + W2* pnorm(expression_values[i],param_p[2,n[2]],param_p[3,n[2]])), (1-(W1*pnorm(expression_values[i],param_p[2,n[1]],param_p[3,n[1]])+ W2*pnorm(expression_values[i],param_p[2,n[2]],param_p[3,n[2]])))))
        null_mean=(W1* param_p[2,n[1]]) + (W2 *param_p[2,n[2]])
      }
      else{
        ## only one distribution from the 3 possible is kept
        n=which(null_p==1)
        pv=sapply(1:length(expression_values), function(i) min(pnorm(expression_values[i],param_p[2,n],param_p[3,n]),(1-pnorm(expression_values[i],param_p[2,n],param_p[3,n]))))
        null_mean=param_p[2,n]
      }
    }
  }

  if(G>2){
    ##       Combination of the normal distribution
    n=which(null_p==1)
    W_sum=sum(param_p[1,n])
    W=sapply(1:length(n), function(i) param_p[1,n[i]]/W_sum)
    
    pv1=sapply(1:length(expression_values), function(i) sum(W*pnorm(expression_values[i],param_p[2,n],param_p[3,n])))
    pv2=sapply(1:length(expression_values), function(i) 1-sum(W*pnorm(expression_values[i],param_p[2,n],param_p[3,n])))
    pv=sapply(1:length(pv1),function(i) min(pv1[i],pv2[i]))
    null_mean=sum(W* param_p[2,n])
  }

  mean_pv=c(null_mean,pv)
  return(mean_pv)
}

getLoglikelihoodFromBIC <- function(v_BIC,nb_condition){

  v_nb_param=c(2,5,8,11,15)
  v_nb_param=v_nb_param[1:length(v_BIC)]
  v_lk=(v_BIC+(v_nb_param*log(nb_condition)))/2
  return(v_lk)
}

getLambdaBIC <- function(lambda,v_loglik,nb_condition){

  v_nb_param=c(2,5,8,11,15)
  v_nb_param=v_nb_param[1:length(v_loglik)]
  v_bic=(2*v_loglik)-(lambda*v_nb_param*log(nb_condition))
  G=which(v_bic==max(v_bic,na.rm=TRUE))
  return(G,v_bic)
}

getScaleMAD <- function(var_param,values,G,d){

  mad_scale=var_param*(sqrt(mad(values))/G^(2/d))
  return(mad_scale)
}

getMinLoglikelihoodNull <- function(expression_values,p_param,min_loglike,id_f1,id_f2){
  l=abs(log(dnorm(expression_values,mean=p_param[2,id_f1],sd=p_param[3,id_f1]))-log(dnorm(expression_values,mean=p_param[2,id_f2],sd=p_param[3,id_f2])))
  min_l=min(l)
  min_l[is.na(min_l)]=0
  
  in_null=0
  if(min_l<min_loglike){
    in_null=1
  }
  lk_r=c(in_null,min_l)
  return(lk_r)
}

getRsdNull<- function(p_param,min_rsd,id_prop_max){

  if(ncol(p_param)==1){
    rsd_separated=0
  }
  else{
    rsd_separated=matrix(nrow=0,ncol=3)
    for(i in 1:ncol(p_param)){
      if(id_prop_max!=i){
        rsd_separated=rbind(rsd_separated,c(id_prop_max,i,(p_param[3,id_prop_max]/p_param[3,i])))
       }
    }
    in_null_boolean=rsd_separated[,3]>min_rsd
    in_null=rep(0,nrow(rsd_separated))
    in_null[in_null_boolean] <- 1
    rsd_separated=cbind(rsd_separated,in_null)
    colnames(rsd_separated)=c("position_normal_max_proportion","position_normal2","ratio_sd","in_null")
  }
  return(rsd_separated)
}


getNullLoglikelihoodRsdMd<- function(expression_values,p_param,diff_median_p,min_loglike,min_rsd,percent_selective_threshold,nb_s_Step1){

  ##combine all normal distribution
  null=rep(1,ncol(p_param))
  min_lk_values=NULL
  lk_names=vector()
  percent_Step1=(nb_s_Step1)/length(expression_values)
  percent_selective=percent_selective_threshold-percent_Step1
  small_percent=p_param[1,]<percent_selective
  max_p_id=which(p_param[1,]==max(p_param[1,]))
  
  if(sum(small_percent)>0){
    p_sort=order(p_param[2,],decreasing=FALSE)
    if(which(p_sort==max_p_id)<length(p_sort)){
      id_outliers=p_sort[which(p_sort==max_p_id)+1]
      lk_result=getMinLoglikelihoodNull(expression_values,p_param,min_loglike,max_p_id,id_outliers)
      in_null=lk_result[1]
      min_lk_values=c(min_lk_values,round(lk_result[2],3))
      lk_names=c(lk_names,paste("Normal",max_p_id,"vs Normal",id_outliers))
      null[id_outliers]=in_null
      id_out=id_outliers

      while(which(p_sort==id_outliers)<length(p_sort)){
        id_outliers=p_sort[which(p_sort==id_outliers)+1]
        lk_result=getMinLoglikelihoodNull(expression_values,p_param,min_loglike,max_p_id,id_outliers)
        in_null=lk_result[1]
        min_lk_values=c(min_lk_values,round(lk_result[2],3))
        lk_names=c(lk_names,paste("Normal",max_p_id,"vs Normal",id_outliers))
        null[id_outliers]=in_null
      }
    }
    
    if(which(p_sort==max_p_id)!=1){
      id_outliers=p_sort[which(p_sort==max_p_id)-1]
      lk_result=getMinLoglikelihoodNull(expression_values,p_param,min_loglike,max_p_id,id_outliers)
      in_null=lk_result[1]
      min_lk_values=c(min_lk_values,round(lk_result[2],3))
      lk_names=c(lk_names,paste("Normal",max_p_id,"vs Normal",id_outliers))
      null[id_outliers]=in_null
      id_out=id_outliers

      while(which(p_sort==id_outliers)>1){
        id_outliers=p_sort[which(p_sort==id_outliers)-1]
        lk_result=getMinLoglikelihoodNull(expression_values,p_param,min_loglike,max_p_id,id_outliers)
        in_null=lk_result[1]
        lk_names=c(lk_names,paste("Normal",max_p_id,"vs Normal",id_outliers))
        min_lk_values=c(min_lk_values,round(lk_result[2],3))
        null[id_outliers]=in_null
      }
    }
  }
##   Do not put null=0 to normal with a proportion larger that the percent_selective_threshold
  null[!small_percent]=1

  ## //////////////////////////////////////////////////////////////////////////////////
  ## Remove normal if ratio(sd)<min_rsd and small percentage

  rsd_values=NULL
  if(length(null)>1 && sum(small_percent>=1)){
    rsd_separated=getRsdNull(p_param,min_rsd,max_p_id)
    rsd_values=round(as.vector(rsd_separated[,3]),3)
    names(rsd_values)=paste("sd(Normal",rsd_separated[,1],")/sd(Normal",rsd_separated[,2],")",sep="")
    rsd_small_id=intersect(which(small_percent==1),rsd_separated[which(rsd_separated[,4]!=1),2])
    if(length(rsd_small_id)>1){
      rsd_percent=sum(p_param[1,rsd_small_id])
      
      while(rsd_percent>percent_selective){
        s=which(p_param[2,rsd_small_id]==min(p_param[2,rsd_small_id]))
        rsd_small_id=rsd_small_id[-s]
        rsd_percent=sum(p_param[1,rsd_small_id])
      }
    }
    null[rsd_small_id]=0
  }

  ## //////////////////////////////////////////////////////////////////////////////
  ##   Remove normal not separated if diff_median<0.75

  id_m=which(diff_median_p[,3]<=0.75)
  for(i in id_m){
    if(max_p_id %in% diff_median_p[id_m,]){
      null[diff_median_p[id_m,1]]=1
      null[diff_median_p[id_m,2]]=1
    }
  }
  ## //////////////////////////////////////////////////////////////////////////////////

  if(!is.null(min_lk_values)){
    names(min_lk_values)=lk_names
  }
  l_null_min_lk_values_rsd=list(null,min_lk_values,rsd_values)
  names(l_null_min_lk_values_rsd)=c("null","min_lk_values","sd_ratio")
  
  return(l_null_min_lk_values_rsd)
}

getDifferenceMedian <- function(expression_values,fit2_p){

  diff_median=matrix(nrow=0,ncol=3)
  if(ncol(fit2_p$NorMixParam)>1){
    for(i in 1:(ncol(fit2_p$NorMixParam)-1)){
      for(j in (i+1):ncol(fit2_p$NorMixParam)){
        diff_median=rbind(diff_median,c(i,j,abs(median(expression_values[which(fit2_p$classification==i)])-median(expression_values[which(fit2_p$classification==j)]))))
      }
    }
  }
  else{
    diff_median=matrix(c(1,1,100),1,3)
  }
  return(diff_median)
}

detectSpecificFromPV <- function(M_pv,M_signe,pv_threshold){
  
  M_s=matrix(0,nrow(M_pv),ncol(M_pv))
  M_s[M_pv<pv_threshold] <- 1 ##0,1,-1!!!!
  M_s_sum_row=apply(abs(M_s),1,sum)
  M_s_sum_column=apply(abs(M_s),2,sum)

  M_s=M_s*M_signe

##   To sort by expression values
  L_pv_s=lapply(1:nrow(M_pv), function(i) cbind(sort(M_pv[i,M_pv[i,]<pv_threshold])))
  names(L_pv_s)=rownames(M_pv)
  L_s_id_full=lapply(1:nrow(M_s), function(i) as.vector(which(!M_s[i,]==0)))
  names(L_s_id_full)=rownames(M_pv)
  s_id_length=sapply(1:length(L_s_id_full), function(i) length(L_s_id_full[[i]]))
  L_s_id=L_s_id_full[which(s_id_length!=0)]
  
  ##   Find the condition selective name if the number of condition selective is only one
  length_pv_s=sapply(1:length(L_pv_s), function(i) length(L_pv_s[[i]]))
  l1=which(length_pv_s==1)
  l1_names=sapply(1:length(l1), function(i) colnames(M_pv)[M_pv[l1[i],]<pv_threshold])
  if(length(l1)>0){
    for(j in 1:length(l1)){
      row.names(L_pv_s[[l1[j]]])=l1_names[j]
    }
  }
  ##   --
  selective=rep("Not specific",nrow(M_pv))
  s=sapply(1:length(length_pv_s), function(p) length_pv_s[[p]]>0)
  selective[s] <- "Specific"
  
  row.names(M_s)=row.names(M_pv)
  names(M_s_sum_row)=row.names(M_pv)
  names(M_s_sum_column)=colnames(M_pv)

  colnames(M_s)=colnames(M_pv)
  
  L_s=list(M_s,M_s_sum_row,M_s_sum_column,L_pv_s,selective,L_s_id)
  names(L_s)=c("M_s","M_s_sum_row","M_s_sum_column","L_pv_s","selective","L_s_id")
  return(L_s)
}

getListSelectiveResult <- function(M_s,M_s_sum_row){
  
##   M_selective contains only the probe sets detected as specific (+ or -)
  selective_probeset=rownames(M_s[M_s_sum_row>0,])
  if(length(selective_probeset)==0){
    print("No probe set is selective")
    return(NULL)
  }
  else{
    M_selective=M_s[rownames(M_s) %in% selective_probeset,]
    rownames(M_selective)=rownames(M_s)[rownames(M_s) %in% selective_probeset]
    colnames(M_selective)=colnames(M_s)
    M_selective_sum_row=apply(abs(M_selective),1,sum)
    ##   M_selective_sum_column=apply(abs(M_selective),2,sum)

    M_selective_positive=(M_selective==1)
    M_selective_positive[M_selective_positive] <- 1
    M_selective_positive_sum=apply(M_selective_positive,2,sum)
    M_selective_negative=(M_selective==-1)
    M_selective_negative[M_selective_negative] <- 1
    M_selective_negative_sum=apply(M_selective_negative,2,sum)

    M_selective_sum_column=rbind(M_selective_positive_sum,M_selective_negative_sum,(M_selective_positive_sum+M_selective_negative_sum))
    rownames(M_selective_sum_column)=c("M_specific_positive_sum","M_specific_negative_sum","Both")
    L_selective=list(M_selective,M_selective_sum_row,M_selective_sum_column)
    names(L_selective)=c("M_selective","M_selective_sum_row","M_selective_sum_column")
    return(L_selective)
  }
  print("fin getListSelectiveResult")

}

mclust_step2 <- function(probeset_name,expression_values,G_initial,fit_p){
  
  G_initial=G_initial
  G=fit_p$G
  classification=fit_p$classification
  
  if(fit_p$G==1){
    NorMixParam=matrix(c(1,fit_p$parameters$mean,sqrt(fit_p$parameters$variance$sigmasq)),nrow=3,ncol=fit_p$G,byrow=TRUE)
        colnames(NorMixParam)=c("Normal 1")
  }
  else{
    NorMixParam=matrix(c(fit_p$parameters$pro,fit_p$parameters$mean,sqrt(fit_p$parameters$variance$sigmasq)),nrow=3,ncol=fit_p$G,byrow=TRUE)
  }
  colnames_NorMix=c("Normal 1","Normal 2","Normal 3","Normal 4","Normal 5")
  colnames(NorMixParam)=colnames_NorMix[1:G]
  rownames(NorMixParam)=c("proportion","mean","sd")
  L_p=list(G_initial,G,NorMixParam,classification)
  names(L_p)=c("G_initial","G","NorMixParam","classification")
  return(L_p)
}

callMclustInStep2<- function(probeset_name,colnames,expression_values,G_initial,fit_p,beta_value=0,specific_outlier_step1){

  G_initial=G_initial
  classification=rep(10,length(expression_values))
  names(classification)=colnames
  specific_outlier_step1=specific_outlier_step1

  if(is.null(specific_outlier_step1)){
    G=fit_p$G
    classification=fit_p$classification
    
    if(fit_p$G==1){
      NorMixParam=matrix(c(1,fit_p$parameters$mean,sqrt(fit_p$parameters$variance$sigmasq)),nrow=3,ncol=fit_p$G,byrow=TRUE)
      colnames(NorMixParam)=c("Normal 1")
    }
    else{
      NorMixParam=matrix(c(fit_p$parameters$pro,fit_p$parameters$mean,sqrt(fit_p$parameters$variance$sigmasq)),nrow=3,ncol=fit_p$G,byrow=TRUE)
    }
  }
  else{
    if(beta_value==0){
      mclust2=Mclust(expression_values[-specific_outlier_step1], G = 1:3,modelNames = "V" )
    }
    else{
      print("beta!=0 strep2!!!!!!")
      mclust2=Mclust(expression_values[-specific_outlier_step1],G = 1:3, prior=priorControl(shrinkage=0,scale=getScaleMAD(beta_value,expression_values[-specific_outlier_step1],G=1:2,1)),modelNames = "V" )
    }
    
    if(mclust2$G==1){
      NorMixParam=matrix(c(1,mclust2$parameters$mean,sqrt(mclust2$parameters$variance$sigmasq)),nrow=3,ncol=mclust2$G,byrow=TRUE)
      colnames(NorMixParam)=c("Normal 1")
      classification[-specific_outlier_step1]=1 
    }
    else{
      NorMixParam=matrix(c(mclust2$parameters$pro,mclust2$parameters$mean,sqrt(mclust2$parameters$variance$sigmasq)),nrow=3,ncol=mclust2$G,byrow=TRUE)
    }
    classification[names(mclust2$classification)]=mclust2$classification
    G=mclust2$G
  }
  colnames_NorMix=c("Normal 1","Normal 2","Normal 3","Normal 4","Normal 5")
  colnames(NorMixParam)=colnames_NorMix[1:G]
  rownames(NorMixParam)=c("proportion","mean","sd")
  
  L_p=list(G_initial,G,NorMixParam,classification,specific_outlier_step1)
  names(L_p)=c("G_initial","G","NorMixParam","classification","specific_outlier_step1")
  return(L_p)
}


get_normal_separated <- function(param){
  
  L_result=list(param)
  names(L_result)=c("param")

  return(L_result)
}

changeColorClassification <- function(M_col){

  M_col_class=M_col
  M_col_class[M_col==1] <- 4
  M_col_class[M_col==2] <- 3
  M_col_class[M_col==3] <- "orange"
  M_col_class[M_col==4] <- "purple"
  M_col_class[M_col==5] <- "cyan"
  M_col_class[M_col==10] <- "black"
  return(M_col_class)

}

##Creation of the outdir directory
createOutdir <- function(outdir = getwd(), force = FALSE){

  if(file.exists(outdir)){
    if(!file.info(outdir)$isdir)
      stop(sprintf("'%s' must be a directory.", outdir))
    
    outdirContents = dir(outdir, all.files = TRUE)
    outdirContents = setdiff(outdirContents, c(".", ".."))
    
    if(!force && length(outdirContents)>0)
      stop(sprintf("'%s' is not empty.", outdir))
    if(force && length(outdirContents)>0){
      unlink(paste(outdir,"/.",sep=""), recursive=TRUE)
      print(paste("force=TRUE, Delete all files in ",outdir,sep="")) 
    }
    setwd(outdir)
    
  }
  else {
    dir.create(outdir, recursive=TRUE)
    message(sprintf("The directory '%s' has been created.", outdir))
    ##     will work in the new created directory
    setwd(outdir)
  }
}


getSelectiveResultTable <- function(M_selective,M_selective_sum_row,names_probeset=NULL){

##   M_selective contains only the probeset with at least one condition specific detection

  tissue_up=sapply(1:nrow(M_selective), function(i) paste(names(which(M_selective[i,]==1)),collapse=','))
  tissue_up[tissue_up==""] <- "-"
  names(tissue_up)=rownames(M_selective)
  tissue_up_nb=sapply(1:nrow(M_selective), function(i) length(names(which(M_selective[i,]==1))))
  
  tissue_down=sapply(1:nrow(M_selective), function(i) paste(names(which(M_selective[i,]==-1)),collapse=','))
  tissue_down[tissue_down==""] <- "-"
  names(tissue_down)=rownames(M_selective)
  tissue_down_nb=sapply(1:nrow(M_selective), function(i) length(names(which(M_selective[i,]==-1))))

  M_result_probeset=cbind(rownames(M_selective),rep("S",nrow(M_selective)),M_selective_sum_row,tissue_up_nb,tissue_down_nb,tissue_up,tissue_down)
  colnames(M_result_probeset)=c("gene","type","Nb_condition_specific","Nb_condition_up","Nb_condition_down","Condition_up","Condition_down")
  
  M_not_selective=c("N","0","0","0","-","-")
  if(!is.null(names_probeset)){
    M_result_selection_probeset=t(sapply(1:length(names_probeset), function(i) if(names_probeset[i]%in% rownames(M_selective)){ M_result_probeset[which(rownames(M_selective)==names_probeset[i]),]} else{c(names_probeset[i],M_not_selective)}))
    colnames(M_result_selection_probeset)=c("gene","type","Nb_condition_specific","Nb_condition_up","Nb_condition_down","Condition_up","Condition_down")
    return(M_result_selection_probeset)
  }
  else{
    return(M_result_probeset)
  }
}

getProfile <- function(M.specific){

  M.specific.profile=matrix(0,nrow(M.specific),2)
  rownames(M.specific.profile)=rownames(M.specific)
  colnames(M.specific.profile)=c("profile","sum")
  profile=sapply(1:nrow(M.specific), function(i) as.character(paste(M.specific[i,],collapse=",")))
  sum.row=apply(abs(M.specific),1,sum)
  M.specific.profile=cbind(profile,sum.row)

  table.profile=table(M.specific.profile[,1])
  M.specific.profile.table=cbind(names(table.profile),as.vector(table.profile))
  colnames(M.specific.profile.table)=c("profile","nb.gene")
  
  M.specific.profile.unique=matrix(as.numeric(unlist(strsplit(M.specific.profile.table[,1],','))),nrow=nrow(M.specific.profile.table),ncol=ncol(M.specific),byrow=TRUE)
  colnames(M.specific.profile.unique)=colnames(M.specific)
  L.specific=list(M.specific.profile,M.specific.profile.unique,M.specific.profile.table)
  names(L.specific)=c("M.specific.profile","M.specific.profile.unique","M.specific.profile.table")
  ##   the rows in M.specific.profile.unique and M.specific.profile.table should correspond
  return(L.specific)
}

getGeneHtmlPage <- function(expressionMatrix,specificResult,name.index.html="index.html",prefix.file=NULL, outdir="Single_result_pages", force=TRUE,gene.html=NULL,gene.html.ids=c(1:10)){

  ## specificResult contains:
##   names(specificResult)=c("prfix.file","fit","param.detection","L.specific.result","L.null","L.mlk","L.rsd","identic.row.ids")
  olddir=getwd()
  on.exit(setwd(olddir))
  
  if(is.null(prefix.file)){
    if(is(specificResult,"sp_list")){
      prefix.file=specificResult$prefix.file
    }
    else{
      stop("Error you need to enter a prefix.file value or to use a specificResult attribute of classe sp_list containing a prefix.file value (object inside the SpeCond() result value or the result of the getSpecificResult() function)")
    }
  }

  name.index.html=paste(prefix.file,"_",name.index.html,sep="")
  message(sprintf("The index html page '%s' will be created in the current directory: '%s'",name.index.html,olddir))

  outdir=paste(prefix.file,outdir,sep="_")
  
  if(nrow(expressionMatrix)<10 && !is.null(gene.html.ids)){
    gene.html.ids=c(1:nrow(expressionMatrix))
  }

  ##   Deal with the row that have not been considered as they contain only identical values
  if(!is.null(specificResult$identic.row.ids)){
    p_identic_row_names=row.names(expressionMatrix)[specificResult$identic.row.ids]
    if(is.null(gene.html) && is.null(gene.html.ids)){
      gene.html=row.names(expressionMatrix)
    }
    if(!is.null(gene.html.ids)){
      gene.html=row.names(expressionMatrix)[gene.html.ids]
    }
    probeset_html_identic=gene.html[gene.html %in% p_identic_row_names]
    gene.html=gene.html[!gene.html %in% p_identic_row_names]
    
    expressionMatrix=expressionMatrix[-(specificResult$identic.row.ids),]
    gene.html.ids=which(rownames(expressionMatrix) %in% gene.html)
  }

  M_class_col=t(sapply(1:length(specificResult$fit), function(i) changeColorClassification(specificResult$fit[[i]]$classification)))

  if(is.null(gene.html.ids) && is.null(gene.html)){
    ## ----------------------------------------
    ##     Generate html page results for all genes
    n_page=sapply(1:nrow(expressionMatrix), function(i) paste(i,"- ID ",row.names(expressionMatrix)[i],sep=""))
    n_link=sapply(1:nrow(expressionMatrix), function(i) paste(prefix.file,"_",row.names(expressionMatrix)[i],".html",sep=""))
    p=openPage(name.index.html)
    h_p=hwrite(n_page,link=paste(outdir,"/",n_link,sep=""),dim=c(length(n_page),1),border=0)
    h_s=hwrite(specificResult$L.specific.result$specific,dim=c(length(specificResult$L.specific.result$specific),1),border=0)
    hwrite(c(h_p,h_s),p,dim=c(1,2),border=0)
    rversion = sessionInfo()$R.version$version.string
    hwrite(paste("Page created by the ", hwrite("SpeCond",link="http://www.bioconductor.org/packages/release/bioc/")," package using hwriter under ",rversion,collapse=""), p, style='font-size:8pt')
    closePage(p)
    index.html.link=paste(getwd(),"/",name.index.html,sep="")

    createOutdir(outdir,force)
    print(paste("The gene html page(s) will be created in the ",outdir," directory",sep=""))

    l=sapply(1:nrow(expressionMatrix), function(i) createSingleGeneHtmlPage(index.html.link,prefix.file,i,row.names(expressionMatrix)[i],expressionMatrix[i,],specificResult$param.detection,specificResult$fit[[i]],specificResult$L.specific.result$specific[i],specificResult$fit[[i]]$NorMixParam,M_class_col[i,],specificResult$L.specific.result$M.specific.all[i,],specificResult$L.specific.result$L.pv[[i]],specificResult$L.null[[i]],specificResult$L.mlk[[i]],specificResult$L.rsd[[i]],n_link=n_link))
  }
  
  else{
    if(!is.null(gene.html)){
      gene.html.ids=which(rownames(expressionMatrix)%in% gene.html)
    }
    
    n_page=sapply(1:length(gene.html.ids), function(i) paste(i,"- ID ",row.names(expressionMatrix)[gene.html.ids[i]],sep=""))
    n_link=sapply(1:length(gene.html.ids), function(i) paste(prefix.file,"_",row.names(expressionMatrix)[gene.html.ids[i]],".html",sep=""))
    p=openPage(name.index.html)
    h_p=hwrite(n_page,link=paste(outdir,"/",n_link,sep=""),dim=c(length(n_page),1),border=0)
    h_s=hwrite(specificResult$L.specific.result$specific[gene.html.ids],dim=c(length(gene.html.ids),1),border=0)
    hwrite(c(h_p,h_s),p,dim=c(1,2),border=0)
    rversion = sessionInfo()$R.version$version.string
    hwrite(paste("Page created by the ", hwrite("SpeCond",link="http://www.bioconductor.org/packages/release/bioc/")," package using hwriter under ",rversion,collapse=""), p, style='font-size:8pt')
    closePage(p)
    index.html.link=paste(getwd(),"/",name.index.html,sep="")
    
    createOutdir(outdir,force)
    print(paste("The gene html page(s) will be created in the ",outdir," directory",sep=""))

    l=sapply(1:length(gene.html.ids), function(i) createSingleGeneHtmlPage(index.html.link,prefix.file,i,row.names(expressionMatrix)[gene.html.ids[i]],expressionMatrix[gene.html.ids[i],],specificResult$param.detection,specificResult$fit[[gene.html.ids[i]]],specificResult$L.specific.result$specific[gene.html.ids[i]],specificResult$fit[[gene.html.ids[i]]]$NorMixParam,M_class_col[gene.html.ids[i],],specificResult$L.specific.result$M.specific.all[gene.html.ids[i],],specificResult$L.specific.result$L.pv[[gene.html.ids[i]]],specificResult$L.null[[gene.html.ids[i]]],specificResult$L.mlk[[gene.html.ids[i]]],specificResult$L.rsd[[gene.html.ids[i]]],n_link=n_link))

##     To do:
##       html page for genes with identical values
  }
  currentdir=getwd()
  geneLink=cbind(rownames(expressionMatrix)[gene.html.ids],n_link)
  genePageLink=list(currentdir, geneLink)
  names(genePageLink)=c("genePageDirectory","geneLink")
  setwd("../.")
  return(genePageLink)
}

getIdenticRow <- function(expressionMatrix){
  ##   Detect rows which have identical values for each condition
  M_identic=t(sapply(1:nrow(expressionMatrix), function(i) expressionMatrix[i,]==expressionMatrix[i,1]))
  M_identic_sum=apply(M_identic,1,sum)
  identic_ids=which(M_identic_sum==ncol(expressionMatrix))
  
  if(length(identic_ids)>0){
    return(identic_ids)
  }
  else{
    return(NULL)
  }
}

getDefaultParameter <- function(){
  
  param.detection=matrix(0,nrow=2,ncol=7)
  rownames(param.detection)=c("Step 1","Step 2")
  colnames(param.detection)=c("beta","lambda","per","md","mlk","rsd","pv")
  param.detection[1,]=c(6,1,0.1,0.75,5,0.1,0.05)
  param.detection[2,]=c(0,1,0.3,0.75,25,0.1,0.05)
  return(param.detection)
}

createParameterMatrix <- function(param.detection=NULL,beta.1=NULL,beta.2=NULL,lambda.1=NULL,lambda.2=NULL,per.1=NULL,per.2=NULL,md.1=NULL,md.2=NULL,mlk.1=NULL,mlk.2=NULL,rsd.1=NULL,rsd.2=NULL,pv.1=NULL,pv.2=NULL){

  if(is.null(param.detection)){
    Mp=getDefaultParameter()
  }
  else{
    Mp=param.detection
  }
  
  if(!is.null(beta.1)){
    Mp[1,"beta"]=beta.1
  }
  
  if(!is.null(beta.2)){
    if(!beta.2==0){print("Attention: beta 2 should be equal to 0")}
    Mp[2,"beta"]=beta.2
  }
  
  if(!is.null(lambda.1)){
    Mp[1,"lambda"]=lambda.1
  }
  
  if(!is.null(lambda.2)){
    Mp[2,"lambda"]=lambda.2
  }

  if(!is.null(per.1)){
    Mp[1,"per"]=per.1
  }

    if(!is.null(per.2)){
    Mp[2,"per"]=per.2
  }
  
    if(!is.null(md.1)){
    Mp[1,"md"]=md.1
  }

    if(!is.null(md.2)){
    Mp[2,"md"]=md.2
  }

    if(!is.null(mlk.1)){
    Mp[1,"mlk"]=mlk.1
  }

    if(!is.null(mlk.2)){
    Mp[2,"mlk"]=mlk.2
  }

    if(!is.null(rsd.1)){
    Mp[1,"rsd"]=rsd.1
  }

    if(!is.null(rsd.2)){
    Mp[2,"rsd"]=rsd.2
  }

    if(!is.null(pv.1)){
    Mp[1,"pv"]=pv.1
  }

    if(!is.null(pv.2)){
    Mp[2,"pv"]=pv.2
  }
  
  if(!Mp[1,"md"]==Mp[2,"md"]){print("Attention: the two md values are different")}
  if(Mp[1,"per"]>=Mp[2,"per"]){print("Attention: per in step1 should not be larger than per in step2 must be  values are different")
                               print("per.1  per.2")
                               print(c(Mp[1,"per"],Mp[2,"per"]))
                             }
  if(!Mp[1,"pv"]==Mp[2,"pv"]){print("Attention: the two p-value thresholds are different")}
  
  return(Mp)
}

fitPrior <- function(expressionMatrix,param.detection=NULL,lambda=1,beta=6, evaluation.lambda.beta=FALSE){

  print("Step1, fitting")
  ##   Remove the identic rows
  identic_row_ids=getIdenticRow(expressionMatrix)
  if(!is.null(identic_row_ids)){
    expressionMatrix=expressionMatrix[-(identic_row_ids),]
    print(dim(expressionMatrix))
  }
  if(!is.null(param.detection)){
    lambda=param.detection[1,"lambda"]
    beta=param.detection[1,"beta"]
  }
  
  if(beta<0){
    print("Error, beta cannot be negative!!!!!!")
  }
  else{
    fit_beta_0=lapply(1:nrow(expressionMatrix), function(i) Mclust( expressionMatrix[i,], G = 1:3,modelNames = "V" ))
    if(beta==0){
      ##     No prior:
      fit=lapply(1:nrow(expressionMatrix), function(i) Mclust( expressionMatrix[i,], G = 1:3,modelNames = "V" ))
    }
    else{
      if(beta>0){
          ##     no mean prior, prior on variance
        fit=lapply(1:nrow(expressionMatrix), function(i) Mclust( expressionMatrix[i,], G = 1:3, prior=priorControl(shrinkage=0,scale=getScaleMAD(beta,expressionMatrix[i,],G=1:3,1)),modelNames = "V" ))
      }
    }
    if(lambda==1){
      fitlambda=fit
      }
    else{
      if(lambda<0){
        print("Error, lambda cannot be negative!!!!!!")
        fitlambda=fit
      }
      else{
        M_loglikelihood=t(sapply(1:length(fit), function(i) getLoglikelihoodFromBIC(fit[[i]]$BIC,ncol(expressionMatrix))))
        L_G_lambdaBIC=lapply(1:nrow(M_loglikelihood), function(i) getLambdaBIC(lambda,M_loglikelihood[i,],ncol(expressionMatrix)))
        if(beta==0){
            ##   no prior on variance
          fitlambda=lapply(1:nrow(expressionMatrix), function(i) Mclust(expressionMatrix[i,],G=L_G_lambdaBIC[[i]]$G,modelNames="V"))
        }
        else{
          fitlambda=lapply(1:nrow(expressionMatrix), function(i) Mclust(expressionMatrix[i,],G=L_G_lambdaBIC[[i]]$G,prior=priorControl(shrinkage=0,scale=getScaleMAD(beta,expressionMatrix[i,],G=L_G_lambdaBIC[[i]]$G,1)),modelNames="V"))
        }
      }
    }
##     In Step1:
    fit2=lapply(1:nrow(expressionMatrix), function(i) mclust_step2(rownames(expressionMatrix)[i],expressionMatrix[i,],G_initial=fit[[i]]$G,fitlambda[[i]]))
    names(fit2)=row.names(expressionMatrix)
    G_lambda_change=table(sapply(1:length(fit), function(i) fit[[i]]$G==fit2[[i]]$G))
    ##       print("nb of G value that did not changed from initial fit to the second one witout values from first fit with b!=0 or lambda!=1 (previous if lambda=1 it should be all true!)")
##       print(G_lambda_change)
    lambda_effect=G_lambda_change["FALSE"]
    lambda_effect[is.na(lambda_effect)]<-0
    
    G_beta_change=table(sapply(1:length(fit), function(i) fit_beta_0[[i]]$G==fit[[i]]$G))
    ##       print("nb of G value that did not changed from an initial fit with beta=0")
    ##       print(G_beta_change)
    beta_effect=G_beta_change["FALSE"]
    beta_effect[is.na(beta_effect)]<-0

    lambda_beta_effect=matrix(c(lambda_effect,beta_effect),1,2)
    colnames(lambda_beta_effect)=c("lambda","beta")
    rownames(lambda_beta_effect)="nb of genes with a change in G"
    
    ##     Create first step parameter html page

    if(evaluation.lambda.beta){
      fitResult=list(fit2,lambda_beta_effect)
      names(fitResult)=c("fit1","G.lambda.beta.effect")
    }
    else{
      fitResult=list(fit2)
      names(fitResult)=c("fit1")
    }
    return(fitResult)
  }
}

fitNoPriorWithExclusion <- function(expressionMatrix,specificOutlierStep1=FALSE,param.detection=NULL,lambda=1,beta=0){

##   To change fit2 to fit??????????????
  
  identic_row_ids=getIdenticRow(expressionMatrix)
  if(!is.null(identic_row_ids)){
    expressionMatrix=expressionMatrix[-(identic_row_ids),]
  }
  
  if(!is.null(param.detection)){
    lambda=param.detection[2,"lambda"]
    beta=param.detection[2,"beta"]
  }

  print("Step2, fitting")
  if(beta==0){
    ##     No prior:
      fit=lapply(1:nrow(expressionMatrix), function(i) Mclust( expressionMatrix[i,], G = 1:3,modelNames = "V" ))
    }
##     else{
##       if(beta>0){
##           ##     no mean prior, prior on variance
##         fit=lapply(1:nrow(expressionMatrix), function(i) Mclust( expressionMatrix[i,], G = 1:3, prior=priorControl(shrinkage=0,scale=getScaleMAD(beta,expressionMatrix[i,],G=1:3,1)),modelNames = "V" ))
##       }

  fit2=lapply(1:nrow(expressionMatrix), function(i)
    if(length(which(names(specificOutlierStep1)==rownames(expressionMatrix)[i]))>0){
      callMclustInStep2(rownames(expressionMatrix)[i],colnames(expressionMatrix),expressionMatrix[i,],G_initial=fit[[i]]$G,fit[[i]],beta_value=beta,specificOutlierStep1[[which(names(specificOutlierStep1)==rownames(expressionMatrix)[i])]])
    }
    else{
      callMclustInStep2(rownames(expressionMatrix)[i],colnames(expressionMatrix),expressionMatrix[i,],G_initial=fit[[i]]$G,fit[[i]],beta_value=beta, NULL)
    }
    )
  names(fit2)=row.names(expressionMatrix)
 
  return(fit2)

}

getSpecificOutliersStep1 <- function(expressionMatrix,fit1=NULL,param.detection=NULL, multitest.correction.method="BY", prefix.file=NULL, print.hist.pv=FALSE){
 
  if(is.null(prefix.file)){
    stop("Error: no prefix is link to this analysis")
  }

  ##   min_loglike=1,min_rsd=6,median_diff=0.75,percent_selective_threshold=0.3,pv_threshold=0.05,lambda=1,beta=0,param_Step1=NULL
  ##     param.detection[1,]=c(beta,lambda,percent_selective_threshold,median_diff,min_loglike,min_rsd,pv_threshold)
  if(is.null(fit1)){
   print("fit1 is NULL, you nust use the function fitPrior() or fitNoPriorWithExlusion() to fit the expression values before using this function")
  }
  else{
    if(is.null(param.detection)){
      param.detection=matrix(0,nrow=2,ncol=7)
      rownames(param.detection)=c("Step 1","Step 2")
      colnames(param.detection)=c("beta","lambda","per","md","mlk","rsd","pv")
      param.detection[1,]=c(6,1,0.1,0.75,5,0.1,0.05)
    }
  
    if(nrow(param.detection)==2){
      param_d=param.detection[1,]
    }
    resultSpecific=getSpecific(expressionMatrix,fit=fit1,param_d,specificOutlierStep1=NULL,multitest.correction.method="BY",prefix.file=paste(prefix.file,"_step1",sep=""),print.hist.pv=FALSE)
    L_specific_Step1=lapply(1:nrow(expressionMatrix), function(i) if(length(which(names(resultSpecific$L.specific.result$L.condition.specific.id)==rownames(expressionMatrix)[i]))>0) {as.vector(unlist(resultSpecific$L.specific.result$L.condition.specific.id[rownames(expressionMatrix)[i]]))}
      else{NULL}
      )
    names(L_specific_Step1)=rownames(expressionMatrix)
    
    return(L_specific_Step1)
  }
}

getSpecificResult<- function(expressionMatrix,fit2=NULL,param.detection=NULL,specificOutlierStep1=NULL,multitest.correction.method="BY",prefix.file=NULL,print.hist.pv=FALSE){

  if(is.null(prefix.file)){
    stop("Error: no prefix is associated to this analysis")
  }

  if(is.null(fit2)){
    print("fit2 is NULL, you nust use the function fitPrior() or fitNoPriorWithExlusion() to fit the expression values before using this function")
  }
  else{
    if(is.null(param.detection)){
      param.detection=getDefaultParameter()
    }

    if(nrow(param.detection)==2){
      param_d=param.detection[2,]
    }
    
    resultSpecific=getSpecific(expressionMatrix,fit=fit2,param.detection=param_d,specificOutlierStep1=specificOutlierStep1,multitest.correction.method="BY",prefix.file=prefix.file,print.hist.pv=print.hist.pv)
    resultSpecific$param.detection=param.detection
    return(resultSpecific)
  }
}

getSpecific <- function(expressionMatrix,fit,param.detection,specificOutlierStep1=NULL,multitest.correction.method="BY",prefix.file=NULL,print.hist.pv=FALSE,identic_row_ids=NULL){

  if(is.null(prefix.file)){
    stop("Error: no prefix is associated to this analysis:")
  }
                         
  percent_selective_threshold=param.detection["per"]
  median_diff=param.detection["md"]
  min_loglike=param.detection["mlk"]
  min_rsd=param.detection["rsd"]
  pv_threshold=param.detection["pv"]
      
  identic_row_ids=getIdenticRow(expressionMatrix)
  if(!is.null(identic_row_ids)){
    expressionMatrix=expressionMatrix[-(identic_row_ids),]
    }
  
  L_diff_median=lapply(1:nrow(expressionMatrix), function(i) getDifferenceMedian(expressionMatrix[i,],fit[[i]]))
    
  print("start: get null distributions")
  L_null_min_lk_values_rsd=lapply(1:nrow(expressionMatrix), function(p) if(length(which(names(specificOutlierStep1)==rownames(expressionMatrix)[p]))>0){
    getNullLoglikelihoodRsdMd(expressionMatrix[p,],fit[[p]]$NorMixParam,L_diff_median[[p]],min_loglike,min_rsd,percent_selective_threshold,length(specificOutlierStep1[[which(names(specificOutlierStep1)==rownames(expressionMatrix)[p])]]))
  }
  else{
    getNullLoglikelihoodRsdMd(expressionMatrix[p,],fit[[p]]$NorMixParam,L_diff_median[[p]],min_loglike,min_rsd,percent_selective_threshold,0)
  }
    )
  print("end: get null distributions")
  L_null=lapply(1:length(L_null_min_lk_values_rsd), function(p) unlist(L_null_min_lk_values_rsd[[p]]$null))
  names(L_null)=rownames(expressionMatrix)
  L_min_lk_values=lapply(1:length(L_null_min_lk_values_rsd), function(p) unlist(L_null_min_lk_values_rsd[[p]]$min_lk_values))
  names(L_min_lk_values)=rownames(expressionMatrix)
  L_rsd=lapply(1:length(L_null_min_lk_values_rsd), function(p) unlist(L_null_min_lk_values_rsd[[p]]$sd_ratio))
  names(L_rsd)=rownames(expressionMatrix)  
    
  M_null_mean_pv=t(sapply(1:nrow(expressionMatrix), function(i) getPValueMean(expressionMatrix[i,],fit[[i]]$NorMixParam,L_null[[i]])))
  M_null_mean=M_null_mean_pv[,1]
  M_signe=t(sapply(1:nrow(expressionMatrix), function(i) expressionMatrix[i,]<M_null_mean[i]))
  M_signe[M_signe] <- (-1)
  M_signe[!M_signe] <- 1
  
  M_pv_i=M_null_mean_pv[,-1]
  
  ##     /////////////////////  Histogramme pv ////////////////////////
  if(print.hist.pv){
    pdf(paste(prefix.file,"_hist_pv_notadjust.pdf",sep=""))
    hist(M_pv_i)
    dev.off()
  }
    ##     /////////////////////                ////////////////////////
  
  ## Need to correct for multiple testing
  M_pv=p.adjust(M_pv_i,method=multitest.correction.method)  
  
  row.names(M_pv)=row.names(expressionMatrix)
  colnames(M_pv)=colnames(expressionMatrix)
  
  print("start: specific detection from p-values")
  L_NorMix_s=detectSpecificFromPV(M_pv,M_signe,pv_threshold)
  print("end: specific detection from p-values")
  
  L_selective_result=getListSelectiveResult(L_NorMix_s$M_s,L_NorMix_s$M_s_sum_row)
    
  if(is.null(L_selective_result)){
    L_specific_result=list(L_NorMix_s$M_s,NULL,L_NorMix_s$M_s_sum_row, L_NorMix_s$M_s_sum_column,L_NorMix_s$L_pv_s,L_NorMix_s$selective,L_NorMix_s$L_s_id)
    names(L_specific_result)=c("M.specific.all","M.specific","M.specific.sum.row","M.specific.sum.column","L.pv","specific","L.condition.specific.id")
  }
  else{
    L_specific_result=list(L_NorMix_s$M_s,L_selective_result$M_selective,L_NorMix_s$M_s_sum_row, L_selective_result$M_selective_sum_column,L_NorMix_s$L_pv_s,L_NorMix_s$selective,L_NorMix_s$L_s_id)
    names(L_specific_result)=c("M.specific.all","M.specific","M.specific.sum.row","M.specific.sum.column","L.pv","specific","L.condition.specific.id")
  }

  specificResult=list(prefix.file,fit,param.detection,L_specific_result,L_null,L_min_lk_values,L_rsd,identic_row_ids)
  names(specificResult)=c("prefix.file","fit","param.detection","L.specific.result","L.null","L.mlk","L.rsd","identic.row.ids")
  specificResult=new("sp_list",specificResult)
  return(specificResult)
}

## Transform an ExpressionSet object to a matrux, summarising the value by using the method (mean or max) for each levels of the condition
getMatrixFromExpressionSet <- function(expSet,condition.factor=NULL,condition.method=c("mean","median","max")){
  
  if(!is(expSet,"ExpressionSet"))
    {
      stop("the expSet argument must be of type expressionSet")
    }
  else{
    Mexp=exprs(expSet)
    if(is.null(condition.factor)){
      message(sprintf("WARNING, if you have several samples for the same condition, you should consider to obtain one value by conditions (mean or maximum of the samples for example) using the argument factor.condition and method.condition to perform the SpeCond analysis"))
      print("The expressionMatrix argument that you entered has been coverted to a matrix using the exprs() function of the Biobase package")
    }
    else{
      M=Mexp
      if(!is.factor(condition.factor)){
        stop("condition must be of class factor")
      }
      else{
        if(!length(condition.factor)==ncol(M)){
          stop("the condition.factor vector size is not equal to the number of conditions in your expressionSet")
        }
        else{
          Mexp=sapply(levels(condition.factor), function(l) if(length(which(condition.factor==l))>1){
            apply(M[,which(condition.factor==l)],1,condition.method)} else{M[,which(condition.factor==l)]})
          rownames(Mexp)=rownames(M)
        }
      }
    }
  }
  return(Mexp)
}
 
SpeCond<- function(expressionMatrix,param.detection=NULL,multitest.correction.method="BY",prefix.file="A",print.hist.pv=FALSE,fit1=NULL,fit2=NULL,specificOutlierStep1=NULL,condition.factor=NULL,condition.method=c("mean","max")){

   if(is(expressionMatrix,"ExpressionSet"))
    {
##       expressionMatrix=exprs(expressionMatrix)
      expressionMatrix=getMatrixFromExpressionSet(expressionMatrix,condition.factor,condition.method)
      print("The expressionMatrix argument that you entered has been coverted to a matrix using the getMatrixFromExpressionSet() function")
    }
   else{
     if(!is(expressionMatrix,"matrix")){
       stop("The expressionMatrix argument used must be a matrix or an ExpressionSet object")
     }
   }

  if(ncol(expressionMatrix)<8){
    message(sprintf("WARNING, you are studying less than 8 conditions, this condition-specific detection tool may not be appropiate for your analysis"))
  }
  
  print("Step1")
  if(is.null(param.detection)){
##     Default parameters
    param.detection=getDefaultParameter()
  }
  if(is.null(fit1)){
    ##     With a prior b>0 to catch the single outliers
    L_b_fit=fitPrior(expressionMatrix,lambda=param.detection[1,"lambda"],beta=param.detection[1,"beta"],evaluation.lambda.beta=FALSE)
    fit1=L_b_fit$fit1
  }
  if(is.null(specificOutlierStep1)){
    specificOutlierStep1=getSpecificOutliersStep1(expressionMatrix,fit=fit1,param.detection=param.detection,multitest.correction.method="BY",prefix.file=prefix.file,print.hist.pv=FALSE)
  }
  else{
    print("Use specificOutlierStep1 as result of step1")
  }
  print("Step2")
  ## Step2
  ## Without prior to obtain normal distribution very well fitting the data
  if(is.null(fit2)){
    fit2=fitNoPriorWithExclusion(expressionMatrix,specificOutlierStep1=specificOutlierStep1,lambda=param.detection[2,"lambda"],beta=param.detection[2,"beta"])
  }
  
  specificResult=getSpecificResult(expressionMatrix,fit=fit2,specificOutlierStep1=specificOutlierStep1,param.detection=param.detection,multitest.correction.method=multitest.correction.method,prefix.file=prefix.file,print.hist.pv=print.hist.pv)

  generalResult=list(prefix.file,fit1,fit2,specificOutlierStep1,specificResult)
  names(generalResult)=list("prefix.file","fit1","fit2","specificOutliersStep1","specificResult")
  generalResult=new("sp_list",generalResult)
  return(generalResult)
}


