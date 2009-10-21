setMethod("show", signature( object= "gene_list"),
          function(object) {
            cat("An object of class \"",class(object),"\"\n",sep="")
            printHeadListGene(object)
            cat("\n")
          }
          )

setMethod("show", signature( object= "sp_list"),
          function(object) {
            if("specificResult" %in% names(object)){
              ##if the generalResult object from SpeCond() is called
              cat("An object of class \"",class(object),"\"\n",sep="")              
              for (what in names(object)) {
                cat("$",what,"\n",sep="")
              }
              print(object$specificResult)
            }
            else{
              cat("An object of class \"",class(object),"\"\n Only the main values are presented here see show.sp_list() function for \na comprehensive view of this object\n",sep="")
              for (what in names(object)) {
                if(what=="prefix.file"){
                  cat("$",what,"\n",sep="")
                  print(object[[what]])
                }
              if(what=="param.detection"){
                cat("$",what,"\n",sep="")
                print(object[[what]])
              }
                if(what=="L.specific.result"){
                  x <- object[[what]]
                  for (what2 in names(x)) {
                    if(what2=="M.specific"){
                      cat("$",what,"$",what2,"\n",sep="")
                      print(what2)
                      y <- x[[what2]]
                      printHead(y)
                      cat("\n")
                  }
                    if(what2=="M.specific.sum.row" | what2=="M.specific.sum.column"){
                      cat("$",what,"$",what2,"\n",sep="")
                      y <- x[[what2]]
                      printHead(y)
                      cat("\n")
                    }
                  }
                }
              }              
              ##             General results
              cat("\n")
              cat("Number of genes evaluated:",length(object$L.specific.result$specific),"\n",sep=" ")
              cat("Number of genes specific:",nrow(object$L.specific.result$M.specific),"\n",sep=" ")
              cat("Range of the number of conditions for which a gene have been detected as specific:",range(object$L.specific.result$M.specific.sum.row),"\n",sep=" ")
              ##             ???
              cat("Number of genes specific by number of specific conditions\n")
              M.specific.sum.row=apply(abs(object$L.specific.result$M.specific),1,sum)
              tab=t(as.matrix(table(M.specific.sum.row)))
              specific_table_nb_tissues=matrix(c(colnames(tab),as.vector(tab)),nrow=2,ncol=ncol(tab),byrow=TRUE)
              rownames(specific_table_nb_tissues)=c("# conditions","# genes")
              specific_table_nb_tissues=cbind(specific_table_nb_tissues,c("sum",sum(tab)))
              print(specific_table_nb_tissues)
            }
          }
          )

show.sp_list <- function(object) {
  cat("An object of class \"",class(object),"\"\n",sep="")
  for (what in names(object)) {
    x <- object[[what]]
    cat("$",what,"\n",sep="")
    if(is(x,"gene_list")){
      print(x)
    }
    else{
      printHead(x)
      cat("\n")
    }
  }
}

printHeadList <- function(l){
  ## print the first four element of the attributes of a list
  if(length(l)<10){
    for (what in names(l)) {
      x <- l[[what]]
      cat("$",what,"\n",sep="")
      printHead(x)
      cat("\n")
    }
  }
  else{
    printHeadListGene(l)
  }
}

printHeadListGene <- function(lg){
  print("Values linked to the list of genes")
  if(length(lg)>10)
    {
      print("display the first 3 attributes")
      
      for (what in names(lg[1:3])){
        x <- lg[[what]]
        cat("$",what,"\n",sep="")
        print(x)
      }
      cat(length(lg)-3,"more elements ...\n") 
    }
  else{
    printHeadList(lg)
  }
}

printHeadVector <- function(v){
  n <- length(v)
  if(n > 20) {
    print(v[1:5])
    cat(n-5,"more elements ...\n")  
  }
  else{
    print(v)
  }
}

printHeadMatrix <- function(m){

  nr=nrow(m)
  nc=ncol(m)
  if(nr > 20 & nc>5) {
    print(m[1:5,1:5])
    cat(nr-5,"more rows and",nc-5,"more columns ...\n")  
  }
  else{
    if(nr > 20){
      print(m[1:5,])
      cat(nr-5,"more rows ...\n")
    }
    else{
      if(nc > 5){
        print(m[,1:5])
        cat(nc-5,"more columns ...\n")
      }
      else{
        print(m)
      }
    }
  }
}

printHead <- function(object){
  
  if(is.list(object)){
    printHeadList(object)
  }
  else{
    if(is.matrix(object)){
      printHeadMatrix(object)
    }
    else{
      if(is.vector(object)){
        printHeadVector(object)
      }
      else{
        if(is.null(object)){
          print("NULL")
        }
      }
    }
  }
}
