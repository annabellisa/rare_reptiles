# Polygon CI function
# Author: Annabel Smith

pg.ci<-function(x,data,x.subset=NULL,colour){
  
  # DEFINITION:
  # x, data and x.subset must be in quotes
  # x: a vector of data on the x-axis
  # data: a data frame which includes x 
  # currently, the confidence intervals must be specified as $lci and $uci, need to fix this... 
  # x.subset: a vector of data for subsetting x. If there no subsets, omit this argument or set it at NULL to plot a single polygon. If there are subsets, enter the colname of data frame x to be used as the subset
  
  # update 13th Feb 2017: I haven't fixed the $uci and $lci problem, but I've just added in an if to deal with two different cases that distinguish $uci.resp, etc.
  # update 14th Feb 2024: This code is a mess with lots of hacks; I've added more ifs to deal with the CI notation in iNEXT. I need to add extra arguments for the CI notation... but don't have time to test it on different data sets.
  
  # No subsets:
  if(is.null(x.subset)==T){
    xx<-paste(data,"$",x,sep="")
    # lci<- paste(data,"$lci",sep="")
    # uci<- paste(data,"$uci",sep="")
    if(length(grep(".resp",colnames(get(data))))>0) lci<-paste(data,"$lci.resp",sep="") else lci<-paste(data,"$lci",sep="")
    if(length(grep(".resp",colnames(get(data))))>0) uci<-paste(data,"$uci.resp",sep="") else uci<-paste(data,"$uci",sep="")
    
    xvec <- c(eval(parse(text=xx)), tail(eval(parse(text=xx)), 1), rev(eval(parse(text=xx))), eval(parse(text=xx))[1])
    yvec <- c(eval(parse(text=lci)), tail(eval(parse(text=uci)), 1), rev(eval(parse(text=uci))), eval(parse(text=lci))[1])
    polygon(xvec, yvec, col=colour, border=NA)
  } # close if no subsets
  
  # with subsets
  if(is.null(x.subset)==F){
    
    if(length(grep(".LCL",colnames(get(data))))>0) lci<-paste(data,"$qD.LCL",sep="") else uci<-paste(data,"$uci",sep="")
    if(length(grep(".LCL",colnames(get(data))))>0) uci<-paste(data,"$qD.UCL",sep="") else uci<-paste(data,"$uci",sep="")
    
    # Get data and vector that is used for subsetting:
    data.withsubset<-get(data)
    subset.all<-data.withsubset[,x.subset]
    
    
    # Specify subs.levs: levels for factors, unique numbers for binary variables, and first and third quartiles for continuous variables
    
    if(is.factor(subset.all)) subs.levs<-levels(subset.all)
    
    if(is.factor(subset.all)==F) {
      
      if(length(unique(subset.all))==2) subs.levs<-unique(subset.all) 

    } # close if subset is not a factor
    
    for (i in 1:length(subs.levs)){
      
      sub.thisrun<-subs.levs[i]
      x.thisrun<-data.withsubset[which(subset.all==sub.thisrun),x]
      # lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"lci"]
      # uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"uci"]
     
      if(length(grep("qD.LCL",colnames(data.withsubset)))>0) lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"qD.LCL"] else lci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"lci"]
      if(length(grep("qD.UCL",colnames(data.withsubset)))>0) uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"qD.UCL"] else uci.thisrun<-data.withsubset[which(subset.all==sub.thisrun),"uci"]
      
      xvec <- c(x.thisrun, tail(x.thisrun, 1), rev(x.thisrun), x.thisrun[1])
      yvec <- c(lci.thisrun, tail(uci.thisrun, 1), rev(uci.thisrun), lci.thisrun[1])
      polygon(xvec, yvec, col=colour, border=NA)
      
    } # close for sub levels i
    
  } # close if subset present
  
} # close polygon function
