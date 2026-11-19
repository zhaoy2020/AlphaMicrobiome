metaset<-function(comm=NULL,meta.group=NULL,meta.com=NULL,meta.frequency=NULL,meta.ab=NULL,tree=NULL)
{
  # version 20241222A
  if(is.null(comm))
  {
    if(!is.null(meta.com))
    {
      if(!is.list(meta.com)){meta.com=list(Meta=meta.com)}
      if(is.null(names(meta.com))){names(meta.com)=paste0("Meta",1:length(meta.com))}
      spall=unique(unlist(lapply(meta.com,colnames)))
      if(!is.null(meta.frequency)){warning("Both meta.com and meta.frequency are given. meta.frequence is re-calculated based on meta.com, and the given meta.frequence is discarded.")}
      if(!is.null(meta.ab)){warning("Both meta.com and meta.ab are given. meta.ab is re-calculated based on meta.com, and the given meta.ab is discarded.")}
      meta.frequency=t(sapply(1:length(meta.com),function(i){outi=rep(0,length(spall));outi[match(colnames(meta.com[[i]]),spall)]=colSums(meta.com[[i]]>0);outi}))
      meta.ab=t(sapply(1:length(meta.com),function(i){outi=rep(0,length(spall));outi[match(colnames(meta.com[[i]]),spall)]=colSums(meta.com[[i]]);outi}))
      rownames(meta.frequency)<-rownames(meta.ab)<-names(meta.com)
      colnames(meta.frequency)<-colnames(meta.ab)<-spall
    }
  }else{
    if(!is.null(meta.group))
    {
      sampc=iCAMP::match.name(rn.list = list(comm=comm,meta.group=meta.group))
      comm=sampc$comm
      meta.group=sampc$meta.group
      meta.lev=unique(meta.group[,1])
    }else{
      meta.group=data.frame(metagrp=rep("Meta",nrow(comm)),stringsAsFactors = FALSE)
      rownames(meta.group)=rownames(comm)
      meta.lev="Meta"
    }
    comms=lapply(meta.lev,function(mi){sampi=rownames(meta.group)[which(meta.group[,1]==mi)];comi=comm[which(rownames(comm) %in% sampi),,drop=FALSE];comi[,colSums(comi)>0,drop=FALSE]})
    names(comms)=meta.lev
    ##########################
    if((!is.null(meta.com)) & ((is.null(meta.frequency)) | (is.null(meta.ab))))
    {
      if(!is.list(meta.com))
      {
        if(length(meta.lev)>1){stop("meta.group has more than one metacommunities, but meta.com only defines one.")}
        meta.com=list(Meta=meta.com)
      }
      if(sum(!(meta.lev %in% names(meta.com)))>0)
      {
        stop('meta.com names must be the same as metacommunity names in meta.group.')
      }else{
        meta.com=meta.com[match(meta.lev,names(meta.com))]
      }
      if(sum(sapply(1:length(meta.lev),function(i){sum(!(colnames(comms[[i]]) %in% colnames(meta.com[[i]])))}))>0)
      {
        stop('comm has some species not included the meta.com of its metacommunity.')
      }
    }
    
    if(!is.null(meta.frequency))
    {
      if(sum(!(colnames(comm) %in% colnames(meta.frequency)))){stop('comm has some species not included in meta.frequence.')}
      if(!is.null(tree)){if(sum(!(colnames(meta.frequency) %in% tree$tip.label))>0){stop('meta.frequency has some species not included in tree.')}}
      if(sum(!(meta.lev %in% rownames(meta.frequency)))>0)
      {
        stop('meta.frequency rownames must be the same as metacommunity names in meta.group.')
      }else{
        meta.frequency=meta.frequency[match(meta.lev,rownames(meta.frequency)),,drop=FALSE]
      }
    }else{
      if(is.null(meta.com))
      {
        meta.frequency=matrix(0,nrow = length(meta.lev),ncol = ncol(comm))
        rownames(meta.frequency)=meta.lev
        colnames(meta.frequency)=colnames(comm)
        for(i in 1:length(comms))
        {
          meta.frequency[i,match(colnames(comms[[i]]),colnames(comm))]=colSums(comms[[i]]>0)
        }
      }else{
        spall=unique(unlist(lapply(meta.com,colnames)))
        meta.frequency=t(sapply(1:length(meta.com),
                                function(i)
                                {
                                  outi=rep(0,length(spall))
                                  outi[match(colnames(meta.com[[i]]),spall)]=colSums(meta.com[[i]]>0)
                                }))
        rownames(meta.frequency)=meta.lev
        colnames(meta.frequency)=spall
      }
    }
    
    # 
    if(!is.null(meta.ab))
    {
      if(sum(!(colnames(comm) %in% colnames(meta.ab)))){stop('comm has some species not included in meta.ab')}
      if(!is.null(tree)){if(sum(!(colnames(meta.ab) %in% tree$tip.label))>0){stop('meta.ab has some species not included in tree.')}}
      if(sum(!(colnames(meta.frequency) %in% colnames(meta.ab)))){stop('meta.frequency has some species not included in meta.ab')}
      if(sum(!(meta.lev %in% rownames(meta.ab)))>0)
      {
        stop('meta.ab rownames must be the same as metacommunity names in meta.group.')
      }else{
        meta.ab=meta.ab[match(meta.lev,rownames(meta.ab)),match(colnames(meta.frequency),colnames(meta.ab)),drop=FALSE]
      }
    }else{
      if(is.null(meta.com))
      {
        meta.ab=matrix(0,nrow = length(meta.lev),ncol = ncol(comm))
        rownames(meta.ab)=meta.lev
        colnames(meta.ab)=colnames(comm)
        for(i in 1:length(comms))
        {
          meta.ab[i,match(colnames(comms[[i]]),colnames(comm))]=colMeans(comms[[i]]/rowSums(comms[[i]]))
        }
      }else{
        spall=unique(unlist(lapply(meta.com,colnames)))
        meta.ab=t(sapply(1:length(meta.com),
                         function(i)
                         {
                           outi=rep(0,length(spall))
                           outi[match(colnames(meta.com[[i]]),spall)]=colMeans(meta.com[[i]]/rowSums(meta.com[[i]]))
                         }))
        rownames(meta.ab)=meta.lev
        colnames(meta.ab)=spall
      }
    }
  }
  if(!is.null(meta.frequency))
  {
    meta.spool=lapply(1:nrow(meta.frequency),function(i){colnames(meta.frequency)[which(meta.frequency[i,]>0)]})
    names(meta.spool)=rownames(meta.frequency)
  }else if(!is.null(meta.ab)){
    meta.spool=lapply(1:nrow(meta.ab),function(i){colnames(meta.ab)[which(meta.ab[i,]>0)]})
    names(meta.spool)=rownames(meta.ab)
  }else{meta.spool=NULL}
  list(comm=comm,meta.group=meta.group,meta.com=meta.com,meta.frequency=meta.frequency,meta.ab=meta.ab,meta.spool=meta.spool,comms=comms)
}