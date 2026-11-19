mntd.big<-function(comm, pd.desc="pd.desc", pd.spname=NULL, pd.wd=getwd(),
                   spname.check=FALSE, abundance.weighted = TRUE, 
                   silent=TRUE,memory.G=50)
{
  requireNamespace("bigmemory")
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memory.G*1024)
    {
      memotry=try(utils::memory.limit(size=memory.G*1024),silent = TRUE)
      if(inherits(memotry,"try-error")){warning(memotry[1])}
    }
  }
  if(spname.check)
  {
    if(is.null(pd.spname)) pd.spname=lazyopen(file=paste0(pd.wd,"/pd.taxon.name.csv"))[,1]
    check.sp=match.name(name.check = pd.spname,cn.list = list(comm=comm))
    comm=check.sp$comm
  }
  pdbig.id=match(colnames(comm),pd.spname)
  #pd=bigmemory::attach.big.matrix(dget(pd.desc))
  pd=try(bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc))))
  if(inherits(pd,"try-error")){pd=bigmemory::attach.big.matrix(paste0(pd.wd,"/",pd.desc))}
  N=nrow(comm)
  gc()
  if(abundance.weighted)
  {
    min.d=matrix(0,nrow = N,ncol = ncol(comm))
    for(i in 1:N)
    {
      if(!silent) message("Now MNTD community i=",i," in ",N,". ",date())
      id=which(comm[i,]>0)
      pdx=pd[pdbig.id[id],pdbig.id[id]]
      diag(pdx)=NA
      min.d[i,id]=apply(pdx,2,min,na.rm=TRUE)
      gc()
    }
    comm.p=comm/rowSums(comm)
    res=min.d * comm.p
    res=rowSums(res)
  }else{
    res=comm[,1]
    for(i in 1:N)
    {
      if(!silent) message("Now MNTD community i=",i," in ",N,". ",date())
      id=which(comm[i,]>0)
      pdx=pd[pdbig.id[id],pdbig.id[id]]
      diag(pdx)=NA
      res[i]=mean(apply(pdx,2,min,na.rm=TRUE))
    }
  }
  names(res)=rownames(comm)
  res
}