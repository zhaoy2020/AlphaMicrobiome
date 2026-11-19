mpd.big<-function(comm, pd.desc="pd.desc", pd.spname=NULL, pd.wd=getwd(),
                  sp.limit=10000,spname.check=FALSE,
                  abundance.weighted = TRUE)
{
  requireNamespace("bigmemory")
  if(spname.check)
  {
    if(is.null(pd.spname)) pd.spname=lazyopen(file=paste0(pd.wd,"/pd.taxon.name.csv"))[,1]
    check.sp=match.name(name.check = pd.spname,cn.list = list(comm=comm))
    comm=check.sp$comm
  }
  pdbig.id=match(colnames(comm),pd.spname)
  if(!abundance.weighted)
  {
    comm[comm>0]=1
    num=rowSums(comm)
  }
  comm=comm/rowSums(comm)
  comm=as.matrix(comm)
  #pd=bigmemory::attach.big.matrix(dget(pd.desc))
  pd=try(bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc))))
  if(inherits(pd,"try-error")){pd=bigmemory::attach.big.matrix(paste0(pd.wd,"/",pd.desc))}
  if(nrow(pd)<=sp.limit)
  {
    comd=comm %*% pd[pdbig.id,pdbig.id]
  }else{
    N.int=floor(sp.limit*sp.limit/ncol(comm))
    N.num=ceiling(ncol(comm)/N.int)
    N.ser=((1:N.num)-1)*N.int+1
    N.ser=c(N.ser,ncol(pd)+1)
    comd=comm
    for(i in 1:N.num)
    {
      message("Now begin i=",i," in ",N.num,", with memory of ",memory.size()," Mb. ",date())
      comd[,N.ser[i]:(N.ser[i+1]-1)]=comm %*% pd[pdbig.id,pdbig.id[N.ser[i]:(N.ser[i+1]-1)]]
      message("temp memory is ",memory.size()," Mb.")
      gc()
    }
  }
  gc()
  res=comd * comm
  res=rowSums(res)
  if(!abundance.weighted)
  {
    res=res*(num/(num-1))
  }
  res
}