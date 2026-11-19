bmpd.big<-function(comm, pd.desc="pd.desc", pd.spname, pd.wd,
                   spname.check=FALSE, abundance.weighted = TRUE,
                   na.zero=TRUE, unit.sum=NULL)
{
  requireNamespace("bigmemory")
  if(spname.check)
  {
    check.sp=match.name(name.check = pd.spname,cn.list = list(comm=comm))
    comm=check.sp$comm
  }
  pdbig.id=match(colnames(comm),pd.spname)
  comt=comm
  if(!abundance.weighted){comt[comt>0]=1}
  N=nrow(comm)
  
  if(is.null(unit.sum))
  {
    rs.comt=rowSums(comt)
    if(na.zero){rs.comt[rs.comt==0]=1}
  }else{
    rs.comt=unit.sum
  }
  comt=comt/rs.comt
  comt=as.matrix(comt)
  #pd=bigmemory::attach.big.matrix(dget(pd.desc))
  pd=try(bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc))))
  if(inherits(pd,"try-error")){pd=bigmemory::attach.big.matrix(paste0(pd.wd,"/",pd.desc))}
  comd=comt
  if(nrow(pd)<=45000)
  {
    comd=comt %*% pd[pdbig.id,pdbig.id]
  }else{
    N.int=floor(45000*45000/nrow(pd))
    N.num=ceiling(ncol(pd)/N.int)
    N.ser=((1:N.num)-1)*N.int+1
    N.ser=c(N.ser,ncol(pd)+1)
    comd=comt
    for(i in 1:N.num)
    {
      comd[,N.ser[i]:(N.ser[i+1]-1)]=comt %*% pd[pdbig.id,pdbig.id[N.ser[i]:(N.ser[i+1]-1)]]
    }
  }
  gc()
  res=comd %*% t(comt)
  res=(res+t(res))/2
  stats::as.dist(res)
}
