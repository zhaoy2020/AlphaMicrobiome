dniche<-function(env,comm,method=c("ab.overlap","niche.value","prefer.overlap"),
                 nworker=4,memory.G=50,out.dist=FALSE,bigmemo=TRUE,nd.wd=getwd(),
                 nd.spname.file="nd.names.csv",detail.file="ND.res.rda",overalldist=FALSE)
{
  checksamp=match.name(name.check = rownames(env),rn.list = list(env=env,comm=comm))
  env=checksamp$env
  comm=checksamp$comm
  envnum=ncol(env)
  spname=colnames(comm)
  if(is.null(colnames(env))){colnames(env)=paste0("env",1:ncol(env))}
  if(is.null(nd.wd)){nd.wd=getwd()}else{if(!dir.exists(nd.wd)){dir.create(nd.wd)}}
  
  res=list()
  res$bigmemo=bigmemo
  res$nd=list()
  if(bigmemo)
  {
    requireNamespace("bigmemory")
    res$nd.wd=nd.wd
    res$names=spname
    utils::write.csv(data.frame(names=spname),file = paste0(nd.wd,"/",nd.spname.file))
    ndbig=list()
    for(j in 1:ncol(env))
    {
      ndbig[[j]] = bigmemory::big.matrix(nrow = length(spname), ncol = length(spname), type = "double",
                                         backingfile = paste0(colnames(env)[j],".ND.bin"),
                                         backingpath = nd.wd,
                                         descriptorfile = paste0(colnames(env)[j],".ND.desc"),
                                         shared = TRUE)
      ndbig[[j]][]=0
      res$nd[[j]]=paste0(colnames(env)[j],".ND.desc")
    }
    if(overalldist)
    {
      j=ncol(env)+1
      ndbig[[j]]=bigmemory::big.matrix(nrow = length(spname), ncol = length(spname), type = "double",
                                       backingfile = "ALL.ND.bin",backingpath = nd.wd,
                                       descriptorfile = "ALL.ND.desc",shared = TRUE)
      ndbig[[j]][]=0
      res$nd[[j]]="ALL.ND.desc"
    }
  }
  if(overalldist)
  {
    env.old=env
    env=data.frame(scale(env.old))
  }
  if(method[1]=="niche.value")
  {
    comts<-(t(comm)/colSums(comm))
    nv<-(as.matrix(comts) %*% as.matrix(env))
    for(j in 1:ncol(env))
    {
      ndjm=as.matrix(stats::dist(nv[,j]))
      if(bigmemo)
      {
        maxndjm=max(ndjm)
        gc()
        if(nrow(ndjm)>20000)
        {
          ndjmxx=ceiling((20000^2)/nrow(ndjm))
          ndjmxsq=seq(from=1,to=nrow(ndjm),by=ndjmxx)
          for(x in 1:(length(ndjmxsq)-1))
          {
            message("ND matrix is too big, now saving j=",j," in ",ncol(env),", x=",x," in ",length(ndjmxsq),". ",date())
            ndbig[[j]][,ndjmxsq[x]:(ndjmxsq[x+1]-1)] = ndjm[,ndjmxsq[x]:(ndjmxsq[x+1]-1)]/maxndjm
            if(overalldist){ndbig[[ncol(env)+1]][,ndjmxsq[x]:(ndjmxsq[x+1]-1)]=ndbig[[ncol(env)+1]][,ndjmxsq[x]:(ndjmxsq[x+1]-1)]+((ndjm[,ndjmxsq[x]:(ndjmxsq[x+1]-1)])^2)}
            gc()
          }
          ndbig[[j]][,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)] = ndjm[,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)]/maxndjm
          if(overalldist){ndbig[[ncol(env)+1]][,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)]=ndbig[[ncol(env)+1]][,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)]+((ndjm[,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)])^2)}
          gc()
        }else{
          ndbig[[j]][] = ndjm/maxndjm
          if(overalldist){ndbig[[ncol(env)+1]][]=ndbig[[ncol(env)+1]][]+(ndjm^2)}
        }
      }else{
        res$nd[[j]]=ndjm/max(ndjm)
        if(overalldist){if(j==1){res$nd[[ncol(env)+1]]=(ndjm^2)}else{res$nd[[ncol(env)+1]]=res$nd[[ncol(env)+1]]+(ndjm^2)}}
      }
    }
  }else if(method[1] %in% c("ab.overlap","prefer.overlap")){
    if(method[1]=="ab.overlap")
    {
      comp=comm
    }else{
      comp=t(t(comm)/colSums(comm))
    }
    dens<-function(i,envj,comp)
    {
      stats::density(envj,weights = comp[,i],from=min(envj),to=max(envj))$y
    }
    dio<-function(i,den,res,j)
    {
      den.b<-den.max<-den[,(i+1):ncol(den),drop=FALSE]
      den.a<-den.am<-matrix(den[,i],nrow=nrow(den),ncol=ncol(den.b))
      den.max[den.max<den.a]=0
      den.am[den.am<=den.b]=0
      den.max=den.max+den.am
      den.dif=abs(den.b-den.a)
      ndio=colSums(den.dif)/colSums(den.max)
      if(res$bigmemo)
      {
        requireNamespace("bigmemory")
        ndj=bigmemory::attach.big.matrix(dget(res$nd[[j]]))
        ndj[i,i]=0
        ndj[i,(i+1):ncol(den)]=ndio
        ndj[(i+1):ncol(den),i]=ndio
        out=i
      }else{
        out=rep(0,ncol(den))
        out[(i+1):ncol(den)]=ndio
      }
      out
    }
    
    requireNamespace("parallel")
    if(.Platform$OS.type=="windows")
    {
      if(utils::memory.limit()<memory.G*1024)
      {
        memotry=try(utils::memory.limit(size=memory.G*1024),silent = TRUE)
        if(inherits(memotry,"try-error")){warning(memotry[1])}
      }
    }
    
    for(j in 1:envnum)
    {
      dens1=stats::density(env[,j],weights = comp[,1],from=min(env[,j]),to=max(env[,j]))
      den1=data.frame(dens1$y)
      rownames(den1)<-dens1$x
      if(nworker==1)
      {
        message("Now computing density model. j=",j," in ",envnum,". begin at ", date(),". Please wait...")
        den<-lapply(2:ncol(comm),dens,envj=env[,j],comp=comp)
      }else{
        c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
        if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
        if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
        message("Now parallel computing density model. j=",j," in ",envnum,". begin at ", date(),". Please wait...")
        den<-parallel::parLapply(c1,2:ncol(comm),dens,envj=env[,j],comp=comp)
        parallel::stopCluster(c1)
      }
      
      den=matrix(unlist(den),nrow=length(den[[1]]))
      den=as.matrix(cbind(den1,den))
      gc()
      
      if(nworker==1)
      {
        message("Now computing niche distance. j=",j," in ",envnum,". begin at ", date(),". Please wait...")
        dis<-lapply(1:(ncol(den)-1),dio,den=den,res=res,j=j)
      }else{
        c0<-try(parallel::makeCluster(nworker,type="PSOCK"))
        if(inherits(c0,"try-error")){c0 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
        if(inherits(c0,"try-error")){c0 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
        message("Now parallel computing niche distance. j=",j," in ",envnum,". begin at ", date(),". Please wait...")
        dis<-parallel::parLapply(c0,1:(ncol(den)-1),dio,den=den,res=res,j=j)
        parallel::stopCluster(c0)
      }
      if(!bigmemo)
      {
        dis=data.frame(matrix(unlist(dis),nrow=length(dis[[1]])))
        dis=cbind(dis,rep(0,ncol(comm)))
        dis=dis+t(dis)
        colnames(dis)<-rownames(dis)<-spname
        if(out.dist){dis=stats::as.dist(dis)}
        res$nd[[j]]=dis
      }
    }
    
    if(overalldist)
    {
      for(j in 1:envnum)
      {
        if(bigmemo)
        {
          if(envnum>20000)
          {
            ndjmxx=ceiling((20000^2)/envnum)
            ndjmxsq=seq(from=1,to=envnum,by=ndjmxx)
            for(x in 1:(length(ndjmxsq)-1))
            {
              message("ND matrix is too big, now saving j=",j," in ",ncol(env),", x=",x," in ",length(ndjmxsq),". ",date())
              ndbig[[ncol(env)+1]][,ndjmxsq[x]:(ndjmxsq[x+1]-1)]=ndbig[[ncol(env)+1]][,ndjmxsq[x]:(ndjmxsq[x+1]-1)]+((ndbig[[j]][,ndjmxsq[x]:(ndjmxsq[x+1]-1)])^2)
              gc()
            }
            ndbig[[ncol(env)+1]][,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)]=ndbig[[ncol(env)+1]][,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)]+((ndbig[[j]][,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)])^2)
            gc()
          }else{
            ndbig[[ncol(env)+1]][]=ndbig[[ncol(env)+1]][]+(ndbig[[j]][]^2)
          }
        }else{
          if(j==1){res$nd[[ncol(env)+1]]=(res$nd[[j]]^2)}else{res$nd[[ncol(env)+1]]=res$nd[[ncol(env)+1]]+(res$nd[[j]]^2)}
        }
      }
    }
  }
  
  if(overalldist)
  {
    if(bigmemo)
    {
      if(nrow(ndjm)>20000)
      {
        for(x in 1:(length(ndjmxsq)-1))
        {
          ndbig[[ncol(env)+1]][,ndjmxsq[x]:(ndjmxsq[x+1]-1)]=(ndbig[[ncol(env)+1]][,ndjmxsq[x]:(ndjmxsq[x+1]-1)])^0.5
          gc()
        }
        ndbig[[ncol(env)+1]][,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)]=(ndbig[[ncol(env)+1]][,ndjmxsq[length(ndjmxsq)]:nrow(ndjm)])^0.5
        gc()
      }else{
        ndbig[[ncol(env)+1]][]=(ndbig[[ncol(env)+1]][])^0.5
      }
    }else{
      res$nd[[ncol(env)+1]]=(res$nd[[ncol(env)+1]])^0.5
    }
  }
  if(overalldist){names(res$nd)=c(colnames(env),"ALL")}else{names(res$nd)=colnames(env)}
  res$method=method[1]
  if(bigmemo){if(!is.null(detail.file)){save(res,file = paste0(nd.wd,"/",detail.file))}}
  res
}