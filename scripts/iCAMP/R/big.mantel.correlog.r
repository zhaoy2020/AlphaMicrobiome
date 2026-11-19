big.mantel.correlog<-function(x.desc,x.name,x.wd,y.desc,y.name,y.wd,
                              nworker=4,interval=0.02,break.pts=NULL,
                              permutations=999,strata=NULL,
                              diag.in=FALSE,temp.wd=getwd(),stepw=10^7,
                              padjust.method="fdr",cleanup=FALSE,parallelset=TRUE)
{
  ############################
  # 1 # load x and y
  ############################
  requireNamespace("bigmemory")
  requireNamespace("parallel")
  xm = bigmemory::attach.big.matrix(dget(paste0(x.wd, "/", x.desc)))
  ym = bigmemory::attach.big.matrix(dget(paste0(y.wd, "/", y.desc)))
  nrck=c(nrow(xm),nrow(ym),ncol(xm),ncol(ym))
  if(length(unique(nrck))!=1){stop("row number and col number should be the same for x and y")}
  nxy=nrow(xm)
  if(diag.in){nm=nxy*(nxy-1)*0.5+nxy}else{nm=nxy*(nxy-1)*0.5}
  idcxy=match(x.name,y.name)
  ############################
  # 2 # set big matrix to save data temporarily
  ############################
  (tstp=format(Sys.time(),"%Y%m%d%H%M%S"))
  bkfn=paste0("temp",tstp,".bin")
  dsfn=paste0("temp",tstp,".desc")
  bmc = bigmemory::big.matrix(nrow = nm, ncol = 8, 
                              type = "double", backingfile = bkfn, backingpath = temp.wd, 
                              descriptorfile = dsfn, shared = TRUE)
  # 2.1 # row and col ids
  k=1
  trsq=seq(from=1,to=(nxy-1),by=500)
  for(i in 1:(nxy-1))
  {
    if(i %in% trsq){message("Now set ids i=",i," in ",nxy-1,". ",date())}
    if(diag.in)
    {
      bmc[k:(k+nxy+1-i-1),1]=rep(i,nxy+1-i)
      bmc[k:(k+nxy+1-i-1),2]=i:nxy
      k=k+nxy+1-i
    }else{
      bmc[k:(k+nxy-i-1),1]=rep(i,nxy-i)
      bmc[k:(k+nxy-i-1),2]=(i+1):nxy
      k=k+nxy-i
    }
  }
  if(diag.in){bmc[nm,1]=nxy;bmc[nm,2]=nxy}
  
  # 2.2 # values
  iseq=seq(from=1,to=nm,by=stepw)
  
  if(parallelset)
  {
    setvaluex<-function(i,temp.wd,dsfn,x.wd,x.desc,y.wd,y.desc,stepw,nm,idcxy)
    {
      requireNamespace("bigmemory")
      bmc = bigmemory::attach.big.matrix(dget(paste0(temp.wd, "/", dsfn)))
      xm = bigmemory::attach.big.matrix(dget(paste0(x.wd, "/", x.desc)))
      ym = bigmemory::attach.big.matrix(dget(paste0(y.wd, "/", y.desc)))
      idi=(i:min((i+stepw-1),nm))
      bmc[idi,3]=xm[cbind(bmc[idi,1],bmc[idi,2])]
      xmax=max(bmc[idi,3])
      xmin=min(bmc[idi,3])
      bmc[idi,4]=ym[cbind(idcxy[bmc[idi,1]],idcxy[bmc[idi,2]])]
      bmc[idi,5]=bmc[idi,3]*bmc[idi,4]
      bmc[idi,6]=bmc[idi,3]^2
      bmc[idi,7]=bmc[idi,4]^2
      c(xmax,xmin)
    }
    
    message("Now set values. ",date())
    t1=Sys.time()
    setclasses=setvaluex(1,temp.wd,dsfn,x.wd,x.desc,y.wd,y.desc,stepw,nm,idcxy)
    (t2=as.numeric(difftime(Sys.time(),t1,units = "hours")))
    xmax=setclasses[1]
    xmin=setclasses[2]
    if(length(iseq)>1)
    {
      nworker1=min(nworker,length(iseq)-1)
      message("Now set values via parallel threads. may take ",round(t2*ceiling((length(iseq)-1)/nworker1),4)," hours. ",date())
      c1 <- try(parallel::makeCluster(nworker1, type = "PSOCK"))
      if (class(c1)[1] == "try-error") {c1 <- try(parallel::makeCluster(nworker1, setup_timeout = 0.5))}
      if (class(c1)[1] == "try-error") {c1 <- parallel::makeCluster(nworker1, setup_strategy = "sequential")}
      setclasses<-parallel::parSapply(c1, iseq, setvaluex,temp.wd,dsfn,x.wd,x.desc,y.wd,y.desc,stepw,nm,idcxy)
      parallel::stopCluster(c1)
      gc()
      xmax=max(setclasses[1,])
      xmin=min(setclasses[2,])
    }
  }else{
    xmax<-xmin<-xm[1,2]
    for(i in iseq)
    {
      message("Now set values i=",i," in ",nm,". ",date())
      idi=(i:min((i+stepw-1),nm))
      bmc[idi,3]=xm[cbind(bmc[idi,1],bmc[idi,2])]
      xmax=max(xmax,bmc[idi,3])
      xmin=min(xmin,bmc[idi,3])
      bmc[idi,4]=ym[cbind(idcxy[bmc[idi,1]],idcxy[bmc[idi,2]])]
      bmc[idi,5]=bmc[idi,3]*bmc[idi,4]
      bmc[idi,6]=bmc[idi,3]^2
      bmc[idi,7]=bmc[idi,4]^2
    }
  }
  # 2.3 # break points and classes
  if(is.null(break.pts))
  {
    break.pts=seq(from=(xmin %/% interval)*interval,to=xmax,by=interval)
    if(break.pts[length(break.pts)]<xmax){break.pts=c(break.pts,break.pts[length(break.pts)]+interval)}
  }else{
    break.pts=break.pts[max(which(break.pts<=xmin)):min(which(break.pts>=xmax))]
  }
  
  bkptm=matrix(break.pts[1:(length(break.pts)-1)],nrow=stepw,ncol=(length(break.pts)-1),byrow = TRUE)
  
  if(parallelset)
  {
    setclassx<-function(i,stepw,nm,temp.wd,dsfn,bkptm)
    {
      requireNamespace("bigmemory")
      bmcx = bigmemory::attach.big.matrix(dget(paste0(temp.wd, "/", dsfn)))
      idi=(i:min((i+stepw-1),nm))
      bmcx[idi,8]=rowSums(bkptm[1:length(idi),,drop=FALSE]<=bmcx[idi,3])
      i
    }
    message("Now set classes.",date())
    t1=Sys.time()
    setclassx(1,stepw,nm,temp.wd,dsfn,bkptm)
    (t2=as.numeric(difftime(Sys.time(),t1,units = "hours")))
    
    if(length(iseq)>1)
    {
      nworker2=min(nworker,length(iseq)-1)
      message("Now set classes via parallel threads. may take ",round(t2*ceiling((length(iseq)-1)/nworker2),4)," hours. ",date())
      c1 <- try(parallel::makeCluster(nworker2, type = "PSOCK"))
      if (class(c1)[1] == "try-error") {c1 <- try(parallel::makeCluster(nworker2, setup_timeout = 0.5))}
      if (class(c1)[1] == "try-error") {c1 <- parallel::makeCluster(nworker2, setup_strategy = "sequential")}
      setclasses<-parallel::parSapply(c1, iseq[2:length(iseq)], setclassx,stepw,nm,temp.wd,dsfn,bkptm)
      parallel::stopCluster(c1)
      gc()
    }
  }else{
    for(i in iseq)
    {
      message("Now set classes i=",i," in ",nm,". ",date())
      idi=(i:min((i+stepw-1),nm))
      bmc[idi,8]=rowSums(bkptm[1:length(idi),,drop=FALSE]<=bmc[idi,3])
    }
  }
  ############################
  # 3 # observed r value
  ############################
  # three methods: correlog, within-class correlation, before break point correlation
  #crlxy<-wccxy<-bbcxy<-matrix(0,nrow=length(break.pts)-1,ncol=6)
  #stepw=10^6
  #iseq=seq(from=1,to=nm,by=stepw)
  #bkptm=matrix(break.pts[1:(length(break.pts)-1)],nrow=stepw,ncol=(length(break.pts)-1),byrow = TRUE)
  bkptmc=col(bkptm)
  
  robsx<-function(i,temp.wd,dsfn,stepw,nm,break.pts,bkptmc)
  {
    requireNamespace("bigmemory")
    bmc = bigmemory::attach.big.matrix(dget(paste0(temp.wd, "/", dsfn)))
    idi=(i:min((i+stepw-1),nm))
    wc01i<-matrix(0,nrow=length(break.pts)-1,ncol=length(idi))
    wc01i[cbind(bmc[idi,8],1:length(idi))]=1
    bc01i=(bkptmc[1:length(idi),,drop=FALSE]>=bmc[idi,8])+0L
    wccxyi=cbind(sapply(3:7,function(j){wc01i %*% bmc[idi,j,drop=FALSE]}),rowSums(wc01i))
    crlxyi= cbind(wccxyi[,6],rep(sum(bmc[idi,4]),nrow(wc01i)),wccxyi[,2],
                  wccxyi[,6],rep(sum(bmc[idi,7]),nrow(wc01i)),rep(length(idi),nrow(wc01i)))
    bbcxyi=cbind(sapply(3:7,function(j){colSums(bc01i*bmc[idi,j])}),colSums(bc01i))
    cbind(crlxyi,wccxyi,bbcxyi)
  }
  
  message("Now calculate observed r values.",date())
  t1=Sys.time()
  robsm1=robsx(1,temp.wd,dsfn,stepw,nm,break.pts,bkptmc)
  (t2=as.numeric(difftime(Sys.time(),t1,units = "hours")))
  
  if(length(iseq)>1)
  {
    nworker2=min(nworker,length(iseq)-1)
    message("Now calculate observed r values via parallel threads. may take ",round(t2*ceiling((length(iseq)-1)/nworker2),4)," hours. ",date())
    c1 <- try(parallel::makeCluster(nworker2, type = "PSOCK"))
    if (class(c1)[1] == "try-error") {c1 <- try(parallel::makeCluster(nworker2, setup_timeout = 0.5))}
    if (class(c1)[1] == "try-error") {c1 <- parallel::makeCluster(nworker2, setup_strategy = "sequential")}
    robsm<-parallel::parSapply(c1, iseq[2:length(iseq)], robsx,temp.wd,dsfn,stepw,nm,break.pts,bkptmc,simplify = "array")
    parallel::stopCluster(c1)
    gc()
    robsm2=apply(robsm,c(1,2),sum)
  }else{robsm2=0}  
  robsm=robsm1+robsm2
  dim(robsm)
  crlxy=robsm[,1:6,drop=FALSE]
  wccxy=robsm[,7:12,drop=FALSE]
  bbcxy=robsm[,13:18,drop=FALSE]
  #for(i in iseq)
  #{
  #  message("Now calculate observed r values i=",i,". ",date())
  #  idi=(i:min((i+stepw-1),nm))
  #  wc01i<-matrix(0,nrow=length(break.pts)-1,ncol=length(idi))
  #  wc01i[cbind(bmc[idi,8],1:length(idi))]=1
  #  bc01i=(bkptmc[1:length(idi),,drop=FALSE]>=bmc[idi,8])+0L
  #  wccxyi=cbind(sapply(3:7,function(j){wc01i %*% bmc[idi,j,drop=FALSE]}),rowSums(wc01i))
  #  crlxyi= cbind(wccxyi[,6],rep(sum(bmc[idi,4]),nrow(wc01i)),wccxyi[,2],
  #                wccxyi[,6],rep(sum(bmc[idi,7]),nrow(wc01i)),rep(length(idi),nrow(wc01i)))
  #  bbcxyi=cbind(sapply(3:7,function(j){colSums(bc01i*bmc[idi,j])}),colSums(bc01i))
  #  crlxy=crlxy+crlxyi
  #  wccxy=wccxy+wccxyi
  #  bbcxy=bbcxy+bbcxyi
  #}
  rv<-function(mv){((mv[,6]*mv[,3])-(mv[,1]*mv[,2]))/((((mv[,6]*mv[,4])-(mv[,1]^2))^0.5)*(((mv[,6]*mv[,5])-(mv[,2]^2))^0.5))}
  r.obs=cbind(crl.obs=-rv(crlxy),
              wcc.obs=rv(wccxy),
              bcc.obs=rv(bbcxy))
  
  x.means=wccxy[,1]/wccxy[,6]
  x.sds=((1/(wccxy[,6]-1))*(wccxy[,4]-((wccxy[,1]^2)/wccxy[,6])))^0.5
  x.sds[which(wccxy[,6]<=1)]=NA
  y.means=wccxy[,2]/wccxy[,6]
  y.sds=((1/(wccxy[,6]-1))*(wccxy[,5]-((wccxy[,2]^2)/wccxy[,6])))^0.5
  y.sds[which(wccxy[,6]<=1)]=NA
  stxys=cbind(x.mean=x.means,x.stdev=x.sds,y.mean=y.means,y.stdev=y.sds)
  if(max(wccxy[,6])<(9.5*10^8))
  {
    y.qts=t(sapply(1:nrow(wccxy),
                   function(i)
                   {
                     message("Now calculating quartiles i=",i,". ",date())
                     idi=bigmemory::mwhich(bmc,8,i,'eq')
                     quantile(bmc[idi,4])
                   }))
    colnames(y.qts)=c("y.min","y.q25","y.median","y.q75","y.max")
    stxys=cbind(stxys,y.qts)
  }
  
  ############################
  # 4 # permutation
  ############################
  if(is.null(strata))
  {
    permat=NULL # if no strata, just random sampling. for big dataset, no need to use getPermuteMatrix, to save time and space.
  }else{
    permat <- vegan:::getPermuteMatrix(permutations, nxy, strata = strata)
    permutations=nrow(permat)
  }
  
  cwb.rand<-function(pj,temp.wd,dsfn,iseq,stepw,nm,nxy,break.pts,
                     permat,crlxy,wccxy,bbcxy,idcxy,y.wd,y.desc,bkptmc)
  {
    requireNamespace("bigmemory")
    bmcx = bigmemory::attach.big.matrix(dget(paste0(temp.wd, "/", dsfn)))
    ym = bigmemory::attach.big.matrix(dget(paste0(y.wd, "/", y.desc)))
    if(is.null(permat))
    {
      perid=sample(nxy,size = nxy,replace = FALSE)
    }else{
      perid=permat[pj,]
    }
    crlyr<-wccyr<-bbcyr<-matrix(0,nrow=length(break.pts)-1,ncol=3)
    
    for(i in iseq)
    {
      message("Now calculate randomized r values i=",i,". ",date())
      idix=(i:min((i+stepw-1),nm))
      idr1=perid[bmcx[idix,1]]
      idr2=perid[bmcx[idix,2]]
      yri=ym[cbind(idcxy[idr1],idcxy[idr2])]
      
      #idiy=sapply(1:length(idr1),function(k){
      #  out=bigmemory::mwhich(bmcx,1:2,list(idr2[k],idr1[k]),list('eq','eq'),'AND')
      #  if(length(out)==0){out=bigmemory::mwhich(bmcx,1:2,list(idr1[k],idr2[k]),list('eq','eq'),'AND')}
      #  out
      #})
      #wc01i<-matrix(0,nrow=length(idix),ncol=length(break.pts)-1)
      #wc01i[cbind(1:length(idix), bmcx[idix,8])]=1
      wc01i<-matrix(0,nrow=length(break.pts)-1,ncol=length(idix))
      wc01i[cbind(bmcx[idix,8],1:length(idix))]=1
      #for(j in 1:ncol(bc01i)){bc01i[which(bmcx[idix,8]<=j),j]=1}
      bc01i=(bkptmc[1:length(idix),,drop=FALSE]>=bmcx[idix,8])+0L
      
      #xyri=bmcx[idix,3]*bmcx[idiy,4]
      #y2ri=bmcx[idiy,4]^2
      xyri=bmcx[idix,3]*yri
      y2ri=yri^2
      
      #wccyri=cbind(colSums(wc01i*bmcx[idiy,4]),colSums(wc01i*xyri),colSums(wc01i*y2ri))
      #crlyri=cbind(rep(sum(bmcx[idiy,4]),ncol(wc01i)),colSums(wc01i*bmcx[idiy,4]),rep(sum(y2ri),ncol(wc01i)))
      #bbcyri=cbind(colSums(bc01i*bmcx[idiy,4]),colSums(bc01i*xyri),colSums(bc01i*y2ri))
      #wccyri=cbind(colSums(wc01i*yri),colSums(wc01i*xyri),colSums(wc01i*y2ri))
      wccyri=cbind(wc01i %*% yri,wc01i %*% xyri,wc01i %*% y2ri)
      crlyri=cbind(rep(sum(yri),nrow(wc01i)),wccyri[,1],rep(sum(y2ri),nrow(wc01i)))
      bbcyri=cbind(colSums(bc01i*yri),colSums(bc01i*xyri),colSums(bc01i*y2ri))
      
      crlyr=crlyr+crlyri
      wccyr=wccyr+wccyri
      bbcyr=bbcyr+bbcyri
    }
    rv<-function(mv){((mv[,6]*mv[,3])-(mv[,1]*mv[,2]))/((((mv[,6]*mv[,4])-(mv[,1]^2))^0.5)*(((mv[,6]*mv[,5])-(mv[,2]^2))^0.5))}
    r.rd=cbind(crl.rd=-rv(cbind(crlxy[,1],crlyr[,1:2],crlxy[,4],crlyr[,3],crlxy[,6])),
               wcc.rd=rv(cbind(wccxy[,1],wccyr[,1:2],wccxy[,4],wccyr[,3],wccxy[,6])),
               bcc.rd=rv(cbind(bbcxy[,1],bbcyr[,1:2],bbcxy[,4],bbcyr[,3],bbcxy[,6])))
  }
  
  t1=Sys.time()
  r.rand1=cwb.rand(1,temp.wd,dsfn,iseq,stepw,nm,nxy,break.pts,
                   permat,crlxy,wccxy,bbcxy,idcxy,y.wd,y.desc,bkptmc)
  tr1=difftime(Sys.time(),t1,units = "hours")
  
  c1 <- try(parallel::makeCluster(nworker, type = "PSOCK"))
  if (class(c1)[1] == "try-error") {c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
  if (class(c1)[1] == "try-error") {c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
  message("Now parallel computing. may take ",round(tr1*permutations/nworker,4)," hours. Begin at ", date(), ". Please wait...")
  r.rand <- parallel::parSapply(c1, 2:permutations, cwb.rand, 
                                temp.wd,dsfn,iseq,stepw,nm,nxy,break.pts,
                                permat,crlxy,wccxy,bbcxy,idcxy,y.wd,y.desc,bkptmc,
                                simplify = "array")
  parallel::stopCluster(c1)
  gc()
  r.rand=array(c(r.rand1,r.rand),dim=(dim(r.rand)+c(0,0,1)))
  
  EPS <- sqrt(.Machine$double.eps)
  pvm=r.obs;pvm[]=NA;pvad=pvm
  for(i in 1:3)
  {
    pvm[,i]=(rowSums((r.rand[,i,]-r.obs[,i])>=(-EPS),na.rm=TRUE)+1)/(rowSums(!is.na(r.rand[,i,]))+1)
    pvm[which(is.na(r.obs[,i])),i]=NA
  }
  pv1n=(rowSums((r.rand[,1,]-r.obs[,1])<=EPS,na.rm=TRUE)+1)/(rowSums(!is.na(r.rand[,1,]))+1)
  pvm[which(r.obs[,1]<0),1]=pv1n[which(r.obs[,1]<0)]
  for(i in 1:3){pvad[,i]=stats::p.adjust(pvm[,i],method = padjust.method)}
  if(cleanup){message("clean up temprary big files.");gc();file.remove(bkfn);file.remove(dsfn)}
  
  list(mantel.res=data.frame(Break.point.start=break.pts[1:(length(break.pts)-1)],Break.point.end=break.pts[2:length(break.pts)],
                             Class.midpoint=0.5*(break.pts[1:(length(break.pts)-1)]+break.pts[2:length(break.pts)]),
                             n.withinclass=wccxy[,6],n.beforeend=bbcxy[,6],stxys,
                             r.correlog=r.obs[,1],r.within.class=r.obs[,2],r.before.endpt=r.obs[,3],
                             p.correlog=pvm[,1],p.within.class=pvm[,2],p.before.endpt=pvm[,3],
                             p.adj.correlog=pvad[,1],p.adj.within.class=pvad[,2],p.adj.before.endpt=pvad[,3]),
       n.class=length(break.pts)-1,
       break.pts=break.pts,
       padjust.method=padjust.method,
       permutations=permutations,
       strata=strata,diag.in=diag.in,
       temp.wd=temp.wd)
}
