calRg.new<-function(x1, x2, ldsc, N1, N2, Ns=min(N1, N2)){
  z1=x1;z2=x2
  if(length(z1)!=length(z2)){ stop('ldsr.calRg:lengths of z1 and z2 are not matched!')}
  if(length(z1)!=length(ldsc)){ stop('ldsr.calRg:lengths of z and ldsc are not matched!')}
  updateFunc_bivariate1 <- function(coeff){
    temp1 <- (coeff[1]>sqrt(N1*N2)/M)*sqrt(N1*N2)/M + (coeff[1]<(-sqrt(N1*N2)/M))*-sqrt(N1*N2)/M + ((coeff[1]<=sqrt(N1*N2)/M & coeff[1]>=-sqrt(N1*N2)/M)*coeff[1])
    w1 <- 1/((N1*h2est1*ldsc/M+1)*(N2*h2est2*ldsc/M+1)+(temp1*ldsc+coeff[2])^2)
    w2 <- pmin(1, ldsc)
    return(w1*w2)}
  M <- length(z1)
  y <- z1*z2
  h2estObj1<-calH2.new1(z1, ldsc, N1, NULL)
  h2estObj2<-calH2.new1(z2, ldsc, N2, NULL)
  h2est1=h2est11 <- as.numeric(h2estObj1$h2)
  h2est1=ifelse(h2est1>0,h2est1,0)
  h2est2=h2est22 <- as.numeric(h2estObj2$h2)
  h2est2=ifelse(h2est2>0,h2est2,0)
  x <- cbind(ldsc, rep(1, length(ldsc)))
  coeff <- irwls(x, y, updateFunc_bivariate1)
  rhoNs <- coeff[[2]]*sqrt(N1*N2)
  rhog <- coeff[[1]]*M/sqrt(N1*N2)
  rg<- rhog/sqrt(h2est1*h2est2)
  return(list(rhog=rhog, rg=rg, rho=rhoNs/Ns, h2_1=h2est11, h2_2=h2est22))
}

calRg <- function(z1, z2, ldsc, N1, N2, Ns=min(N1, N2)){
  if(length(z1)!=length(z2)){ stop('ldsr.calRg:lengths of z1 and z2 are not matched!')}
  if(length(z1)!=length(ldsc)){ stop('ldsr.calRg:lengths of z and ldsc are not matched!')}
  updateFunc_bivariate <- function(coeff){
    temp1 <- (coeff[1]>sqrt(N1*N2)/M)*sqrt(N1*N2)/M + (coeff[1]<(-sqrt(N1*N2)/M))*-sqrt(N1*N2)/M + ((coeff[1]<=sqrt(N1*N2)/M & coeff[1]>=-sqrt(N1*N2)/M)*coeff[1])
    w1 <- 1/((N1*h2est1*ldsc/M+1)*(N2*h2est2*ldsc/M+1)+(temp1*ldsc+coeff[2])^2)
    w2 <- 1/pmax(ldsc, 1)
    return(w1*w2)}
  M <- length(z1)
  y <- z1*z2
  h2estObj1<-calH2(z1, ldsc, N1, NULL)
  h2estObj2<-calH2(z2, ldsc, N2, NULL)
  h2est1=h2est11 <- as.numeric(h2estObj1$h2)
  h2est1=ifelse(h2est1>0,h2est1,0)
  h2est2=h2est22 <- as.numeric(h2estObj2$h2)
  h2est2=ifelse(h2est2>0,h2est2,0)
  x <- cbind(ldsc, rep(1, length(ldsc)))
  coeff <- irwls(x, y, updateFunc_bivariate)
  rhoNs <- coeff[[2]]*sqrt(N1*N2)
  rhog <- coeff[[1]]*M/sqrt(N1*N2)
  rg<- rhog/sqrt(h2est1*h2est2)
  return(list(rhog=rhog, rg=rg, rho=rhoNs/Ns, h2_1=h2est11, h2_2=h2est22))
}


define_block_unit <- function(size_num=200,mj,sizedf){
  sizes_perb=c()
  for(b in seq(0.5,2.5,by=0.01)){
    target_size=sum(mj)/size_num/b
    grouped_block=list()
    size_each=c()
    size_now=0
    block_group=c()
    for(i in 1:nrow(sizedf)){
      if(size_now<target_size){
        size_now=size_now+sizedf$mj[i]
        block_group=c(block_group,sizedf$V2[i])}
      else{
        grouped_block=append(grouped_block,list(block_group))
        size_each=c(size_each,size_now)
        size_now=sizedf$mj[i]
        block_group=c(sizedf$V2[i])}}
    sizes_perb=c(sizes_perb,length(grouped_block))}
  b=seq(0.5,2.5,by=0.01)[which(abs(sizes_perb-size_num)==min(abs(sizes_perb-size_num)))][1]
  target_size=sum(mj)/size_num/b
  grouped_block=list()
  size_each=c()
  size_now=0
  block_group=c()
  for(i in 1:nrow(sizedf)){
    if(size_now<target_size){
      size_now=size_now+sizedf$mj[i]
      block_group=c(block_group,sizedf$V2[i])}
    else{
      grouped_block=append(grouped_block,list(block_group))
      size_each=c(size_each,size_now)
      size_now=sizedf$mj[i]
      block_group=c(sizedf$V2[i])}}
  sizes_perb=c(sizes_perb,length(grouped_block))
  return(return(grouped_block))}

bivariate_do.jack <- function(bivariate_stats,j,N1,N2,NS=min(N1,N2)){
  bivariate_stats <- bivariate_stats[-j]
  x1 <- unlist(sapply(bivariate_stats,'[[',1))
  lam <- unlist(sapply(bivariate_stats,'[[',5))
  x2 <- unlist(sapply(bivariate_stats,'[[',3))
  res <- calRg.new(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  res$rhog <- res$rhog/length(x1)
  res$h2_1 <- res$h2_1/length(x1)
  res$h2_2 <- res$h2_2/length(x1)
  return(res)}

bivariate_do.jack.ldsc <- function(bivariate_stats,j,N1,N2,NS=min(N1,N2)){
  bivariate_stats <- bivariate_stats[-j]
  x1 <- unlist(sapply(bivariate_stats,'[[',2))
  lam <- unlist(sapply(bivariate_stats,'[[',6))
  x2 <- unlist(sapply(bivariate_stats,'[[',4))
  res <- calRg(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  res$rhog <- res$rhog/length(x1)
  res$h2_1 <- res$h2_1/length(x1)
  res$h2_2 <- res$h2_2/length(x1)
  return(res)}


bivariate_lder <- function(bivariate_stats,N1,N2,NS=min(N1,N2),size_num=200,R2=0){
  C=sqrt(1-R2)
  temp1 <- unlist(lapply(bivariate_stats,length))
  bivariate_stats[temp1==1] <- NULL
  x1 <- unlist(sapply(bivariate_stats,'[[',1))
  lam <- unlist(sapply(bivariate_stats,'[[',5))
  x2 <- unlist(sapply(bivariate_stats,'[[',3))
  res <- calRg.new(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  mj=sapply(1:length(bivariate_stats),FUN=function(x)(length(bivariate_stats[[x]]$ldsc)))
  sizedf=as.data.frame(cbind(mj,1:length(bivariate_stats)))
  grouped_block=define_block_unit(size_num=size_num,mj=mj,sizedf=sizedf)
    pb <- txtProgressBar(0,length(grouped_block),style=3)
    temp1.jack=list()
    print("Estimating SE with delete-block-jackknife")
    for(t in 1:length(grouped_block)){temp1.jack=append(temp1.jack,list(bivariate_do.jack(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N1=N1,N2=N2,NS=min(N1,N2))))
    setTxtProgressBar(pb, t)}
    print("delete-block-jackknife done")
    close(pb)
    rrhog <-  unlist(sapply(temp1.jack,'[[','rhog'))
    rrho <- unlist(sapply(temp1.jack,'[[','rho'))
    h2_1 <- unlist(sapply(temp1.jack,'[[','h2_1'))
    h2_2 <- unlist(sapply(temp1.jack,'[[','h2_2'))
    n.block <- length(grouped_block)
    rrhog33.jack <- sd(rrhog)*sqrt(n.block)*length(x1)*C
    rrho33.jack <- sd(rrho)*sqrt(n.block)*C
    h2_133.jack <- sd(h2_1)*sqrt(n.block)*length(x1)*C^2
    h2_233.jack <- sd(h2_2)*sqrt(n.block)*length(x1)
    cov_rrhog_h2_1 <- cov(rrhog,h2_1)*length(x1)^2*n.block*C^3
    return(list(genecov=res$rhog*C,rho=res$rho*C,h2I=res$h2_1*C^2,h2g=res$h2_2,genecov.se=rrhog33.jack,rho.se=rrho33.jack,h2I.se=h2_133.jack,h2g.se=h2_233.jack,cov_rrhog_h2_1=cov_rrhog_h2_1,genecov.p=pchisq((res$rhog/rrhog33.jack)^2,df=1,lower.tail = F),h2I.p=pchisq((res$h2_1/h2_133.jack)^2,df=1,lower.tail = F),h2g.p=pchisq((res$h2_2/h2_233.jack)^2,df=1,lower.tail = F)))
  }



bivariate_ldsc <- function(bivariate_stats,N1,N2,NS=min(N1,N2),size_num=200,R2=0){
  C=sqrt(1-R2)
  temp1 <- unlist(lapply(bivariate_stats,length))
  bivariate_stats[temp1==1] <- NULL
  x1 <- unlist(sapply(bivariate_stats,'[[',2))
  lam <- unlist(sapply(bivariate_stats,'[[',6))
  x2 <- unlist(sapply(bivariate_stats,'[[',4))
  res <- calRg(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  mj=sapply(1:length(bivariate_stats),FUN=function(x)(length(bivariate_stats[[x]]$ldsc)))
  sizedf=as.data.frame(cbind(mj,1:length(bivariate_stats)))
  grouped_block=define_block_unit(size_num=size_num,mj=mj,sizedf=sizedf)
    pb <- txtProgressBar(0,length(grouped_block),style=3)
    temp1.jack=list()
    print("Estimating SE with delete-block-jackknife")
    for(t in 1:length(grouped_block)){temp1.jack=append(temp1.jack,list(bivariate_do.jack.ldsc(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N1=N1,N2=N2,NS=min(N1,N2))))
    setTxtProgressBar(pb, t)}
    print("delete-block-jackknife done")
    close(pb)
    rrhog <-  unlist(sapply(temp1.jack,'[[','rhog'))
    rrho <- unlist(sapply(temp1.jack,'[[','rho'))
    h2_1 <- unlist(sapply(temp1.jack,'[[','h2_1'))
    h2_2 <- unlist(sapply(temp1.jack,'[[','h2_2'))
    n.block <- length(grouped_block)
    rrhog33.jack <- sd(rrhog)*sqrt(n.block)*length(x1)*C
    rrho33.jack <- sd(rrho)*sqrt(n.block)*C
    h2_133.jack <- sd(h2_1)*sqrt(n.block)*length(x1)*C^2
    h2_233.jack <- sd(h2_2)*sqrt(n.block)*length(x1)
    cov_rrhog_h2_1 <- cov(rrhog,h2_1)*length(x1)^2*n.block*C^3
    return(list(genecov=res$rhog*C,rho=res$rho*C,h2I=res$h2_1*C^2,h2g=res$h2_2,genecov.se=rrhog33.jack,rho.se=rrho33.jack,h2I.se=h2_133.jack,h2g.se=h2_233.jack,cov_rrhog_h2_1=cov_rrhog_h2_1,genecov.p=pchisq((res$rhog/rrhog33.jack)^2,df=1,lower.tail = F),h2I.p=pchisq((res$h2_1/h2_133.jack)^2,df=1,lower.tail = F),h2g.p=pchisq((res$h2_2/h2_233.jack)^2,df=1,lower.tail = F)))
}



test_bivariate_GE <- function(bivariate_result,options="BVN"){
  if(! options %in% c("Fisher","BVN","All")){
    stop('Please specify the correct test options: "Fisher" or "BVN" or "All"') }
  c=bivariate_result$cov_rrhog_h2_1
  a11=bivariate_result$h2_1.sd^2
  a22=bivariate_result$genecov.sd^2
  estimates=c(bivariate_result$h2_1,bivariate_result$genecov)
  S=matrix(c(a11,c,c,a22),2,2)
  EI=eigen(S)
  S.5=EI$vectors %*% diag(sqrt(EI$values)) %*% t(EI$vectors)
  decor_es= as.numeric(solve(S.5) %*% estimates)
  ps=(1-pnorm(abs(decor_es)))*2
  chis=-2*sum(log(ps))
  fisher_P=pchisq(chis,df=4,lower.tail = F)
  MHN_P=exp(-sum(decor_es^2)/2)
  a1=c("Fisher","BVN")
  a2=c(fisher_P,MHN_P)
  a=as.data.frame(cbind(a1,a2))
  if(options=="All"){return(as.numeric(a$a2))}
  else{return(as.numeric(a[a$a1==options,2]))}
}

bivariate_do.jackIe <- function(bivariate_stats,j,N1,N2,NS=min(N1,N2)){
  bivariate_stats <- bivariate_stats[-j]
  x1 <- unlist(sapply(bivariate_stats,'[[',1))
  lam <- unlist(sapply(bivariate_stats,'[[',5))
  x2 <- unlist(sapply(bivariate_stats,'[[',7))
  res <- calRg.new(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  res$rhog <- res$rhog/length(x1)
  res$h2_1 <- res$h2_1/length(x1)
  res$h2_2 <- res$h2_2/length(x1)
  res$rhog <- res$rhog/(1+res$h2_2)
  return(res)}

bivariate_do.jackIe1 <- function(bivariate_stats,j,N1,N2,NS=min(N1,N2)){
  bivariate_stats <- bivariate_stats[-j]
  x1 <- unlist(sapply(bivariate_stats,'[[',2))
  lam <- unlist(sapply(bivariate_stats,'[[',6))
  x2 <- unlist(sapply(bivariate_stats,'[[',8))
  res <- calRg(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  res$rhog <- res$rhog/length(x1)
  res$h2_1 <- res$h2_1/length(x1)
  res$h2_2 <- res$h2_2/length(x1)
  res$rhog <- res$rhog/(1+res$h2_2)
  return(res)}

bivariate_do.jackge <- function(bivariate_stats,j,N1,N2,NS=min(N1,N2)){
  bivariate_stats <- bivariate_stats[-j]
  x1 <- unlist(sapply(bivariate_stats,'[[',3))
  lam <- unlist(sapply(bivariate_stats,'[[',5))
  x2 <- unlist(sapply(bivariate_stats,'[[',7))
  res <- calRg.new(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  res$rhog <- res$rhog/length(x1)
  res$h2_1 <- res$h2_1/length(x1)
  res$h2_2 <- res$h2_2/length(x1)
  return(res)}

bivariate_do.jackge1 <- function(bivariate_stats,j,N1,N2,NS=min(N1,N2)){
  bivariate_stats <- bivariate_stats[-j]
  x1 <- unlist(sapply(bivariate_stats,'[[',2))
  lam <- unlist(sapply(bivariate_stats,'[[',6))
  x2 <- unlist(sapply(bivariate_stats,'[[',8))
  res <- calRg(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  res$rhog <- res$rhog/length(x1)
  res$h2_1 <- res$h2_1/length(x1)
  res$h2_2 <- res$h2_2/length(x1)
  return(res)}

lderge_correct <- function(bivariate_stats,N1,N2,size_num=200,R2=0){
  # the bivariate_stats includes the GWIS for Y and GWAS for E
  C=sqrt(1-R2)
  temp1 <- unlist(lapply(bivariate_stats,length))
  bivariate_stats[temp1==1] <- NULL
  x1 <- unlist(sapply(bivariate_stats,'[[',1))
  lam <- unlist(sapply(bivariate_stats,'[[',5))
  x2 <- unlist(sapply(bivariate_stats,'[[',3))
  res2 <- calRg.new(x1, x2, lam, N1, N3, Ns=min(N1, N3))
  res2$rhog=res2$rhog/(1+res2$h2_2)

  mj=sapply(1:length(bivariate_stats),FUN=function(x)(length(bivariate_stats[[x]]$ldsc)))
  sizedf=as.data.frame(cbind(mj,1:length(bivariate_stats)))
  grouped_block=define_block_unit(size_num=size_num,mj=mj,sizedf=sizedf)
  n.block <- length(grouped_block)
  pb <- txtProgressBar(0,length(grouped_block),style=3)
  temp2.jack=list()
  print("Estimating SE with delete-block-jackknife")
  for(t in 1:length(grouped_block)){
    temp2.jack=append(temp2.jack,list(bivariate_do.jackIe(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N1,N3,NS=min(N1,N3))))
    setTxtProgressBar(pb, t)}
  print("delete-block-jackknife done")
  close(pb)
  h2_1 <- unlist(sapply(temp2.jack,'[[','h2_1'))
  rrhog_Ie <-  unlist(sapply(temp2.jack,'[[','rhog'))

  h2_correct_blocks = h2_1 - 2*rrhog_Ie^2
  h2_133.jack <- sd(h2_correct_blocks)*sqrt(n.block)*length(x1)*C^2
  h2_1_corrected = (res$h2_1-2*res2$rhog^2)*C^2

  return(list(h2I=h2_1_corrected,h2I.se=h2_133.jack,h2I.p=pchisq((h2_1_corrected/h2_133.jack)^2,df=1,lower.tail = F)))
}


ldscge_correct <- function(bivariate_stats,N1,N2,size_num=200,R2=0){
  # the bivariate_stats includes the GWIS for Y and GWAS for E
  C=sqrt(1-R2)
  temp1 <- unlist(lapply(bivariate_stats,length))
  bivariate_stats[temp1==1] <- NULL
  x1 <- unlist(sapply(bivariate_stats,'[[',2))
  lam <- unlist(sapply(bivariate_stats,'[[',6))
  x2 <- unlist(sapply(bivariate_stats,'[[',4))
  res2 <- calRg(x1, x2, lam, N1, N3, Ns=min(N1, N3))
  res2$rhog=res2$rhog/(1+res2$h2_2)
  
  mj=sapply(1:length(bivariate_stats),FUN=function(x)(length(bivariate_stats[[x]]$ldsc)))
  sizedf=as.data.frame(cbind(mj,1:length(bivariate_stats)))
  grouped_block=define_block_unit(size_num=size_num,mj=mj,sizedf=sizedf)
  n.block <- length(grouped_block)
  pb <- txtProgressBar(0,length(grouped_block),style=3)
  temp2.jack=list()
  print("Estimating SE with delete-block-jackknife")
  for(t in 1:length(grouped_block)){
    temp2.jack=append(temp2.jack,list(bivariate_do.jackIe1(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N1,N2,NS=min(N1,N2))))
    setTxtProgressBar(pb, t)}
  print("delete-block-jackknife done")
  close(pb)
  h2_1 <- unlist(sapply(temp2.jack,'[[','h2_1'))
  rrhog_Ie <-  unlist(sapply(temp2.jack,'[[','rhog'))
  
  h2_correct_blocks = h2_1 - 2*rrhog_Ie^2
  h2_133.jack <- sd(h2_correct_blocks)*sqrt(n.block)*length(x1)*C^2
  h2_1_corrected = (res$h2_1-2*res2$rhog^2)*C^2
  
  return(list(h2I=h2_1_corrected,h2I.se=h2_133.jack,h2I.p=pchisq((h2_1_corrected/h2_133.jack)^2,df=1,lower.tail = F)))
}  


bivariate_lder_correct <- function(trivariate_stats,N1,N2,N3,size_num=200,R2=0,print_GIE=F){
  C=sqrt(1-R2)
  bivariate_stats=trivariate_stats
  temp1 <- unlist(lapply(bivariate_stats,length))
  bivariate_stats[temp1==1] <- NULL
  x1 <- unlist(sapply(bivariate_stats,'[[',1))
  lam <- unlist(sapply(bivariate_stats,'[[',5))
  x2 <- unlist(sapply(bivariate_stats,'[[',3))
  res <- calRg.new(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  
  x1 <- unlist(sapply(bivariate_stats,'[[',1))
  lam <- unlist(sapply(bivariate_stats,'[[',5))
  x2 <- unlist(sapply(bivariate_stats,'[[',7))
  res2 <- calRg.new(x1, x2, lam, N1, N3, Ns=min(N1, N3))
  res2$rhog=res2$rhog/(1+res2$h2_2)
  
  x1 <- unlist(sapply(bivariate_stats,'[[',3))
  lam <- unlist(sapply(bivariate_stats,'[[',5))
  x2 <- unlist(sapply(bivariate_stats,'[[',7))
  res3 <- calRg.new(x1, x2, lam, N2, N3, Ns=min(N2, N3))
  
  mj=sapply(1:length(bivariate_stats),FUN=function(x)(length(bivariate_stats[[x]]$ldsc)))
  sizedf=as.data.frame(cbind(mj,1:length(bivariate_stats)))
  grouped_block=define_block_unit(size_num=size_num,mj=mj,sizedf=sizedf)
  n.block <- length(grouped_block)
  pb <- txtProgressBar(0,length(grouped_block),style=3)
  temp1.jack=temp2.jack=temp3.jack=list()
  print("Estimating SE with delete-block-jackknife")
  for(t in 1:length(grouped_block)){
    temp1.jack=append(temp1.jack,list(bivariate_do.jack(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N1,N2,NS=min(N1,N2))))
    temp2.jack=append(temp2.jack,list(bivariate_do.jackIe(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N1,N3,NS=min(N1,N3))))
    temp3.jack=append(temp3.jack,list(bivariate_do.jackge(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N2,N3,NS=min(N2,N3))))
    setTxtProgressBar(pb, t)}
  print("delete-block-jackknife done")
  close(pb)
  rrhog <-  unlist(sapply(temp1.jack,'[[','rhog'))
  h2_1 <- unlist(sapply(temp1.jack,'[[','h2_1'))

  rrhog_Ie <-  unlist(sapply(temp2.jack,'[[','rhog'))
  rrhog_ge <-  unlist(sapply(temp3.jack,'[[','rhog'))
  
  h2_correct_blocks = h2_1 - 2*rrhog_Ie^2
  rho_Ig_correct_blocks = rrhog - rrhog_Ie*rrhog_ge
  
  rrhog33.jack <- sd(rho_Ig_correct_blocks)*sqrt(n.block)*length(x1)*C
  h2_133.jack <- sd(h2_correct_blocks)*sqrt(n.block)*length(x1)*C^2
  cov_rrhog_h2_1 <- cov(rho_Ig_correct_blocks,h2_correct_blocks)*length(x1)^2*n.block*C^3
  rrhogGIE.jack <- sd(rrhog_Ie)*sqrt(n.block)*length(x1)*C
  
  
  h2_1_corrected = (res$h2_1-2*res2$rhog^2)*C^2
  rho_Ig_corrected = res$rhog*C - res2$rhog*res3$rhog*C 
  
  if(print_GIE==T){return(list(genecov=rho_Ig_corrected,h2I=h2_1_corrected,genecov.se=rrhog33.jack,h2I.se=h2_133.jack,cov_rrhog_h2_1=cov_rrhog_h2_1,genecov.p=pchisq((rho_Ig_corrected/rrhog33.jack)^2,df=1,lower.tail = F),h2I.p=pchisq((h2_1_corrected/h2_133.jack)^2,df=1,lower.tail = F),gcov_IE=res2$rhog*C,gcov_IE.se=rrhogGIE.jack,gcov_IE.p=pchisq((res2$rhog*C/rrhogGIE.jack)^2,df=1,lower.tail = F)))}
  return(list(genecov=rho_Ig_corrected,h2I=h2_1_corrected,genecov.se=rrhog33.jack,h2I.se=h2_133.jack,cov_rrhog_h2_1=cov_rrhog_h2_1,genecov.p=pchisq((rho_Ig_corrected/rrhog33.jack)^2,df=1,lower.tail = F),h2I.p=pchisq((h2_1_corrected/h2_133.jack)^2,df=1,lower.tail = F)))
}


bivariate_ldsc_correct <- function(trivariate_stats,N1,N2,N3,size_num=200,R2=0,print_GIE=F){
  C=sqrt(1-R2)
  bivariate_stats=trivariate_stats
  temp1 <- unlist(lapply(bivariate_stats,length))
  bivariate_stats[temp1==1] <- NULL
  x1 <- unlist(sapply(bivariate_stats,'[[',2))
  lam <- unlist(sapply(bivariate_stats,'[[',6))
  x2 <- unlist(sapply(bivariate_stats,'[[',4))
  res <- calRg(x1, x2, lam, N1, N2, Ns=min(N1, N2))
  
  x1 <- unlist(sapply(bivariate_stats,'[[',2))
  lam <- unlist(sapply(bivariate_stats,'[[',6))
  x2 <- unlist(sapply(bivariate_stats,'[[',8))
  res2 <- calRg(x1, x2, lam, N1, N3, Ns=min(N1, N3))
  res2$rhog=res2$rhog/(1+res2$h2_2)
  
  x1 <- unlist(sapply(bivariate_stats,'[[',4))
  lam <- unlist(sapply(bivariate_stats,'[[',6))
  x2 <- unlist(sapply(bivariate_stats,'[[',8))
  res3 <- calRg(x1, x2, lam, N1, N3, Ns=min(N2, N3))
  
  mj=sapply(1:length(bivariate_stats),FUN=function(x)(length(bivariate_stats[[x]]$ldsc)))
  sizedf=as.data.frame(cbind(mj,1:length(bivariate_stats)))
  grouped_block=define_block_unit(size_num=size_num,mj=mj,sizedf=sizedf)
  n.block <- length(grouped_block)
  pb <- txtProgressBar(0,length(grouped_block),style=3)
  temp1.jack=temp2.jack=temp3.jack=list()
  print("Estimating SE with delete-block-jackknife")
  for(t in 1:length(grouped_block)){
    temp1.jack=append(temp1.jack,list(bivariate_do.jack.ldsc(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N1,N2,NS=min(N1,N2))))
    temp2.jack=append(temp2.jack,list(bivariate_do.jackIe1(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N1,N3,NS=min(N1,N3))))
    temp3.jack=append(temp3.jack,list(bivariate_do.jackge1(bivariate_stats=bivariate_stats,j=grouped_block[[t]],N2,N3,NS=min(N2,N3))))
    setTxtProgressBar(pb, t)}
  print("delete-block-jackknife done")
  close(pb)
  rrhog <-  unlist(sapply(temp1.jack,'[[','rhog'))
  h2_1 <- unlist(sapply(temp1.jack,'[[','h2_1'))

  rrhog_Ie <-  unlist(sapply(temp2.jack,'[[','rhog'))
  rrhog_ge <-  unlist(sapply(temp3.jack,'[[','rhog'))
  
  h2_correct_blocks = h2_1 - 2*rrhog_Ie^2
  rho_Ig_correct_blocks = rrhog - rrhog_Ie*rrhog_ge
  
  rrhog33.jack <- sd(rho_Ig_correct_blocks)*sqrt(n.block)*length(x1)*C
  h2_133.jack <- sd(h2_correct_blocks)*sqrt(n.block)*length(x1)*C^2
  cov_rrhog_h2_1 <- cov(rho_Ig_correct_blocks,h2_correct_blocks)*length(x1)^2*n.block*C^3
  rrhogGIE.jack <- sd(rrhog_Ie)*sqrt(n.block)*length(x1)*C
  
  
  h2_1_corrected = (res$h2_1-2*res2$rhog^2)*C^2
  rho_Ig_corrected = res$rhog*C - res2$rhog*res3$rhog*C 

  if(print_GIE==T){return(list(genecov=rho_Ig_corrected,h2I=h2_1_corrected,genecov.se=rrhog33.jack,h2I.se=h2_133.jack,cov_rrhog_h2_1=cov_rrhog_h2_1,genecov.p=pchisq((rho_Ig_corrected/rrhog33.jack)^2,df=1,lower.tail = F),h2I.p=pchisq((h2_1_corrected/h2_133.jack)^2,df=1,lower.tail = F),gcov_IE=res2$rhog*C,gcov_IE.se=rrhogGIE.jack,gcov_IE.p=pchisq((res2$rhog*C/rrhogGIE.jack)^2,df=1,lower.tail = F)))}
  return(list(genecov=rho_Ig_corrected,h2I=h2_1_corrected,genecov.se=rrhog33.jack,h2I.se=h2_133.jack,cov_rrhog_h2_1=cov_rrhog_h2_1,genecov.p=pchisq((rho_Ig_corrected/rrhog33.jack)^2,df=1,lower.tail = F),h2I.p=pchisq((h2_1_corrected/h2_133.jack)^2,df=1,lower.tail = F)))
}






