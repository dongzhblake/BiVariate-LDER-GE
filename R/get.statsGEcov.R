format_sumstats <- function(assoc,assoc2,assoc3=NULL){
    if(length(which(is.na(assoc$z)))>0){assoc <- assoc[-which(is.na(assoc$z)),]}
    if(length(which(is.na(assoc2$z)))>0){assoc2 <- assoc2[-which(is.na(assoc2$z)),]}
    assoc=assoc[assoc$snp %in% assoc2$snp,];assoc2=assoc2[assoc2$snp %in% assoc$snp,]
    assoc=assoc[order(assoc$snp),];assoc2=assoc2[order(assoc2$snp),]
    order1=ifelse(assoc$a0==assoc2$a0 & assoc$a1==assoc2$a1,T,F)
    order2=ifelse(assoc$a0==assoc2$a1 & assoc$a1==assoc2$a0,T,F)
    order3=(order1+order2)>0
    assoc2$a2=assoc2$a0
    assoc2[order2,]$a0=assoc2[order2,]$a1;assoc2[order2,]$a1=assoc2[order2,]$a2
    assoc2[order2,]$z=-assoc2[order2,]$z
    assoc=assoc[order3,];assoc2=assoc2[order3,]
    assoc2=assoc2[,-6]
    assoc1new=assoc
    assoc2new=assoc2
    if (!is.null(assoc3)){
      if(length(which(is.na(assoc3$z)))>0){assoc3 <- assoc3[-which(is.na(assoc3$z)),]}
      assoc=assoc[assoc$snp %in% assoc3$snp,];assoc3=assoc3[assoc3$snp %in% assoc$snp,]
      assoc=assoc[order(assoc$snp),];assoc3=assoc3[order(assoc3$snp),]
      order1=ifelse(assoc$a0==assoc3$a0 & assoc$a1==assoc3$a1,T,F)
      order2=ifelse(assoc$a0==assoc3$a1 & assoc$a1==assoc3$a0,T,F)
      order3=(order1+order2)>0
      assoc3$a2=assoc3$a0
      assoc3[order2,]$a0=assoc3[order2,]$a1;assoc3[order2,]$a1=assoc3[order2,]$a2
      assoc3[order2,]$z=-assoc3[order2,]$z
      assoc=assoc[order3,];assoc3=assoc3[order3,]
      assoc3=assoc3[,-6]
      return(list(assoc=assoc,assoc2=assoc2new,assoc3=assoc3))}
    return(list(assoc=assoc,assoc2=assoc2new))}

get.stats.new <- function(j,assoc,assoc2=NULL,assoc3=NULL,ldpath,ldpath.shrink,n.ld=1e4,impute_missing=F){
  if(length(which(is.na(assoc$z)))>0){assoc <- assoc[-which(is.na(assoc$z)),]}
  if(is.null(ldpath.shrink)){ldpath.shrink <- ldpath}
  if(file.exists(paste0(ldpath.shrink, '/SNPINFO', j, '.txt'))){
    ld.info <- fread(paste0(ldpath.shrink, '/SNPINFO', j, '.txt'))
    ld.info$order <- 1:nrow(ld.info)
    assoc.sub <- assoc[which(assoc$chr==ld.info$CHR[1]),]
    assoc.sub <- assoc.sub[assoc.sub$snp %in% ld.info$SNP,]
    assoc.new <- merge(assoc.sub,ld.info,by.x='snp',by.y='SNP',all.y=T)
    if(sum(!complete.cases(assoc.new))>0){
    assoc.missing=assoc.new[!complete.cases(assoc.new),]
    assoc.new[!complete.cases(assoc.new),]$chr=assoc.missing$CHR
    assoc.new[!complete.cases(assoc.new),]$a0=assoc.missing$A1
    assoc.new[!complete.cases(assoc.new),]$a1=assoc.missing$A2
    if(impute_missing==T){assoc.new[!complete.cases(assoc.new),]$z=0}
    else(assoc.new=assoc.new[complete.cases(assoc.new),])}
    if(!is.null(assoc2)){
    assoc.sub2 <- assoc2[which(assoc2$chr==ld.info$CHR[1]),]
    assoc.sub2 <- assoc.sub2[assoc.sub2$snp %in% ld.info$SNP,]
    assoc.new2 <- merge(assoc.sub2,ld.info,by.x='snp',by.y='SNP',all.y=T)
    if(sum(!complete.cases(assoc.new2))>0){
      assoc.missing2=assoc.new2[!complete.cases(assoc.new2),]
      assoc.new2[!complete.cases(assoc.new2),]$chr=assoc.missing2$CHR
      assoc.new2[!complete.cases(assoc.new2),]$a0=assoc.missing2$A1
      assoc.new2[!complete.cases(assoc.new2),]$a1=assoc.missing2$A2
      if(impute_missing==T){assoc.new2[!complete.cases(assoc.new2),]$z=0}
      else(assoc.new2=assoc.new2[complete.cases(assoc.new2),])}}
    
    if(!is.null(assoc3)){
      assoc.sub3 <- assoc3[which(assoc3$chr==ld.info$CHR[1]),]
      assoc.sub3 <- assoc.sub3[assoc.sub3$snp %in% ld.info$SNP,]
      assoc.new3 <- merge(assoc.sub3,ld.info,by.x='snp',by.y='SNP',all.y=T)
      if(sum(!complete.cases(assoc.new3))>0){
        assoc.missing3=assoc.new3[!complete.cases(assoc.new3),]
        assoc.new3[!complete.cases(assoc.new3),]$chr=assoc.missing3$CHR
        assoc.new3[!complete.cases(assoc.new3),]$a0=assoc.missing3$A1
        assoc.new3[!complete.cases(assoc.new3),]$a1=assoc.missing3$A2
        if(impute_missing==T){assoc.new3[!complete.cases(assoc.new3),]$z=0}
        else(assoc.new3=assoc.new3[complete.cases(assoc.new3),])}}
    
    if(dim(assoc.new)[1]>0){
      # do not consider ambiguous strand
      #sign <- agtc(assoc.new$a0, assoc.new$a1, assoc.new$A1, assoc.new$A2)
      if(dim(assoc.new)[1]>1){
        assoc.new <- assoc.new[order(assoc.new$BP),]
        if(!is.null(assoc2)){assoc.new2 <- assoc.new2[order(assoc.new2$BP),]}
        if(!is.null(assoc3)){assoc.new3 <- assoc.new3[order(assoc.new3$BP),]}
        if(nrow(ld.info)==1){ld.shrink <- 1 }
        else{
          ld.shrink <- as.matrix(fread(paste0(ldpath.shrink, '/LD', j, '.txt')))
          R0 <- ld.shrink[assoc.new$order,assoc.new$order]
          R0[which(is.na(R0))] <- 0 }
        if(nrow(ld.info)==nrow(assoc.new) & file.exists(paste0(ldpath.shrink,"/eigenresult_",j))){temp <- readRDS(paste0(ldpath.shrink,"/eigenresult_",j))}
        else{temp <- eigen(R0,symmetric = T)}
        temp$values[temp$values<1e-6] <- 0
        U <- temp$vectors
        V <- diag(temp$values)
        V.inv <-  1/V
        V.inv[which(V.inv==Inf)] <- 0
        eff.num <- length(temp$values)
        if(is.null(eff.num)){eff.num <- 1}
        eigen.mat<- sqrt(V.inv[1:eff.num,1:eff.num])%*%(t(U)[1:eff.num,])
        lam <- diag(V)[1:eff.num] ## export
        z1 <- assoc.new$z
        x1 <- eigen.mat%*%z1 ## export
        if(!is.null(assoc2)){
          z2 <- assoc.new2$z
          x2 <- eigen.mat%*%z2 }
        if(!is.null(assoc3)){
          z3 <- assoc.new3$z
          x3 <- eigen.mat%*%z3 }
        if(nrow(ld.info)==1){R0 <- 1}
        else{
          ld <- as.matrix(fread(paste0(ldpath, '/LD', j, '.txt')))
          R0 <- ld[assoc.new$order,assoc.new$order]
          R0[which(is.na(R0))] <- 0 }
        ldsc <- ldscore(R0,N=n.ld) ## export
        l1=list(x=x1,z=z1,lam=lam,ldsc=ldsc)
        if(!is.null(assoc2) & is.null(assoc3)){l1=list(x=x1,z=z1,x2=x2,z2=z2,lam=lam,ldsc=ldsc)}
        else if(!is.null(assoc2) & !is.null(assoc3)){l1=list(x=x1,z=z1,x2=x2,z2=z2,lam=lam,ldsc=ldsc,x3=x3,z3=z3)}
        return(l1)}
      else{print(paste0('0 SNPs in LD ',j))}}
    else{print(paste0('0 SNPs in LD ',j))}}
  else{print(paste0(ldpath.shrink, '/SNPINFO', j, '.txt'))
    print(paste0('SNPINFO',j, ' does not exist'))}}

require("utils")

extract_stats <- function(nblocks,assoc=assoc,assoc2=NULL,assoc3=NULL,ldpath=ldpath,ldpath.shrink=ldpath.shrink,n.ld=n.ld,impute_missing=F){
  stats <- list()
  pb <- txtProgressBar(0,nblocks,style=3)
  print("Matching summary statistics with LD panel")
  for(i in 0:nblocks) {l=get.stats.new(i,assoc=assoc,assoc2=assoc2,assoc3=assoc3,ldpath=ldpath,ldpath.shrink=ldpath.shrink,n.ld=n.ld,impute_missing=impute_missing)
    stats=append(stats,list(l))
  setTxtProgressBar(pb, i) }
  close(pb)
  print("LD matching done")
  return(stats)
}
