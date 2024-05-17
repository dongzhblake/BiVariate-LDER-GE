#' @title Main function
#' @description Run LDER-GE 
#' @param assoc GWIS summary statistics, need to include snp, chr, a0, a1, z (header is necessary)
#' @param n.gwas The sample size of the GWAS summary statistics
#' @param path The path of LD panel directory
#' @param LD.insample T/F, whether the LD reference is estimated with the GWAS cohort (T) or external reference panel (e.g. 1000 Genome Project) (F)
#' @param n.ld The sample size of the LD reference
#' @param cores The number of cores for computation in parallel
#' @param method 'lder', 'ldsc', or 'both'
#' @param size_num Number of blocks for jackknife
#' @param type if 'jack' then conducts delete-wise block jackknife. if 'none' then does not do inference 
#' @import  data.table stats utils
#' @export
#'
#'
runBV_LDER_GE <- function(assoc_gwis, assoc_gwas, n.gwis, n.gwas, n.overlap, path, LD.insample=T,  n.ld,method='lder', type='jack',size_num=200){
  library(parallel)
  refo_sts=format_sumstats(assoc_gwis,assoc_gwas)
  assoc1=refo_sts$assoc;assoc2=refo_sts$assoc2
  A=unlist(strsplit(list.files(path,pattern="INFO"),".txt"));B=unlist(strsplit(A,"SNPINFO"));C=max(as.numeric(B),na.rm=T)
  ldpath <- ldpath.shrink <- path
  print("matching summary statistics with LD panel")
  bivariate_stats <- extract_stats(C,assoc=assoc1,assoc2=assoc2,assoc3=NULL,ldpath=ldpath,ldpath.shrink=ldpath.shrink,n.ld=n.ld,impute_missing=F)
  N1=n.gwis
  N2=n.gwas
  Ns=ifelse(is.null(n.overlap),min(n.gwis,n.gwas),n.overlap)
  if(method=='lder'){
    res <- bivariate_lder(bivariate_stats,N1,N2,NS=min(N1,N2),size_num=size_num)
    bitest=test_bivariate_GE(res,"BVN")
    return(list(lder = res, BV_test_lder = bitest))
  }else if(method=='ldsc'){
    res <- bivariate_ldsc(bivariate_stats,N1,N2,NS=min(N1,N2),size_num=size_num)
    bitest=test_bivariate_GE(res,"BVN")
    return(list(ldsc = res, BV_test_ldsc = bitest))
  }else if(method=='both'){
    res1 <- bivariate_lder(bivariate_stats,N1,N2,NS=min(N1,N2),size_num=size_num)
    res2 <- bivariate_ldsc(bivariate_stats,N1,N2,NS=min(N1,N2),size_num=size_num)
    bitest1=test_bivariate_GE(res1,"BVN")
    bitest2=test_bivariate_GE(res2,"BVN")
    return(list(lder=res1,ldsc=res2, BV_test_lder = bitest1, BV_test_ldsc = bitest2))
  }
}

