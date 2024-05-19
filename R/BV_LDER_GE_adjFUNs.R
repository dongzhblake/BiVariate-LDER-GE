#' @title Main function
#' @description Run BVLDER-GE-adjust
#' @param assoc_gwis_Y GWIS (GE interaction effect) summary statistics for Y, need to include snp, chr, a0, a1, z (header is necessary)
#' @param assoc_gwas_Y GWAS (additive genetic effect) summary statistics for Y, need to include snp, chr, a0, a1, z (header is necessary)
#' @param assoc_gwas_E GWAS (additive genetic effect) summary statistics for E, need to include snp, chr, a0, a1, z (header is necessary)
#' @param n.gwis_Y The sample size of the GWIS (GE interaction effect) summary statistics for Y
#' @param n.gwas_Y The sample size of the GWAS (additive genetic effect) summary statistics for Y
#' @param n.gwas_E The sample size of the GWAS (additive genetic effect) summary statistics for E
#' @param path The path of LD panel directory
#' @param LD.insample T/F, whether the LD reference is estimated with target cohort (T) or external reference panel (e.g. 1000 Genome Project) (F)
#' @param n.ld The sample size of the LD reference
#' @param method 'lder', 'ldsc', or 'both'
#' @param size_num Number of blocks for jackknife
#' @param R2 R2 between E and Y. If focusing on conditional heritability (i.e., Y adjusted for E already) then set to 0.
#' @import  data.table stats utils
#' @export
#'
#'
runBV_LDER_GE_adj <- function(assoc_gwis_Y, assoc_gwas_Y, assoc_gwas_E, n.gwis_Y, n.gwas_Y, n.gwas_E, path, LD.insample=T,  n.ld,method='lder',size_num=200,R2=0){
  library(parallel)
  refo_sts=format_sumstats(assoc_gwis_Y, assoc_gwas_Y, assoc_gwas_E)
  assoc1=refo_sts$assoc;assoc2=refo_sts$assoc2;assoc3=refo_sts$assoc3
  A=unlist(strsplit(list.files(path,pattern="INFO"),".txt"));B=unlist(strsplit(A,"SNPINFO"));C=max(as.numeric(B),na.rm=T)
  ldpath <- ldpath.shrink <- path
  print("matching summary statistics with LD panel")
  trivariate_stats <- extract_stats(C,assoc=assoc1,assoc2=assoc2,assoc3=assoc3,ldpath=ldpath,ldpath.shrink=ldpath.shrink,n.ld=n.ld,impute_missing=F)
  N1 = n.gwis_Y
  N2 = n.gwas_Y
  N3 = n.gwas_E
  if(method=='lder'){
    res <- bivariate_lder_correct(trivariate_stats,N1,N2,N3,size_num=200,R2)
    bitest=test_bivariate_GE(res,"BVN")
    return(list(lder = res, BV_test_lder = bitest))
  }else if(method=='ldsc'){
    res <- bivariate_ldsc_correct(trivariate_stats,N1,N2,N3,size_num=200,R2)
    bitest=test_bivariate_GE(res,"BVN")
    return(list(ldsc = res, BV_test_ldsc = bitest))
  }else if(method=='both'){
    res1 <- bivariate_lder_correct(trivariate_stats,N1,N2,N3,size_num=200,R2)
    res2 <- bivariate_ldsc_correct(trivariate_stats,N1,N2,N3,size_num=200,R2)
    bitest1=test_bivariate_GE(res1,"BVN")
    bitest2=test_bivariate_GE(res2,"BVN")
    return(list(lder=res1,ldsc=res2, BV_test_lder = bitest1, BV_test_ldsc = bitest2))
  }
}


