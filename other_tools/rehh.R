#!/usr/bin/env Rscript


# library("optparse")
# option_list = list(
#   make_option(c("-f", "--file"), type="character", default=NULL, 
#               help="dataset file name", metavar="character"),
#   make_option(c("-o", "--out"), type="character", default="out.txt", 
#               help="output file name [default= %default]", metavar="character")
# ); 
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);




library("rehh")
# wget https://cran.r-project.org/src/contrib/data.table_1.14.8.tar.gz
# install.packages("data.table_1.14.8.tar.gz", repos = NULL, type="source")

args = commandArgs(trailingOnly=TRUE)

vcffile = "a.vcf"
NUM_THREADS=1

if (length(args)<2) {
  stop("At least two arguments must be supplied (input file).n", call.=FALSE)
} else {
  vcffile = args[1]
  NUM_THREADS = as.integer(args[2])
}

hap <- data2haplohh(vcffile, verbose = TRUE,chr.name="1", min_maf=0.05, allele_coding = "01", polarize_vcf=F, )
#recode.allele=T,

#setwd("~/workspace/selscan-bin/data/rehh/")


# hap<-data2haplohh(hap_file="example1.thap",map_file="example1.map",haplotype.in.columns=TRUE,
#                   recode.allele=T,chr.name="chr1")
# res.ihh<-scan_hh_full(hap,discard_integration_at_border = FALSE)
# ihs = ihh2ihs(res.ehh)


# ihs = ihh2ihs(res.ehh)


# takes time to recode
# hap<-data2haplohh(hap_file="o20k.thap",map_file="o20k.rehh2.map",haplotype.in.columns=TRUE,
#                   recode.allele=T,chr.name="10")
# res.ehh<-calc_ehh(hap,mrk=1)


####### IMPO
# hap<-data2haplohh(hap_file="o20k_12.thap",map_file="o20k_12.map",haplotype.in.columns=TRUE,
#                   chr.name="10", min_maf=0.05)
####### IMPO

res.ihh<-scan_hh_full(hap,discard_integration_at_border = FALSE, threads=NUM_THREADS)
#res.ehh<-calc_ehh(hap,mrk=1)
#res.ehh<-calc_ehh(hap,mrk=1,limehh=0.0, limhaplo = 2)
# res = scan_hh(hap, threads=4)
# 
# res.ihs = ihh2ihs(res)

col_a = (res.ihh[["POSITION"]])
col_b = (res.ihh[["FREQ_D"]])
col_c = (res.ihh[["IHH_D"]])
col_d = (res.ihh[["IHH_A"]])
col_e = log10(col_c/col_d)

data_cols = data.frame(col_a, col_b, col_c, col_d, col_e)
write.table(data_cols,file=paste(vcffile,".rehh",sep=""),sep=" ",row.names=FALSE, col.names =FALSE)

# 
# v=(res.ehh[["ehh"]][["EHH_A"]])
# my_log <- file("my_log.txt")
# sink(my_log, append = TRUE, type = "output") 
# print(v)
