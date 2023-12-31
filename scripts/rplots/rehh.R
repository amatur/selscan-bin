library("rehh")
# make.example.files()
# 

setwd("~/workspace/selscan-bin/data/rehh/")
# hap<-data2haplohh(hap_file="example1.thap",map_file="example1.map",haplotype.in.columns=TRUE,
#                   recode.allele=T,chr.name="chr1")
# res.ihh<-scan_hh_full(hap,discard_integration_at_border = FALSE)
# ihs = ihh2ihs(res.ehh)


# ihs = ihh2ihs(res.ehh)


# takes time to recode
# hap<-data2haplohh(hap_file="o20k.thap",map_file="o20k.rehh2.map",haplotype.in.columns=TRUE,
#                   recode.allele=T,chr.name="10")
# res.ehh<-calc_ehh(hap,mrk=1)

hap<-data2haplohh(hap_file="o20k_12.thap",map_file="o20k_12.map",haplotype.in.columns=TRUE,
                  chr.name="10", min_maf=0.05)
res.ihh<-scan_hh_full(hap,discard_integration_at_border = FALSE)
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
write.table(data_cols,file="out_ihh_20k_rehh",sep=" ",row.names=FALSE, col.names =FALSE)

# 
v=(res.ehh[["ehh"]][["EHH_A"]])
my_log <- file("my_log.txt")
sink(my_log, append = TRUE, type = "output") 
print(v)
# write.matrix(v,file="Mat.csv")
# plot(v, type = "o")