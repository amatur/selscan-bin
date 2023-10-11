library("rehh")
# make.example.files()
# 
# hap<-data2haplohh(hap_file="example1.thap",map_file="example1.map",haplotype.in.columns=TRUE,
#                   recode.allele=T,chr.name="chr1")
# res.ehh<-calc_ehh(hap,mrk=1)

hap<-data2haplohh(hap_file="o20k.thap",map_file="o20k.rehh2.map",haplotype.in.columns=TRUE,
                  recode.allele=T,chr.name="10")
res.ehh<-calc_ehh(hap,mrk=1)



# 
v=(res.ehh[["ehh"]][["EHH_A"]])
my_log <- file("my_log.txt")
sink(my_log, append = TRUE, type = "output") 
print(v)
# write.matrix(v,file="Mat.csv")
# plot(v, type = "o")