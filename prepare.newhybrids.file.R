vcf<-read.vcfR("/Users/devder/Desktop/parrots/la.parrots/parrotsnpsunfiltered.vcf")
#read in locality info for samples
locs<-read.table("~/Desktop/parrots/la.parrots/parrotsamples.txt", header=T, sep="\t")
colnames(locs)<-c("id","location","tiss")
head(vcf@fix)
vcf@gt[1:10,1:10]
#check that sampling IDs in text file match sample IDs in the vcf
colnames(vcf@gt)[-1] == locs$id

#remove SNPs less than 75% complete
vcf<-missing.by.snp(vcf, cutoff = .75)
#probably worth doing some cursory filtering for missing data

#extract genotype matrix
mat<-extract.gt(vcf)
#remove loci that aren't biallelic
mat<-mat[nchar(vcf@fix[,5]) == 1,]
table(mat)
dim(mat) #4658 x 30

#convert genotypes into newhybrids format
mat<-t(mat)
mat[mat == "0/0"]<-"00"
mat[mat == "0/1"]<-"01"
mat[mat == "1/1"]<-"11"
mat[is.na(mat)]<-"0"
table(mat)
#mat<-as.data.frame(mat)
mat[1:5,1:5]
inds<-rownames(mat)
loci<-c("LocusNames",colnames(mat))
rownames(mat)<-NULL
colnames(mat)<-NULL
id<-c("","",rep("z1s", times=4),rep("",times=20),rep("z0s", times=4))

#write file
setwd("~/Downloads")
cat("NumIndivs 30\n",file="parrots.newhybs.txt")
cat("NumLoci 4658\n",file="parrots.newhybs.txt",append=TRUE)
cat("Digits 1\n",file="parrots.newhybs.txt",append=TRUE)
cat("Format Lumped\n",file="parrots.newhybs.txt",append=TRUE)
cat("\n",file="parrots.newhybs.txt",append=TRUE)
cat(c(loci,"\n"),file="parrots.newhybs.txt",append=TRUE)
cat("\n",file="parrots.newhybs.txt",append=TRUE)
for (i in 1:nrow(mat)){
  cat(c(i,mat[i,]),file="parrots.newhybs.txt", sep=" ",append=TRUE)
  cat("\n",file="parrots.newhybs.txt",append=TRUE)
}

file.show("parrots.newhybs.txt")









