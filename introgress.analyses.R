library(vcfR)
library(RADstackshelpR)

vcf<-read.vcfR("/Users/devder/Desktop/parrots/la.parrots/parrotsnpsunfiltered.vcf")
#read in locality info for samples
locs<-read.table("~/Desktop/parrots/la.parrots/parrotsamples.txt", header=T, sep="\t")
colnames(locs)<-c("id","location","tiss")
head(vcf@fix)
vcf@gt[1:10,1:10]
#check that sampling IDs in text file match sample IDs in the vcf
colnames(vcf@gt)[-1] == locs$id

#check out missing data
missing.by.snp(vcf)
#probably worth doing some cursory filtering for missing data

#calculate missingness per sample
mat<-extract.gt(vcf)
miss<-c()
for (i in 1:ncol(mat)){
  miss[i]<-sum(is.na(mat[,i]))
}
locs$miss<-miss
locs #elevated number of missing genotypes in toepad samples, as expected

#remove loci that aren't biallelic
mat<-mat[nchar(vcf@fix[,5]) == 1,]
dim(mat)
#identify SNPs that are fixed away from LA
conv.mat[1:5,1:5]
conv.mat<-mat
conv.mat[conv.mat == "0/0"]<-0
conv.mat[conv.mat == "0/1"]<-1
conv.mat[conv.mat == "1/1"]<-2
conv.mat<-as.data.frame(conv.mat)
#convert to numeric
for (i in 1:ncol(conv.mat)){
  conv.mat[,i]<-as.numeric(as.character(conv.mat[,i]))
}

#show colnames to verify you're subsetting correctly
#colnames(conv.mat) 
#calc AF
lilac.af<-(rowSums(conv.mat[,c(3:6)], na.rm=T)/(rowSums(is.na(conv.mat[,c(3:6)]) == FALSE)))/2
red.af<-(rowSums(conv.mat[,c(27:30)], na.rm=T)/(rowSums(is.na(conv.mat[,c(27:30)]) == FALSE)))/2

#find fixed SNPs
diff<-abs(lilac.af - red.af)
table(is.na(diff))
#how many SNPs are fixed
table(is.na(diff) == FALSE & diff == 1)

#list some fixed SNPs
head(vcf@fix[,1][is.na(diff) == FALSE & diff == 1])

#subsample original matrix to only fixed diff SNPs
gen.mat<-mat[is.na(diff) == FALSE & diff == 1,]
dim(gen.mat)
#subsample matrix converted for AF calcs to only fixed SNPS
conv.mat<-conv.mat[is.na(diff) == FALSE & diff == 1,]
dim(conv.mat)
#write a logical test to convert alleles so that a single number represents one parental ancestry
for (i in 1:nrow(gen.mat)){
  #if 1 is the red crowned allele (e.g. absent in lilacs 3-6)
  if((sum(conv.mat[i,c(3:6)], na.rm=T)/(sum(is.na(conv.mat[i,c(3:6)]) == FALSE)))/2 == 0){
    #swap all '0/0' cells with '2/2'
    gen.mat[i,][gen.mat[i,] == "0/0"]<-"2/2"
    #swap all '1/1' cells with '0/0'
    gen.mat[i,][gen.mat[i,] == "1/1"]<-"0/0"
    #finally convert all '2/2' cells (originally 0/0) into '1/1'
    gen.mat[i,][gen.mat[i,] == "2/2"]<-"1/1"
    #no need to touch hets
  }
}

#convert R class NAs to the string "NA/NA"
gen.mat[is.na(gen.mat) == TRUE]<-"NA/NA"

#if it worked correctly this should be only missing or '1/1'
table(gen.mat[,3])
table(gen.mat[,4])
table(gen.mat[,5])
table(gen.mat[,6])

#make locus info df
locus.info<-data.frame(locus=vcf@fix[vcf@fix[,3] %in% rownames(gen.mat),1],
                       type=rep("C", times=nrow(gen.mat)),
                       lg=gsub("uce-","",vcf@fix[vcf@fix[,3] %in% rownames(gen.mat),1]),
                       marker.pos=gsub("SNP_","",rownames(gen.mat)))
#make linkage group numeric
locus.info$lg<-as.numeric(as.character(locus.info$lg))
locus.info$marker.pos<-as.numeric(as.character(locus.info$marker.pos))

#we now have a gt matrix in proper format for introgress
#convert genotype data into a matrix of allele counts
count.matrix<-prepare.data(admix.gen=gen.mat, loci.data=locus.info,
                           parental1="1",parental2="0", pop.id=F,
                           ind.id=F, fixed=T)

#estimate hybrid index values
hi.index.sim<-est.h(introgress.data=count.matrix,loci.data=locus.info,
                    fixed=T, p1.allele="1", p2.allele="0")

#locus.info$locus<-rep("", times=nrow(locus.info))
#LociDataSim1$lg<-c(1:110)
mk.image(introgress.data=count.matrix, loci.data=locus.info,
         marker.order=NULL,hi.index=hi.index.sim, ylab.image="Individuals",
         xlab.h="Red-crowned ancestry", pdf=F,
         col.image=c("red","green","blue"))

#calculate mean heterozygosity
het<-calc.intersp.het(introgress.data=count.matrix)
dev.off()
#make triangle plot
plot(x=hi.index.sim$h, y=het, bg=rgb(0,0,0,alpha=0.3), pch=21, cex=2, col="black",
     xlab="Hybrid Index", ylab="Interspecific heterozygosity",
     ylim=c(0,1))
segments(x0 =0, y0 =0, x1 =.5, y1 =1)
segments(x0 =1, y0 =0, x1 =.5, y1 =1)

#calculate number of missing genotypes in our parental pops for each of the 419 fixed sites
lilac.miss<-c()
red.miss<-c()
for (i in 1:nrow(gen.mat)){
  lilac.miss[i]<-sum(gen.mat[i,3:6] == "NA/NA")
  red.miss[i]<-sum(gen.mat[i,27:30] == "NA/NA")
}

#filter gen.mat to have a minimum of 2 genotypes in both species to call a site fixed
gen.mat.filtered<-gen.mat[lilac.miss <= 2 & red.miss <= 2,]
dim(gen.mat.filtered) #214 SNPs left
locus.info.filtered<-locus.info[lilac.miss <= 2 & red.miss <= 2,]

#revisualize to see if this affects downstream inference
#convert genotype data into a matrix of allele counts
count.matrix<-prepare.data(admix.gen=gen.mat.filtered, loci.data=locus.info.filtered,
                           parental1="1",parental2="0", pop.id=F,
                           ind.id=F, fixed=T)

#estimate hybrid index values
hi.index.sim<-est.h(introgress.data=count.matrix,loci.data=locus.info.filtered,
                    fixed=T, p1.allele="1", p2.allele="0")

#plot
mk.image(introgress.data=count.matrix, loci.data=locus.info.filtered,
         marker.order=NULL,hi.index=hi.index.sim, ylab.image="Individuals",
         xlab.h="Red-crowned ancestry", pdf=F,
         col.image=c("red","green","blue"))

#this looks much better. Individuals 1-4 are Lilacs from Mexico, 5&6 are lilacs from LA,
#7-9 are apparent hybrids, 10-26 are red crowns from LA, 27:30 are red crowns from Mexico.
#white bars are UCE loci
#UCE 7322 is fascinating. Every single bird we sequenced from LA (except a single het) is homozygous lilac
#for the fixed for SNP 489. Strong evidence for adaptive introgression of this locus.
#definitely worth mining out the full sequence of this UCE and mapping it to zebra finch genome to investigate putative function

#calculate mean heterozygosity
het<-calc.intersp.het(introgress.data=count.matrix)
dev.off()
#make triangle plot
plot(x=hi.index.sim$h, y=het, bg=rgb(0,0,0,alpha=0.3), pch=21, cex=2, col="black",
     xlab="Hybrid Index", ylab="Interspecific heterozygosity",
     ylim=c(0,1))
segments(x0 =0, y0 =0, x1 =.5, y1 =1)
segments(x0 =1, y0 =0, x1 =.5, y1 =1)



#filter gen.mat to allow no missing data in either species to call a site fixed
#assess downstream effects
gen.mat.filtered<-gen.mat[lilac.miss == 0 & red.miss == 0,]
dim(gen.mat.filtered) #214 SNPs left
locus.info.filtered<-locus.info[lilac.miss == 0 & red.miss == 0,]

#revisualize to see if this affects downstream inference
#convert genotype data into a matrix of allele counts
count.matrix<-prepare.data(admix.gen=gen.mat.filtered, loci.data=locus.info.filtered,
                           parental1="1",parental2="0", pop.id=F,
                           ind.id=F, fixed=T)

#estimate hybrid index values
hi.index.sim<-est.h(introgress.data=count.matrix,loci.data=locus.info.filtered,
                    fixed=T, p1.allele="1", p2.allele="0")

#plot
mk.image(introgress.data=count.matrix, loci.data=locus.info.filtered,
         marker.order=NULL,hi.index=hi.index.sim, ylab.image="Individuals",
         xlab.h="Red-crowned ancestry", pdf=F,
         col.image=c("red","green","blue"))

#this looks cleaner, but I actually think we lose a lot of the signal showing just how messy this is.
#this reduced dataset makes it seem like almost all of the messiness is exlusively found as het sites,
#which didn't really seem to be the case when we allowed more missing data

#calculate mean heterozygosity
het<-calc.intersp.het(introgress.data=count.matrix)
dev.off()
#make triangle plot
plot(x=hi.index.sim$h, y=het, bg=rgb(0,0,0,alpha=0.3), pch=21, cex=2, col="black",
     xlab="Hybrid Index", ylab="Interspecific heterozygosity",
     ylim=c(0,1))
segments(x0 =0, y0 =0, x1 =.5, y1 =1)
segments(x0 =1, y0 =0, x1 =.5, y1 =1)

#removing the missing data certainly seems to make our triangle plot cleaner though,
#so we just have to decide which trade off we want





