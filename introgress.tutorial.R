#introgress tutorial
install.packages("introgress")
library(introgress)

#setwd
setwd("~/Desktop/parrots")

#read in data for individuals from admixed pop
AdmixDataSim1<-read.csv(file="AdmixDataSim1.txt", header = F)

#read in marker info
LociDataSim1<-read.csv(file="LociDataSim1.txt", header = T)

#look at help
help("AdmixDataSim1")
dim(AdmixDataSim1) #rows are individuals, columns are SNPs
dim(LociDataSim1) #tidy format info about the SNPs

#convert genotype data into a matrix of allele counts
count.matrix<-prepare.data(admix.gen=AdmixDataSim1, loci.data=LociDataSim1,
                           parental1="P1",parental2="P2", pop.id=F,
                           ind.id=F, fixed=T)

#estimate hybrid index values
hi.index.sim<-est.h(introgress.data=count.matrix,loci.data=LociDataSim1,
                    fixed=T, p1.allele="P1", p2.allele="P2")

LociDataSim1$locus<-rep("", times=nrow(LociDataSim1))
LociDataSim1$lg<-c(1:110)
mk.image(introgress.data=count.matrix, loci.data=LociDataSim1,
         marker.order=NULL,hi.index=hi.index.sim, ylab.image="Individuals",
         xlab.h="population 2 ancestry", pdf=F,
         col.image=c("lavender","grey","red"))

#calculate mean heterozygosity across these 110 fixed markers for each sample
#using their function
het<-calc.intersp.het(introgress.data=count.matrix)

#or do it yourself
mat<-count.matrix$Count.matrix
het<-c()
for (i in 1:ncol(mat)){
  het[i]<-length(mat[,i][mat[,i] ==1])/sum(is.na(mat[,i]) == FALSE)
}

#make triangle plot
triangle.plot(hi.index=hi.index.sim, int.het=het, pdf=T, out.file="tri.pdf")
triangle.plot(hi.index=hi.index.sim, int.het=het, pdf = F)

plot(x=hi.index.sim$h, y=het, bg=rgb(0,0,0,alpha=0.3), pch=21, cex=2, col="black",
     xlab="Hybrid Index", ylab="Interspecific heterozygosity",)
segments(x0 =0, y0 =0, x1 =.5, y1 =1)
segments(x0 =1, y0 =0, x1 =.5, y1 =1)


#see if this works for California and Woodhouse's Scrub-Jays
library(vcfR)
vcf<-read.vcfR("/Users/devder/Desktop/aph.data/unzipped.filtered.vcf")
#read in locality info for samples
locs<-read.csv("~/Desktop/aph.data/rad.sampling.csv")
#subset locs to include only samples that passed filtering, and have been retained in the vcf
locs<-locs[locs$id %in% colnames(vcf@gt),]
rownames(locs)<-1:95

#remove sumi, island, florida
locs<-locs[c(1:32,57:95),]
vcf@gt<-vcf@gt[,c(1:33,58:96)]
colnames(vcf@gt)[-1] == locs$id
rownames(locs)<-1:71

#drop 20, 26:32, 43, 50:71
locs<-locs[c(1:19,21:25,33:42,44:49),]
rownames(locs)<-1:40
vcf@gt<-vcf@gt[,c(1:20,22:26,34:43,45:50)]
colnames(vcf@gt)[-1] == locs$id

#now have only california from CA & oregon, and wood from Nevada, CO, utah, NM, and AZ
#identify SNPs that are fixed away from the hybrid zone (Oregon vs. NM/CO)
mat<-extract.gt(vcf)
mat[1:5,1:5]
conv.mat<-mat
conv.mat[conv.mat == "0/0"]<-0
conv.mat[conv.mat == "0/1"]<-1
conv.mat[conv.mat == "1/1"]<-2
conv.mat<-as.data.frame(conv.mat)
#convert to numeric
for (i in 1:ncol(conv.mat)){
  conv.mat[,i]<-as.numeric(as.character(conv.mat[,i]))
}

#calc AF
oregon.af<-(rowSums(conv.mat[,c(1:5)], na.rm=T)/(rowSums(is.na(conv.mat[,c(1:5)]) == FALSE)))/2
nm.co.af<-(rowSums(conv.mat[,c(35:40)], na.rm=T)/(rowSums(is.na(conv.mat[,c(35:40)]) == FALSE)))/2

#find fixed SNPs
diff<-abs(oregon.af - nm.co.af)
table(diff == 1)
vcf@fix[,1][diff == 1]

#subsample original matrix to only fixed diff SNPs
gen.mat<-mat[diff==1,]

#flip rows that aren't allele 1=CA, to be allele 1=CA (2,14,15,23,24,28,29,30)
gen.mat[c(2,14,15,23,24,28,29,30),][gen.mat[c(2,14,15,23,24,28,29,30),] == "0/0"]<-"2/2"
gen.mat[c(2,14,15,23,24,28,29,30),][gen.mat[c(2,14,15,23,24,28,29,30),] == "1/1"]<-"0/0"
gen.mat[c(2,14,15,23,24,28,29,30),][gen.mat[c(2,14,15,23,24,28,29,30),] == "2/2"]<-"1/1"
gen.mat[is.na(gen.mat) == TRUE]<-"NA/NA"
gen.mat<-as.data.frame(gen.mat)

#make locus info df
locus.info<-data.frame(locus=rownames(gen.mat),
                       type=rep("C", times=nrow(gen.mat)),
                       lg=vcf@fix[,1][diff == 1],
                       marker.pos=1:nrow(gen.mat))
#make linkage group numeric
locus.info$lg<-gsub("Pseudochr", "", locus.info$lg)
locus.info$lg[locus.info$lg == "M"]<-10
locus.info$lg[locus.info$lg == "Z"]<-11
locus.info$lg[locus.info$lg == "1A"]<-1.5

#we now have a gt matrix in proper format for introgress
#convert genotype data into a matrix of allele counts
count.matrix<-prepare.data(admix.gen=gen.mat, loci.data=locus.info,
                           parental1="1",parental2="0", pop.id=F,
                           ind.id=F, fixed=T)

#estimate hybrid index values
hi.index.sim<-est.h(introgress.data=count.matrix,loci.data=locus.info,
                    fixed=T, p1.allele="1", p2.allele="0")

locus.info$locus<-rep("", times=nrow(locus.info))
#LociDataSim1$lg<-c(1:110)
mk.image(introgress.data=count.matrix, loci.data=locus.info,
         marker.order=NULL,hi.index=hi.index.sim, ylab.image="Individuals",
         xlab.h="population 2 ancestry", pdf=F,
         col.image=c("green","grey","blue"))

#calculate mean heterozygosity across these 110 fixed markers for each sample
#using their function
het<-calc.intersp.het(introgress.data=count.matrix)
dev.off()
#make triangle plot
triangle.plot(hi.index=hi.index.sim, int.het=het, pdf = F)

plot(x=hi.index.sim$h, y=het, bg=rgb(0,0,0,alpha=0.3), pch=21, cex=2, col="black",
     xlab="Hybrid Index", ylab="Interspecific heterozygosity",
     ylim=c(0,1))
segments(x0 =0, y0 =0, x1 =.5, y1 =1)
segments(x0 =1, y0 =0, x1 =.5, y1 =1)

#make color vector
col.vec<-(c(rep("black", times=5),rep("red", times=4),rep("green", times=5),
            rep("gray", times=5),rep("blue", times=5),rep("yellow", times=2),
            rep("pink", times=3),rep("brown", times=3), rep("orange", times=2),
            rep("black", times=6)))
#plot colored by sampling locality
plot(x=hi.index.sim$h, y=het, col=col.vec, pch=21, cex=2,
     xlab="Hybrid Index", ylab="Interspecific heterozygosity",
     ylim=c(0,1))
segments(x0 =0, y0 =0, x1 =.5, y1 =1)
segments(x0 =1, y0 =0, x1 =.5, y1 =1)
legend("topright", inset=.02,
       c("norcal","bay","centralcal","bigbear","AZ","nevada","eastcal","utah"),
       fill=c("red","green","gray","blue","yellow","pink","brown","orange"), cex=0.6,)

colnames(gen.mat) == locs$id






