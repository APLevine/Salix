#Pedigree recombination simulation
#APL
#12th April 2015

#Chromosome lengths in cM
chr <- c(280,266,227,210,212,194,192,170,164,183,157,176,133,129,138,135,140,124,115,106,70,76)

#HapMap recombination map
#http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/

#This chromosome length
chr.len <- chr[1]
#Marker cM spacing
d <- 1
#Probability of recombination between markers
p <- (1-exp(-2*d/100))/2
# since markers are spaced at 1cM intervals this is simply = .01
p <- .01

#Number of markers
n <- chr.len/d

#perform cross over
meiosis <- function(hap1,hap2){
    #input haplotypes
    y <- cbind(hap1,hap2)
    z <- y #for output
    x <- rbinom(n,1,p) #binary indicator of cross-over at each marker
    co <- which(x==1) #positions of cross-overs
    count = 0
    for (i in co){
    count = count + 1
    if (count %%2 == 1){
    z[i:n,] = y[i:n,c(2,1)]
    }
    if (count %%2 == 0){
    z[i:n,] = y[i:n,]
    }
    }
    keep <- rbinom(1,1,0.5)+1 #independent assortment
    return(z[,keep])
}

library(kinship2)

#Read pedigree
x <- read.table("pedigree_30_member.pro")
# # x=rbind(c(1,1,0,0,1,0), #input pedigree
# # c(1,2,0,0,2,0),
# # c(1,3,1,2,1,0),
# # c(1,4,1,2,2,0),
# # c(1,5,1,2,1,0),
# # c(1,6,1,2,2,0),
# # c(1,7,0,0,1,0),
# # c(1,8,7,6,1,0))
ped=pedigree(id=x[,2],dadid=x[,3],momid=x[,4],sex=x[,5],affected=x[,6],missid="0")
generation <- plot(ped)$y

#Output data frame, each individual has two columns
out <- as.data.frame(matrix(NA,nrow=n,ncol=length(x[,1])*2))
#Row names, m_i for i:n
row.names(out) <- paste("m",1:n,sep="_")
#Column names, each individual has i_1 and i_2
names(out)[1:length(x[,1])*2-1] <- paste(x[,2],1,sep="_")
names(out)[1:length(x[,1])*2] <- paste(x[,2],2,sep="_")

#Identify founders and populate haplotypes with identifier _1 and _2
founders <- which(x[,3]==0)
for (i in founders){
        pos.1 <- i*2-1
        pos.2 <- i*2
        out[,pos.1] <- paste(x[i,2],"1",sep=".")
        out[,pos.2] <- paste(x[i,2],"2",sep=".")
}

#For simulated non-founder haplotypes
sim <- out
nf <- which(x[,3]!=0) #non-founders
nf.generation <- generation[nf] #non-founders, generation
nf.o <- nf[order(nf.generation)] #non-founders ordered

for (i in nf.o){
        i.hap.pos <- c(i*2-1,i*2) #output column numbers
        for (parent in c(3,4)){ #father then mother
            count = parent-2 #1 then 2
            i.parent <- x[i,parent] #parent identifier
            i.parent.pos <- which(x[,2]==i.parent) #parent position
            i.parent.hap.1 <- sim[,i.parent.pos*2-1] #parent haplotype 1
            i.parent.hap.2 <- sim[,i.parent.pos*2] #parent haplotype 2
            parent.out <- meiosis(i.parent.hap.1,i.parent.hap.2) #meiosis result
            sim[,i.hap.pos[count]] <- parent.out #store
        }
}

#Save flow information
flow.out <- row.names(sim)
flow.out <- cbind(flow.out, sim)
write.table(flow.out,file="sim1.hap",quote=F,row.names=F,col.names=T)

#For simulated genotypes
genotypes <- sim
for (i in 1:n){
        u.haps <- unique(as.character(sim[i,])) #unique haplotypes
        for (j in u.haps){ #assign genotype
            this <- rbinom(1,1,0.5) #af=0.5
            genotypes[i,which(genotypes[i,]==j)] <- this
        }
}

#Generate tped
left <- cbind(rep(1,n),row.names(out),1:n,rep(0,n))
tped <- cbind(left,genotypes)

#Save tped
write.table(tped,file="sim1.tped",quote=F,row.names=F,col.names=F,sep="\t")

#Save tfam
write.table(x,file="sim1.tfam",quote=F,row.names=F,col.names=F,sep="\t")

