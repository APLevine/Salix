#Pedigree recombination simulation
#APL
#2nd April 2015

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

#Number of markers
n <- chr.len/d

# 1000 meioses
# #Store output
# x <- c()
# for (i in 1:1000){
	# x <- cbind(x,rbinom(n,1,p))
# }
# #Number of recombinations per simulation
# apply(x,2,sum)

# # # library(kinship2)
# # # fam <- rbind(cbind(1,0,0,1,1),cbind(2,0,0,2,1),cbind(3,1,2,1,2))
# # # ped <- pedigree(id=fam[,1],dadid=fam[,2],momid=fam[,3],sex=fam[,4],aff=fam[,5])
# # # plot(ped)


cross <- function(haplotypes){ #perform cross over
	y <- haplotypes #input haplotypes
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
	return(z)
}

founder <- cbind(rep("A",n),rep("B",n))
cross(founder)
