# data analysis of Phytophthora infestans population in Guatemala and Honduras

################################################

##PCA
library (polysat)

mygen <- read.GeneMapper ("PCA182_H&G.txt")


mygen

#add a description
Description(mygen) <- "Phytophthora infestans G&H"
#add population names
PopNames(mygen) <- c("H.Int","H.Lemp","H.Oco","H.FrcMor","G.SM", "G.CHT","G.SOL","G.HTG","EU13A2","US7A2")

pops <- c("H.Int","H.Lemp","H.Oco","H.FrcMor","G.SM", "G.CHT","G.SOL","G.HTG","EU13A2","US7A2")

#add population info in the same order as the samples (see file Pop.Info.txt), this can also be an additional column in the original dataset (Intibuca=1, Lempira=2, Ocotepeque=3, FranciscoMorazan=4, SanMarcos=5, Chimaltenango=6, Sololá=7, Huehuetenango=8, EU13A2=9, US7A2=10)

PopInfo(mygen) <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,3,1,1,1,1,4,1,1,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1,3,3,1,1,1,1,1,1,1,1,1,4,4,4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,5,5,6,7,8,8,8,8,8,9,9,9,10)

mygen




# indicate type of repeats (di, tri or tetra) in each locus (Pi63=3, Pi70=3, D13=2, G11=2, Pi04=2, Pi4B=2, SSR11=2, SSR2=2, SSR3=2, SSR4=2, SSR6=2, SSR8=2)

Usatnts(mygen) <- c(2,2,2,2,3,3,2,2,2,2,2,2)

# indicate that all samples are triploid
Ploidies(mygen) <- 3

# view ploidies
Ploidies(mygen)
# view the genotype array as it currently is: filled with missing values

summary(mygen)


#Genetic distances between individuals using Bruvo distance and PCA

# make Mean Pairwise Distance Matrix
Bruvomat <- meandistance.matrix(mygen, distmetric = Bruvo.distance)
Bruvomat

pca <-cmdscale(Bruvomat)
pca

#see in three dimensions
pca2 <-cmdscale(Bruvomat, k=3)


# set colors for datapoints/population
mycol <- c("#9D98F6","#0099FF","#33cc99","blue","black","firebrick4", 
           "#117733","#FFA69E","#003399","#ff66ff")




# make PCA plot
plot(pca[,1],pca2[,3], col=mycol[PopInfo(mygen)], 
     main="PCA with Bruvo distance", 
     pch=17, bg="grey96",abline(v=0.0, h=0.0))



# add reference isolate names 
text (x=-0.005948804, y=-0.019277777, "13A2_118", cex=0.8, col="#003399")
text (x=-0.03628418, y=0.00958399218, "13A2_177", cex=0.8, col="#003399")
text (x=-0.0294376121593, y=-0.00834159447, "13A2_180", cex=0.8, col="#003399")
text (x=0.278002132, y=0.13559143311, "US7A2", cex=0.8, col="#ff66ff")


#group1
#add extra isolates 
text ( x=0.2955725378,y=-0.131777780,"H126" ,cex=0.8,col="#9D98F6")
text ( x=0.2755725378,y=-0.135777780,"H155" ,cex=0.8,col="#9D98F6")
text ( x=0.2855725378,y=-0.121117777804,"G53" ,cex=0.8,col="#FFA69E")
text ( x=0.2671590938,y=-0.1394344650,"G29" ,cex=0.8,col="#117733")
text ( x=0.271240073,y=-0.1784208947,"H55" ,cex=0.8,col="#9D98F6")


#group2

text ( x=0.337353148,y=0.0528113222,"H161" ,cex=0.8,col="#9D98F6")
text ( x=0.32502472901,y=0.0877386226,"G49" ,cex=0.8,col="#FFA69E")
text ( x=0.34224254767,y=0.125793092,"G41" ,cex=0.8,col="#FFA69E")
text ( x=0.29110276232,y=0.0082549158,"H102" ,cex=0.8,col="#9D98F6")
text ( x=0.2650871335830,y=0.1207660729,"H149",cex=0.8,col="#9D98F6")
text ( x=0.2452681215,y=0.128840406,"G52" ,cex=0.8,col="#FFA69E")
text ( x=0.215984069411,y=0.06269353158,"H134" ,cex=0.8,col="#9D98F6")
text ( x=0.18515160386549,y=0.05342062586,"H120" ,cex=0.8,col="#9D98F6")


# add legend      
legend("bottomleft", c("H.Int","H.Lemp","H.Oco","H.FrcMor","G.SM", "G.CHT","G.SOL","G.HTG","EU13A2","US7A2")
       , text.col=mycol, bg="grey96", cex=0.55)

dev.off()





##############
#Tree
############################################################
## P. infestans SSR (GenAlEx) -> Bruvo distance -> BioNJ/NJ ##
############################################################

## 0) Packages
library(poppr)
library(adegenet)
library(ape)
library(ggtree)
library(ggplot2)

## 1) Read data (allows up to 3 alleles per locus)
Pinf_pc <- read.genalex("H&G182.csv", ploidy = 3, geo = FALSE, genclone = FALSE)

## 2) Recode/clean data (handles 0s and genotypes with >2 alleles when applicable)
Pinf_pc_rc <- recode_polyploids(Pinf_pc, newploidy = TRUE)

## 3) Basic QC (optional but recommended)
print(Pinf_pc_rc)
table(Pinf_pc_rc@ploidy)                # how many individuals are 2 vs 3
locNames(Pinf_pc_rc)                    # locus names
Pinf_pc_rc@loc.n.all                    # number of alleles per locus

## 4) Bruvo distance
## IMPORTANT: 'replen' must have the SAME length and order as locNames(Pinf_pc_rc)
replen_vec <- c(2,2,2,2,3,2,3,2,2,2,2,2)

stopifnot(length(replen_vec) == nLoc(Pinf_pc_rc))
names(replen_vec) <- locNames(Pinf_pc_rc)  # documents the order

bruv <- bruvo.dist(Pinf_pc_rc,
                   replen = replen_vec,
                   add = TRUE,
                   loss = TRUE,
                   by_locus = FALSE)

## 5) BioNJ tree
tree_bionj <- bionj(as.matrix(bruv))


## Quick plot (for inspection only)
plot(tree_bionj, main = "BioNJ – Bruvo distance")

## Export for editing in MEGA / FigTree / Inkscape
write.tree(tree_bionj, file = "BioNJ_Bruvo_182.nwk")




##########################################################
##Minimum Spanning Network (MSN)
#import data as genind object

library(RColorBrewer)


#add the groups into the original datafile pop slot
Pinf_pop <- read.genalex("H&G182B.csv", ploidy=3, geo = FALSE, genclone = FALSE)
Pinf_pop

#add other information, in this case the geographical info:
demographics <- read.csv("182_other.csv", header=T)
other(Pinf_pop)$xy <- demographics

#add repeats of the SSR loci
repeats <- c(2,2,2,2,3,2,3,2,2,2,2,2) #number of repeats of loci inthe same order: D13,SSR8,SSR4,Pi04,Pi70,SSR6,Pi63,G11,SSR3_Pi02,SSR11,SSR2,4B
repeats
other(Pinf_pop)$repeat_lengths <- repeats
other(Pinf_pop)

#add the dataframe from the other slot into the starta slot
strata(Pinf_pop) <- other(Pinf_pop)$xy[-1]
Pinf_pop

#remove 0s
Pinf_pop_rc <- recode_polyploids (Pinf_pop, newploidy = TRUE)
Pinf_pop_rc

#create RColorBrewer palette           
three_colors <-c("#77AADD","#ff66ff","#660066")

#msn with distance bar and other features
pcmsn <- bruvo.msn(Pinf_pop_rc, replen =repeats)

plot_poppr_msn(Pinf_pop_rc, pcmsn, size.leg = FALSE, inds=0, palette = three_colors)


dev.off()

######################################################################################
##Analysis of Molecular Variance, AMOVA 

#Hierarchical AMOVA by clonal lineage (13A2,US7A2,CA)

amova.result <-poppr.amova(Pinf_pop_rc, ~clonal.lineage/location, clonecorrect =TRUE,
                           method="ade4", within=FALSE)  

amova.result

###################################################################
##Locus stats
library("magrittr")  
locus_table(Pinf_pc_rc)
