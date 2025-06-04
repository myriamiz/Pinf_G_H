# data analysis of Phytophthora infestans population in Guatemala and Honduras

################################################

##PCA
library (polysat)

mygen <- read.GeneMapper ("PCA243_G&H.txt")
mygen

#add a description
Description(mygen) <- "Phytophthora infestans G&H"
#add population names
PopNames(mygen) <- c("H.Int","H.FrcMor","H.Oco","H.Lemp","G.SOL","G.QTZ","G.HTG","G.SM", "G.CHT","13_A2","US8A2","US7A2")

pops <- c("H.Int","H.FrcMor","H.Oco","H.Lemp","G.SOL","G.QTZ","G.HTG","G.SM", "G.CHT","13_A2","US8A2","US7A2")

#add population info in the same order as the samples (see file Pop.Info.txt), this can also be an additional column in the original dataset (H.Int=1, H.FrcMor=2, H.Oco=3, H.Lemp=4, G.SOL=5, G.QTZ=6, G.HTG=7, G.SM=8, G.CHT=9, 13_A2=10, US8A2=11, US7A2=12)

PopInfo(mygen) <- c(1,1,1,2,3,1,1,1,1,1,1,1,3,3,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,4,1,1,1,1,1,1,1,1,1,5,6,6,1,1,1,1,1,1,1,1,1,1,1,7,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,6,1,1,1,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,10,1,1,10,1,1,2,1,1,1,1,1,1,1,6,1,1,1,1,8,6,1,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1,9,1,1,10,1,1,1,1,1,1,1,1,6,1,1,1,7,1,1,1,1,1,1,1,1,1,1,1,1,1,5,1,1,1,1,1,1,1,1,1,10,1,8,1,1,1,1,1,1,1,1,1,9,11,7,8,12,1,8,1,6,1,1,5,7,8,1,6,7,1,9,8,8,7,8,9,6,7,7,9,8,6,6,1,6,7,1,1,6,8,6,5,6,6,7,1,1,1)

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
           "#93E1D8","#117733","#FFA69E","#003399","#ffcc66","#ff66ff")#"#660066","purple4","royalblue1","hotpink1","lightskyblue3" )



# make PCA plot
plot(pca[,1],pca2[,3], col=mycol[PopInfo(mygen)], 
     main="PCA with Bruvo distance", 
     pch=17, bg="grey96",abline(v=0.0, h=0.0))



# add reference isolate names 
text (x=0.048172839, y=0.0122837225, "13A2_1", cex=0.8, col="#003399")
text (x=0.05323101864, y=0.00130165859, "13A2_118", cex=0.8, col="#003399")
text (x=0.059224174799, y=-0.0097606472, "13A2_177", cex=0.8, col="#003399")
text (x=0.085781333, y=-0.0001527002, "13A2_180", cex=0.8, col="#003399")
text (x=-0.12033415597, y=0.03175248, "US8A2", cex=0.8, col="#ffcc66")
text (x=-0.26289181, y=-0.110283815, "US7A2", cex=0.8, col="#ff66ff")



# add legend      
legend("topright", c("H.Int","H.FrcMor","H.Oco","H.Lemp","G.SOL","G.QTZ","G.HTG","G.SM", "G.CHT","13_A2","US8A2","US7A2")
       , text.col=mycol, bg="grey96", cex=0.55)

dev.off()





##########################################################
##Tree


library(poppr)

Pinf_pc <- read.genalex("G&H243.csv", ploidy=3, geo = FALSE, genclone = FALSE)
Pinf_pc
#remove 0s
Pinf_pc_rc <- recode_polyploids (Pinf_pc, newploidy = TRUE)
Pinf_pc_rc


##estimating distances##
##BRUVO DISTANCE##
bruv <- bruvo.dist(Pinf_pc_rc,replen=c(2,2,2,2,3,2,3,2,2,2,2,2),add = TRUE, loss = TRUE, by_locus = FALSE)
bruv

##estimate UPGMA trees with each distance##

##bruvo's UPGMA tree
UPGMA_bruv <- upgma(bruv)

write.tree(UPGMA_bruv, file="UPGMA_bruvoTree.nwk")


##########################################################
##Minimum Spanning Network (MSN)
#import data as genind object

library(RColorBrewer)


#add the groups into the original datafile pop slot
Pinf_pop <- read.genalex("G&H243B.csv", ploidy=3, geo = FALSE, genclone = FALSE)
Pinf_pop

#add other information, in this case the geographical info:
demographics <- read.csv("243_other.csv", header=T)
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
four_colors <-c("#77AADD","#ffcc66","#ff66ff","#660066")

#msn with distance bar and other features
pcmsn <- bruvo.msn(Pinf_pop_rc, replen =repeats)

plot_poppr_msn(Pinf_pop_rc, pcmsn, size.leg = FALSE, inds=0, palette = four_colors)


dev.off()

######################################################################################
#Analysis of Molecular Variance, AMOVA 

#Hierarchical AMOVA by clonal lineage (13A2,US8A2,US7A2,CA)

amova.result <-poppr.amova(Pinf_pop_rc, ~clonal_lineage/location, clonecorrect =TRUE,
                           method="ade4", within=FALSE)  

amova.result

####################################################################################
##Locus stats
library("magrittr")  
locus_table(Pinf_pc_rc)


