# ===================================================
# Script Name: dynamics_of_size_and_shape_variation.R
# Author: Marko Bugarčić
# Year Created: 2025
# Description: This script performs geometric morphometric
# analysis of newt embryos' size and shape using landmark data. Results 
# are published in [too add upon acceptance].
# ===================================================


library(RRPP)
library(geomorph)
library(xlsx)
library(ggplot2)
library(abind)
library(car)

#setting work directory
setwd("C:/Podaci/Rad/mrmoljci dani/dani/r15")

#setting color scales fpr plotting
g1clr <- c('#fdae6b','#fd8d3c','#f16913','#d94801','#8c2d04')
g2clr <-c('#4292c6','#2171b5','#084594')

#data loading
# tailbud g1
g1atps <- readland.tps("g1a.tps", readcurves = F, specID = "imageID")
g1btps <- readland.tps("g1b.tps", readcurves = F, specID = "imageID")
slide1 <- read.xlsx("slajding2.xlsx",sheetIndex =  1, header = T)
slide1 <- as.matrix(slide1 [,-1])
g1otps <- readland.tps("g1.tps", readcurves = F, specID = "imageID")
# prehatching g2
g2atps <- readland.tps("g2a.tps")
g2btps <- readland.tps("g2b.tps") 
levo <- matrix(c( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 22, 23, 24))  
desno <- matrix(c( 21, 20, 19, 18, 17, 16, 15, 14, 13, 12 , 11, 25, 26, 27)) 
parovi <- cbind(levo, desno)
slide2 <- read.xlsx("slajding.xlsx",sheetIndex =  1, header = T)
slide2 <- as.matrix(slide2 [,-1])
g2otps <- readland.tps("g2.tps", readcurves = F, specID = "imageID")
#factors
faktg1  <- read.xlsx("FactorsGroup1.xlsx",sheetIndex =  1, header = T)
faktg2 <- read.xlsx("FactorsGroup2.xlsx",sheetIndex =  1, header = T)
#filtering data
s28s29 <- which(faktg2$stage %in% c(28, 29))
faktg2 <- faktg2[-s28s29, ]
g2otps<-g2otps[,,-s28s29]
g2atps<-g2atps[,,-s28s29]
g2btps<-g2btps[,,-s28s29]

#GPA
#g1
#averaging coordiates
g1tps <- (g1atps+g1btps)/2
#making data frame for tailbud
dfg1 <- geomorph.data.frame(shape= g1tps, stage = faktg1$stage, ind = faktg1$ind)
#performing gpa
gpag1 <- gpagen(dfg1$shape, curves = slide1)
#inspecting outliers whole tailbud
plotOutliers(gpag1$coords)
#inspecting outliers by stage tailbud
for (stage in unique(faktg1$stage)) {
  
  
  plotOutliers(gpag1$coords[,,faktg1$stage == stage], inspect.outliers = T)
}
#g2
#averaging coordinates
g2tps <- ((g2atps+g2btps)/2)
#making data frame for prehatch
dfg2 <- geomorph.data.frame(shape= g2tps, stage = faktg2$stage, ind = faktg2$ind)
#performing gpa with bilateral symmetry
gpag2 <- bilat.symmetry(dfg2$shape, ind = dfg2$ind, object.sym = T, land.pairs = parovi, curves = slide2)
#calculating cs
gpag2cs <- gpagen(dfg2$shape, curves = slide2)
#inspecting outliers whole prehatch
plotOutliers(gpag2$symm.shape, inspect.outliers = T)
#inspecting outliers prehatch by stage
for (stage in unique(faktg2$stage)) {
  
  
  plotOutliers(gpag2$symm.shape[,,faktg2$stage == stage], inspect.outliers = T)
}
#ME assesment by  Collyer 2024

meg2a <- g2atps
meg2b <- g2btps
#combining arrays
g1melmk <- abind(g1atps, g1btps)
g2melmk <- abind(meg2a, meg2b)
#creating factor for landmarking first and second time
tackanje1 <- c(rep(1,100),rep(2,100))
tackanje2 <- c(rep(1,60),rep(2,60))
#performing gpa tailbud
gpag1me <- gpagen(g1melmk, curves = slide1)
#performing gpa prehatch
gpag2me <- bilat.symmetry(g2melmk, ind = c(1:120), object.sym = T, land.pairs = parovi, curves = slide2)
#creating factor for individuals
ind1 <- rep(faktg1$ind,2)
ind2 <- rep(faktg2$ind,2)

#creating data frame for measurement.error
meg1df <- rrpp.data.frame( lmk = two.d.array(gpag1me$coords) , id = ind1, rep = tackanje1, stage = rep(faktg1$stage,2))
#measurement error fit
meg1new <- measurement.error(data = meg1df, Y = "lmk", subjects = "id", replicates = "rep", iter = 999, Parallel = T,
                             print.progress = F, groups = "stage")
#results for measurement error tailbud
anova(meg1new)
#visualization
plot(meg1new)
#creating data frame for measurement.error
meg2df <- rrpp.data.frame( lmk = two.d.array(gpag2me$symm.shape) , id = ind2, rep = tackanje2, stage = rep(faktg2$stage,2))
#measurement error fit
meg2new <- measurement.error(data = meg2df, Y = "lmk", subjects = "id", replicates = "rep", iter = 999, Parallel = T,
                             print.progress = F, groups = "stage")
#results of measurement error prehatch
anova(meg2new)
#visualization
plot(meg2new)

#measurement error by stages in tailbud
#creating list for results
anova_resme1 <- list()
#creating factor for stages
s1 <- rep(faktg1$stage,2)
#looping through each stage and performing measurement.error
for (i in unique(s1)){
  megidf <- rrpp.data.frame( lmk = two.d.array(gpag1me$coords[,,s1==i]) , id = ind1[s1==i],
                             rep = tackanje1[s1==i])
  meginew <- measurement.error(data = megidf, Y = "lmk", subjects = "id", replicates = "rep", iter = 999, Parallel = F,
                               print.progress = F)
  anova_resme1[[i]] <- anova(meginew)
}
anova_resme1
#output to file
writeLines(capture.output(anova_resme1), "g1mes.txt")

#measurement error by stages in prehatch
#creating list for reults
anova_resme2 <- list()
#creating factor for stages
s2 <- rep(faktg2$stage,2)
#looping thru each stage and performing measurement.error
for (i in unique(s2)){
  megidf <- rrpp.data.frame( lmk = two.d.array(gpag2me$symm.shape[,,s2==i]) , id = ind2[s2==i],
                             rep = tackanje2[s2==i])
  meginew <- measurement.error(data = megidf, Y = "lmk", subjects = "id", replicates = "rep", iter = 999, Parallel = F,
                               print.progress = F)
  anova_resme2[[i]] <- anova(meginew)
}
anova_resme2
#output to file
writeLines(capture.output(anova_resme2), "g2mes.txt")

#comparing populations in prehatch
pops <- faktg2$pop
#array to  matrix
popcoords <- two.d.array(gpag2$symm.shape)
#pairwise comparisons between levels of  poulation:stage factor
pop30 <- pairwise(lm.rrpp(popcoords[faktg2$stage==30,]~pops[faktg2$stage==30]), groups = pops[faktg2$stage==30])
pop31 <- pairwise(lm.rrpp(popcoords[faktg2$stage==31,]~pops[faktg2$stage==31]), groups = pops[faktg2$stage==31])
#results
summary_result30<-summary(pop30, confidence = 0.95, test.type = "dist", digits = 5)
summary_result30
summary_result31<-summary(pop31, confidence = 0.95, test.type = "dist", digits = 5)
summary_result31
#cluster analysis
library(HDclassif)
#g1
set.seed(1111)
g1m <- two.d.array(gpag1$coords) 
clust1 <- hddc(g1m, model ="ALL", criterion = "bic",itermax = 10000, mc.cores = 30, scaling = T)
#calculating posterior probabilities
clust1post <- predict(clust1, g1m, cls = faktg1$stage)
clust1post$posterior
probg1 <- round(clust1post$posterior*100, digits = 3)
probg1 <- data.frame(probg1[1:100,])
probg1$Stage <- faktg1$stage
probg1
#export cluster results as % of probability per cluster for each individual
write.xlsx(probg1, "g1 clus.xlsx")

#g2
g2m <- two.d.array(gpag2$symm.shape)
clust2 <- hddc(g2m, model = "ALL", criterion = "bic",itermax = 10000, mc.cores = 30, scaling = T)
#calculatin posterior probabilities
clust2post <- predict(clust2, g2m, cls =faktg2$stage)
probg2 <- round(clust2post$posterior*100, digits = 3)
probg2 <- data.frame(probg2[1:60,])
probg2$Stage <- faktg2$stage
probg2
#export cluster results as % of probability per cluster for each individual
write.xlsx(probg2, "g2 clus.xlsx")

#CS analysis
#g1
csg1 <- gpag1$Csize
#calculating means
csg1mean <- aggregate(csg1, by = list(faktg1$stage), FUN = mean)
csg1mean
#calculating variance
csg1var <- aggregate(csg1, by = list(faktg1$stage), FUN = var)
csg1var

#g2
csg2 <- gpag2cs$Csize
#calculating means
csg2mean <- aggregate(csg2, by = list(faktg2$stage), FUN = mean)
csg2mean
#calculating variance
csg2var <- aggregate(csg2, by = list(faktg2$stage), FUN = var)
csg2var


# cs mean pval
# g1
csvg1 <- aov(gpag1$Csize ~as.factor(faktg1$stage))
csvg1
summary(csvg1)
#pairwise post-hoc analysis
TukeyHSD(csvg1)
#g2
csvg2 <- aov(gpag2cs$Csize ~as.factor(faktg2$stage))
csvg2
summary(csvg2)
#pairwise post-hoc analysis
TukeyHSD(csvg2)

#g1 cs plot export
tiff('graphs/csg1.tiff', width = 2680, height = 2680, units = "px", res = 300)
mean_values <-aggregate(csg1, FUN = mean, by = list(faktg1$stage))
se_values <- aggregate(csg1, FUN = sd, by = list(faktg1$stage))
se_values[,2] <-se_values[,2] / sqrt(20)

summary_data <- data.frame(
  mean = mean_values,
  se = se_values)


p <- ggplot() +
  geom_point(data = data.frame(stage = faktg1$stage, cs = csg1), aes(x = stage, y = cs), size = 4, alpha = 0.15) +
  
  geom_errorbar(data = summary_data, aes(x = summary_data[,1], ymin = summary_data[,2] - summary_data[,4],
                                         ymax = summary_data[,2] + summary_data[,4]), width = 0.2, size = 1) + 
  geom_point(data = csg1mean, aes(x = csg1mean[,1], y = csg1mean[,2]), size = 8) +
  labs(x = "Stage", y = "Centroid size (CS)")+ theme_bw() +
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23),panel.grid = element_blank())
p + scale_x_discrete(labels = c(21,22,23,24,25))
dev.off()

#g2 cs plot export
tiff('graphs/csg2.tiff', width = 2680, height = 2680, units = "px", res = 300)
mean_values <-aggregate(csg2, FUN = mean, by = list(faktg2$stage))
se_values <- aggregate(csg2, FUN = sd, by = list(faktg2$stage))
se_values[,2] <-se_values[,2] / sqrt(20)

summary_data <- data.frame(
  mean = mean_values,
  se = se_values)


p <- ggplot() +
  geom_point(data = data.frame(stage = faktg2$stage, cs = csg2), aes(x = stage, y = cs), size = 4, alpha = 0.15) +
  
  geom_errorbar(data = summary_data, aes(x = summary_data[,1], ymin = summary_data[,2] - summary_data[,4],
                                         ymax = summary_data[,2] + summary_data[,4]), width = 0.2, size = 1) +
  geom_point(data = csg2mean, aes(x = csg2mean[,1], y = csg2mean[,2]), size = 8) +
  labs(x = " ", y = " ")+ theme_bw() +
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23),panel.grid = element_blank())
p + scale_x_discrete(labels = c(30,31,32))
dev.off()

#CS levene test
#g1
ltg1 <- leveneTest(gpag1$Csize~faktg1$stage)
ltg1

#g1 pairwise levene
levene_results <- list()


for (i in 21:24) {
  for (j in (i + 1):25) {
    subset_data <- gpag1$Csize[faktg1$stage %in% c(i, j)]
    subset_dat <- faktg1$stage[faktg1$stage %in% c(i, j)]
    ltg <- leveneTest(subset_data~ subset_dat)
    levene_results[[paste(i, j, sep = "_")]] <- ltg
  }
}
p_values <- matrix(NA, nrow = 5, ncol = 5, dimnames = list(c("21", "22", "23", "24", "25"), 
                                                           c("21", "22", "23", "24", "25")))
#formating results
for (combination in names(levene_results)) {
  stages <- strsplit(combination, "_")[[1]]
  p_value <- levene_results[[combination]]$`Pr(>F)`[1]
  p_values[stages[[1]], stages[[2]]] <- p_value
}
#levene pairwise test results
p_values


#g2
ltg2 <- leveneTest(gpag2cs$Csize~faktg2$stage)
ltg2

#g2 pairwise levene
ltp1 <- leveneTest(gpag2cs$Csize[faktg2$stage %in% c(30, 31)]~faktg2$stage[faktg2$stage %in% c(30, 31)])
ltp1

ltp2 <- leveneTest(gpag2cs$Csize[faktg2$stage %in% c(31, 32)]~faktg2$stage[faktg2$stage %in% c(31, 32)])
ltp2

ltp3 <- leveneTest(gpag2cs$Csize[faktg2$stage %in% c(30, 32)]~faktg2$stage[faktg2$stage %in% c(30, 32)])
ltp3

#PCA

#g1 pca
pcag1 <- gm.prcomp(gpag1$coords)
#export of pc scores
write.xlsx(data.frame(faktg1$ind,faktg1$stage,csg1,pcag1$x),"pcag1.xlsx")
#plot to obtain coordinates
p1 <- plot(pcag1,  pch=23, cex = 1, bg = g1clr[as.factor(faktg1$stage)],  legend = T, labels = T)


p1df <- data.frame(x=p1$plot_args$x,y=p1$plot_args$y, stage = faktg1$stage)
# p1df$x <- p1df$x*(-1)
#drawing polygons for each stage
st21 <- as.matrix(p1df[p1df$stage==21,-3])
st21pg <- chull(st21[,1], st21[,2])
st21poly <- st21[st21pg,]

st22 <- as.matrix(p1df[p1df$stage==22,-3])
st22pg <- chull(st22[,1], st22[,2])
st22poly <- st22[st22pg,]

st23 <- as.matrix(p1df[p1df$stage==23,-3])
st23pg <- chull(st23[,1], st23[,2])
st23poly <- st23[st23pg,]

st24 <- as.matrix(p1df[p1df$stage==24,-3])
st24pg <- chull(st24[,1], st24[,2])
st24poly <- st24[st24pg,]

st25 <- as.matrix(p1df[p1df$stage==25,-3])
st25pg <- chull(st25[,1], st25[,2])
st25poly <- st25[st25pg,]

pcag1plot <- pcag1

#exporting pca plot with polygons
tiff('graphs/pcag1.tiff', width = 2680, height = 2680, units = "px", res = 300)

par(mar = c(5, 5, 5, 5) + 0.1)
plot(pcag1plot,  pch=21, cex = 1.5, bg = g1clr[as.factor(faktg1$stage)],  cex.axis = 1.5, cex.lab = 2)
polygon(x = st21poly[,1], y = st21poly[,2], col = alpha(g1clr[1],  0.5),border = NA)
polygon(x = st22poly[,1], y = st22poly[,2], col = alpha(g1clr[2],  0.5),border = NA)
polygon(x = st23poly[,1], y = st23poly[,2], col = alpha(g1clr[3],  0.5),border = NA)
polygon(x = st24poly[,1], y = st24poly[,2], col = alpha(g1clr[4],  0.5),border = NA)
polygon(x = st25poly[,1], y = st25poly[,2], col = alpha(g1clr[5],  0.5),border = NA)
text(x = c(-0.15,-0.16, -0.01, 0.09,0.2),                           
     y = c(-0.14,0.09,0.121,0.09,0.01),
     labels = levels(as.factor(faktg1$stage)),
     col = "black", cex = 2)
dev.off()
par(mar = c(5, 4, 4, 2) + 0.1)

#g2 pca
pcag2 <- gm.prcomp(gpag2$symm.shape)
#export of pca scores
write.xlsx(data.frame(faktg2$ind,faktg2$stage,faktg2$pop,csg2,pcag2$x),"pcag2.xlsx") 
#pcag2$x[,2] <- pcag2$x[,2]*-1
p2 <- plot(pcag2, pch=23, cex = 1, bg = g2clr[as.factor(faktg2$stage)],  legend = T, labels = T)
p2df <- data.frame(x=p2$plot_args$x,y=p2$plot_args$y, stage = faktg2$stage)
#drawing polygons for each stage
st30 <- as.matrix(p2df[p2df$stage==30,-3])
st30pg <- chull(st30[,1], st30[,2])
st30poly <- st30[st30pg,]

st31 <- as.matrix(p2df[p2df$stage==31,-3])
st31pg <- chull(st31[,1], st31[,2])
st31poly <- st31[st31pg,]

st32 <- as.matrix(p2df[p2df$stage==32,-3])
st32pg <- chull(st32[,1], st32[,2])
st32poly <- st32[st32pg,]


par(mar = c(5, 5, 5, 5) + 0.1)
#export of prehatch pca plot with polygons
tiff('graphs/pcag2.tiff', width = 2680, height = 2680, units = "px", res = 300)
par(mar = c(5, 5, 5, 5) + 0.1)

plot(pcag2,  pch=21, cex = 1.5, bg = g2clr[as.factor(faktg2$stage)],   cex.axis = 1.5, cex.lab = 2)
polygon(x = st30poly[,1], y = st30poly[,2], col = alpha(g2clr[1],  0.5),border = NA)
polygon(x = st31poly[,1], y = st31poly[,2], col = alpha(g2clr[2],  0.5),border = NA)
polygon(x = st32poly[,1], y = st32poly[,2], col = alpha(g2clr[3],  0.5),border = NA)
text(x = c(0.046,-0.04, -0.02),                           
     y = c(0.018,0.002,-0.039)*-1,
     labels = levels(as.factor(faktg2$stage)),          
     col = "black", cex = 2)
dev.off()
par(mar = c(5, 4, 4, 2) + 0.1)


#pca plots by stage with outlier inspection
#g1
#21
coord21 <- gpag1$coords[,,faktg1$stage==21]
PCA21 <- gm.prcomp(coord21)
plot(PCA21, pch=21, cex = 3, bg = g1clr[1],  legend = T, labels = T)
text(PCA21$x, labels = faktg1$ind[faktg1$stage==21],title("s21"))
plotOutliers(coord21, inspect.outliers = T)
plotAllSpecimens(coord21, label = T)
#22
coord22 <- gpag1$coords[,,faktg1$stage==22]
PCA22 <- gm.prcomp(coord22)
plot(PCA22, pch=21, cex = 3, bg = g1clr[1],  legend = T, labels = T)
text(PCA22$x, labels = faktg1$ind[faktg1$stage==22],title("s22"))
plotOutliers(coord22, inspect.outliers = T)
plotAllSpecimens(coord22, label = T)
#23
coord23 <- gpag1$coords[,,faktg1$stage==23]
PCA23 <- gm.prcomp(coord23)
plot(PCA23, pch=21, cex = 3, bg = g1clr[1],  legend = T, labels = T)
text(PCA23$x, labels = faktg1$ind[faktg1$stage==23],title("s23"))
plotOutliers(coord23, inspect.outliers = T)
plotAllSpecimens(coord23, label = T)
#24
coord24 <- gpag1$coords[,,faktg1$stage==24]
PCA24 <- gm.prcomp(coord24)
plot(PCA24, pch=21, cex = 3, bg = g1clr[1],  legend = T, labels = T)
text(PCA24$x, labels = faktg1$ind[faktg1$stage==24],title("s24"))
plotOutliers(coord24, inspect.outliers = T)
plotAllSpecimens(coord24, label = T)
#25
coord25 <- gpag1$coords[,,faktg1$stage==25]
PCA25 <- gm.prcomp(coord25)
plot(PCA25, pch=21, cex = 3, bg = g1clr[1],  legend = T, labels = T)
text(PCA25$x, labels = faktg1$ind[faktg1$stage==25],title("s25"))
plotOutliers(coord25, inspect.outliers = T)
plotAllSpecimens(coord25, label = T)

#g2
#30
coord30 <- gpag2$symm.shape[,,faktg2$stage==30]
PCA30 <- gm.prcomp(coord30)
plot(PCA30, pch=21, cex = 3, bg = g1clr[1],  legend = T, labels = T)
text(PCA30$x, labels = faktg2$ind[faktg2$stage==30],title("s23"))
plotOutliers(coord30, inspect.outliers = T)
plotAllSpecimens(coord30, label = T)
#31
coord31 <- gpag2$symm.shape[,,faktg2$stage==31]
PCA31 <- gm.prcomp(coord31)
plot(PCA31, pch=21, cex = 3, bg = g1clr[1],  legend = T, labels = T)
text(PCA31$x, labels = faktg2$ind[faktg2$stage==31],title("s31"))
plotOutliers(coord31, inspect.outliers = T)
plotAllSpecimens(coord31, label = T)
#32
coord32 <- gpag2$symm.shape[,,faktg2$stage==32]
PCA32 <- gm.prcomp(coord32)
plot(PCA32, pch=21, cex = 3, bg = g1clr[1],  legend = T, labels = T)
text(PCA32$x, labels = faktg2$ind[faktg2$stage==32],title("s32"))
plotOutliers(coord32, inspect.outliers = T)
plotAllSpecimens(coord32, label = T)

#mean shape comparisons
#g1
pairwoblikg1 <- pairwise(lm.rrpp(two.d.array(gpag1$coords)~faktg1$stage), groups = faktg1$stage)
summary(pairwoblikg1, confidence = 0.95, test.type = "dist", digits = 5)


#g2
pairwoblikg2 <- pairwise(lm.rrpp(two.d.array(gpag2$symm.shape)~faktg2$stage), groups = faktg2$stage)
summary(pairwoblikg2, confidence = 0.95, test.type = "dist")


#function for procrustes distance of individuals
p.dist <- function (a) {
  p <- dim (a)[3]
  ms <- apply (a, c(1, 2), mean)
  dists <- vector ("numeric", p)
  for(i in 1:p) {
    dists[i] <- sqrt (sum ((a[,,i] - ms)^2))
  }
  dists
}
#g1
pd21 <- p.dist(coord21)
pd22 <- p.dist(coord22)
pd23 <- p.dist(coord23)
pd24 <- p.dist(coord24)
pd25 <- p.dist(coord25)

#g2
pd30 <- p.dist(coord30)
pd31 <- p.dist(coord31)
pd32 <- p.dist(coord32)

pd <- data.frame(pd21,pd22,pd23,pd24,pd25,
                 pd30,pd31,pd32)
pdg1 <- data.frame(pd21,pd22,pd23,pd24,pd25)
pdg2 <- data.frame(pd30,pd31,pd32)

library("ggplot2")
# export plots of procrustes distance
tiff('graphs/pdg1.tiff', width = 2680, height = 2680, units = "px", res = 300)
mean_values <- colMeans(pdg1)
se_values <- apply(pdg1, 2, sd) / sqrt(nrow(pdg1))


summary_data <- data.frame(
  variable = names(mean_values),
  mean = mean_values,
  se = se_values)


library(reshape2)
pdg1_melted <- melt(pdg1)


p <- ggplot() +
  geom_point(data = pdg1_melted, aes(x = variable, y = value), size = 4,alpha = 0.15) +
  
  geom_errorbar(data = summary_data, aes(x = variable, ymin = mean - se, ymax = mean + se), width = 0.2, size = 1) +
  geom_point(data = summary_data, aes(x = variable, y = mean), size = 8) +
  labs(x = "Stage", y = "Procrustes distance")+ theme_bw() +
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23),panel.grid = element_blank())
p + scale_x_discrete(labels = c(21,22,23,24,25))
dev.off()


tiff('graphs/pdg2.tiff', width = 2680, height = 2680, units = "px", res = 300)
mean_values <- apply(pdg2, 2, mean)
median_values <- apply(pdg2, 2, median)
se_values <- apply(pdg2, 2, sd) / sqrt(nrow(pdg2))


summary_data <- data.frame(
  variable = names(mean_values),
  mean = mean_values,
  median = median_values,
  se = se_values)

pdg2_melted <- melt(pdg2)


p <- ggplot() +
  geom_point(data = pdg2_melted, aes(x = variable, y = value), size = 4, alpha = 0.15) +
  
  geom_errorbar(data = summary_data, aes(x = variable, ymin = mean - se, ymax = mean + se), width = 0.2, size = 1) +
  geom_point(data = summary_data, aes(x = variable, y = mean), size = 8) +
  labs(x = " ", y = " ")+ theme_bw() +
  theme(axis.text.x = element_text(size = 23),
        axis.text.y = element_text(size = 23),
        axis.title.x = element_text(size = 23),
        axis.title.y = element_text(size = 23),panel.grid = element_blank())
p + scale_x_discrete(labels = c(30,31,32))
dev.off()

#regressions of shape on developmental stage
#stages~pc1
#g1 model fit
we1 <- procD.lm(gpag1$coords~faktg1$stage)
#export graph
tiff('graphs/rsg1.tiff', width = 2980, height = 1980, units = "px", res = 300)
par(mar = c(5,5, 5,5) + 0.1)
#obtain scores
sr1 <- plotAllometry(we1, size = as.numeric(as.factor(faktg1$stage)), logsz = F,method = "PredLine", col = g1clr[as.factor(faktg1$stage)],
                     bg =  g1clr[as.factor(faktg1$stage)],pch = 21,xlab="Stage", 
                     ylab="Centroid size", cex = 4, xaxt="n", ylim = c(-0.4, 0.21),cex.lab= 1.5, cex.axis = 1)
#obtain mean values of scores
wey1 <- aggregate(sr1$plot_args$y, FUN = mean, by = list(sr1$plot_args$x))
wey1
#draw individuals
points (x = sr1$plot_args$x, y = sr1[["RegScore"]], pch = 21,
        bg =   g1clr[as.factor(faktg1$stage)], cex = 2)

axis(side = 1, at=c(1,	2,	3,	4, 5), labels = levels(as.factor(faktg1$stage)),cex.lab= 1.5, cex.axis = 1.5)
#draw means for stages
lines(x=c(1:5) , y=wey1$x, lwd = 3)
dev.off()

#g2 model fit
we2 <- procD.lm(gpag2$symm.shape~faktg2$stage)
par(mar = c(5,5, 5,5) + 0.1)
#export graph
tiff('graphs/rsg2.tiff', width = 2980, height = 1980, units = "px", res = 300)
par(mar = c(5,5, 5,5) + 0.1)
#get data for plotting
sr2<- plotAllometry(we2, size = as.numeric(as.factor(faktg2$stage)), logsz = F,method = "PredLine", col =  g2clr[as.factor(faktg2$stage)],
                    bg =  g2clr[as.factor(faktg2$stage)],pch = 21,xlab="Stage",
                    ylab="Centroid size", cex = 4, xaxt="n", ylim = c(-0.09, 0.05),cex.lab= 1.5, cex.axis = 1)
#calculate means
wey2 <- aggregate(sr2$plot_args$y, FUN = mean, by = list(sr2$plot_args$x))
wey2
#draw individuals
points (x = sr2$plot_args$x, y = sr2[["RegScore"]], pch = 21,
        bg =   g2clr[as.factor(faktg2$stage)], cex = 2)

axis(side = 1, at=c(1,	2,	3), labels = levels(as.factor(faktg2$stage)),cex.lab= 1.5, cex.axis = 1.5)
#draw means
lines(x=c(1:3) , y=wey2$x, lwd = 3)
dev.off()
par(mar = c(5, 4, 4, 2) + 0.1)

#shape changes stage by stage
#g1
#s21

#s22
tiff('graphs/Stage 21 22.tiff', width = 2680, height = 2680, units = "px", res = 300)
plotRefToTarget(mshape(coord21),mshape(coord22), method = "vector", mag = 1,
                gridPars=gridPar( pt.size = 2, tar.link.lwd = 3,link.lwd = 3, grid.lwd = 3))
#title("Stage 21 to 22",  cex.main  = 2)
dev.off()
#s23
tiff('graphs/Stage 22 23.tiff', width = 2680, height = 2680, units = "px", res = 300)
plotRefToTarget(mshape(coord22),mshape(coord23), method = "vector", mag = 1,
                gridPars=gridPar( pt.size = 2, tar.link.lwd = 3,link.lwd = 3, grid.lwd = 3))
#title("Stage 23 to 23",  cex.main  = 2)
dev.off()
#s24
tiff('graphs/Stage 23 24.tiff', width = 2680, height = 2680, units = "px", res = 300)
plotRefToTarget(mshape(coord23),mshape(coord24), method = "vector", mag = 1,
                gridPars=gridPar( pt.size = 2, tar.link.lwd = 3,link.lwd = 3, grid.lwd = 3))
#title("Stage 23 to 24",  cex.main  = 2)
dev.off()
#s25
tiff('graphs/Stage 24 25.tiff', width = 2680, height = 2680, units = "px", res = 300)
plotRefToTarget(mshape(coord24),mshape(coord25), method = "vector", mag = 1,
                gridPars=gridPar( pt.size = 2, tar.link.lwd = 3,link.lwd = 3, grid.lwd = 3))
#title("Stage 24 to 25",  cex.main  = 2)
dev.off()

tiff('graphs/Stage 23 25.tiff', width = 2680, height = 2680, units = "px", res = 300)
plotRefToTarget(mshape(coord23),mshape(coord25), method = "vector", mag = 1,
                gridPars=gridPar( pt.size = 2, tar.link.lwd = 3,link.lwd = 3, grid.lwd = 3))
dev.off()

#g2
#s31
tiff('graphs/Stage 30 31.tiff', width = 2680, height = 2680, units = "px", res = 300)
plotRefToTarget(mshape(coord30),mshape(coord31), method = "vector", mag = 5,
                gridPars=gridPar( pt.size = 2, tar.link.lwd = 3,link.lwd = 3, grid.lwd = 3))
dev.off()
#s32
tiff('graphs/Stage 31 32.tiff', width = 2680, height = 2680, units = "px", res = 300)
plotRefToTarget(mshape(coord31),mshape(coord32), method = "vector", mag = 5,
                gridPars=gridPar( pt.size = 2, tar.link.lwd = 3,link.lwd = 3, grid.lwd = 3))
dev.off()
#30 to 32
tiff('graphs/Stage 30 32.tiff', width = 2680, height = 2680, units = "px", res = 300)
plotRefToTarget(mshape(coord30),mshape(coord32), method = "vector", mag = 3,
                gridPars=gridPar( pt.size = 2, tar.link.lwd = 3,link.lwd = 3, grid.lwd = 3))
dev.off()

#shape variation
#g1
vg1 <- morphol.disparity(gpag1$coords~faktg1$stage, groups = faktg1$stage , iter = 1000)
vg1 #
#alpha 0.005



#g2
vg2 <- morphol.disparity(gpag2$symm.shape ~faktg2$stage, groups = faktg2$stage , iter = 1000)
vg2
#alpha 0.0167






#regression shape~cs*stage
#g1
rg1 <- procD.lm(gpag1$coords~log(gpag1$Csize)*as.factor(dfg1$stage))
summary(rg1)

#g2
rg2 <- procD.lm(gpag2$symm.shape ~log(gpag2cs$Csize)*as.factor(dfg2$stage))
summary(rg2)


#regressions by stage
#g1
par(mfrow = c(5, 1))
#s21
rg21 <- procD.lm(gpag1$coords[,,dfg1$stage==21]~log(gpag1$Csize[dfg1$stage==21]))
summary(rg21)
as21<-plotAllometry(rg21, size = gpag1$Csize[dfg1$stage==21], logsz = T, main="Stage 21")

plot(x=as21$plot_args$x, y = as21[["RegScore"]], pch = 21, xlab = "log(CS)", ylab = "Regression scores", cex=2)


#s22
rg22 <- procD.lm(gpag1$coords[,,dfg1$stage==22]~log(gpag1$Csize[dfg1$stage==22]))
summary(rg22)
as22<-plotAllometry(rg22, size = gpag1$Csize[dfg1$stage==22], logsz = T, main="Stage 22")

plot(x=as22$plot_args$x, y = as22[["RegScore"]], pch = 21, xlab = "log(CS)", ylab = "Regression scores", cex=2)

#s23
rg23 <- procD.lm(gpag1$coords[,,dfg1$stage==23]~log(gpag1$Csize[dfg1$stage==23]))
summary(rg23)
as23<-plotAllometry(rg23, size = gpag1$Csize[dfg1$stage==23], logsz = T, main="Stage 23")

plot(x=as23$plot_args$x, y = as23[["RegScore"]], pch = 21, xlab = "log(CS)", ylab = "Regression scores", cex=2)

#s24
rg24 <- procD.lm(gpag1$coords[,,dfg1$stage==24]~log(gpag1$Csize[dfg1$stage==24]))
summary(rg24)
as24<-plotAllometry(rg24, size = gpag1$Csize[dfg1$stage==24], logsz = T, main="Stage 24")

plot(x=as24$plot_args$x, y = as24[["RegScore"]], pch = 21, xlab = "log(CS)", ylab = "Regression scores", cex=2)

#s25
rg25 <- procD.lm(gpag1$coords[,,dfg1$stage==25]~log(gpag1$Csize[dfg1$stage==25]))
summary(rg25)
as25<-plotAllometry(rg25, size = gpag1$Csize[dfg1$stage==25], logsz = T, main="Stage 25")

#summary of results for tailbud
regadjpvalg1 <- data.frame(Stage = c(21:25),
                           Rsq = c(rg21$aov.table$Rsq[1],rg22$aov.table$Rsq[1],
                                   rg23$aov.table$Rsq[1],rg24$aov.table$Rsq[1],rg25$aov.table$Rsq[1]), 
                           Pval = c(rg21$aov.table$`Pr(>F)`[1],rg22$aov.table$`Pr(>F)`[1],
                                    rg23$aov.table$`Pr(>F)`[1],rg24$aov.table$`Pr(>F)`[1],rg25$aov.table$`Pr(>F)`[1]))
regadjpvalg1$Rsq <- round(regadjpvalg1$Rsq, digits = 4)
regadjpvalg1
#alpha 0.01
#logCS reg score g1
par(mfrow = c(5, 1))
plot(x=as21$plot_args$x, y = as21[["RegScore"]], pch = 21,
     xlab = "log(CS)", ylab = "Regression scores", cex=2,
     main=paste("Stage 21 Rsq = ", regadjpvalg1[1,2],"P-value = ", regadjpvalg1[1,3], sep = " "), bg = g1clr[1])
plot(x=as22$plot_args$x, y = as22[["RegScore"]], pch = 21,
     xlab = "log(CS)", ylab = "Regression scores", cex=2,
     main=paste("Stage 22 Rsq = ", regadjpvalg1[2,2],"P-value = ", round(regadjpvalg1[2,3], digits = 3), sep = " "), bg = g1clr[2])
plot(x=as23$plot_args$x, y = as23[["RegScore"]], pch = 21,
     xlab = "log(CS)", ylab = "Regression scores", cex=2,
     main=paste("Stage 23 Rsq = ", regadjpvalg1[3,2],"P-value = ", regadjpvalg1[3,3], sep = " "), bg = g1clr[3])
plot(x=as24$plot_args$x, y = as24[["RegScore"]], pch = 21,
     xlab = "log(CS)", ylab = "Regression scores", cex=2,
     main=paste("Stage 24 Rsq = ", regadjpvalg1[4,2],"P-value = ", regadjpvalg1[4,3], sep = " "), bg = g1clr[4])
plot(x=as25$plot_args$x, y = as25[["RegScore"]], pch = 21,
     xlab = "log(CS)", ylab = "Regression scores", cex=2,
     main=paste("Stage 25 Rsq = ", regadjpvalg1[5,2],"P-value = ", regadjpvalg1[5,3], sep = " "), bg = g1clr[5])


par(mfrow = c(1, 1))
#alpha 0.01

#g2
#s30
rg30 <- procD.lm(gpag2$symm.shape [,,dfg2$stage==30]~log(gpag2cs$Csize[dfg2$stage==30]))
summary(rg30)
as30<-plotAllometry(rg30, size = gpag2cs$Csize[dfg2$stage==30], logsz = T, main="Stage 30")

plot(x=as30$plot_args$x, y = as30[["RegScore"]], pch = 21, xlab = "log(CS)", ylab = "Regression scores", cex=2)


#s31
rg31 <- procD.lm(gpag2$symm.shape[,,dfg2$stage==31]~log(gpag2cs$Csize[dfg2$stage==31]))
summary(rg31)
as31<-plotAllometry(rg31, size = gpag2cs$Csize[dfg2$stage==31], logsz = T, main="Stage 31")

plot(x=as31$plot_args$x, y = as31[["RegScore"]], pch = 21, xlab = "log(CS)", ylab = "Regression scores", cex=2)

#s32
rg32 <- procD.lm(gpag2$symm.shape[,,dfg2$stage==32]~log(gpag2cs$Csize[dfg2$stage==32]))
summary(rg32)
as32<-plotAllometry(rg32, size = gpag2cs$Csize[dfg2$stage==32], logsz = T, main="Stage 32")

plot(x=as32$plot_args$x, y = as32[["RegScore"]], pch = 21, xlab = "log(CS)", ylab = "Regression scores", cex=2)
#summary of results for prehatch
regadjpvalg2 <- data.frame(Stage = c(30:32),
                           Rsq = c(rg30$aov.table$Rsq[1],rg31$aov.table$Rsq[1],
                                   rg32$aov.table$Rsq[1]), 
                           Pval = c(rg30$aov.table$`Pr(>F)`[1],rg31$aov.table$`Pr(>F)`[1],
                                    rg32$aov.table$`Pr(>F)`[1]))
regadjpvalg2$Rsq <- round(regadjpvalg2$Rsq, digits = 4)
regadjpvalg2
#alpha
0.05/3
par(mfrow = c(5, 1))

plot(x=as30$plot_args$x, y = as30[["RegScore"]], pch = 21,
     xlab = "log(CS)", ylab = "Regression scores", cex=2,
     main=paste("Stage 30 Rsq = ", regadjpvalg2[1,2],"P-value = ", round(regadjpvalg2[1,3], digits = 3), sep = " "), bg = g2clr[1])
plot(x=as31$plot_args$x, y = as31[["RegScore"]], pch = 21,
     xlab = "log(CS)", ylab = "Regression scores", cex=2,
     main=paste("Stage 31 Rsq = ", regadjpvalg2[2,2],"P-value = ", round(regadjpvalg2[2,3], digits = 3), sep = " "), bg = g2clr[2])
plot(x=as32$plot_args$x, y = as32[["RegScore"]], pch = 21,
     xlab = "log(CS)", ylab = "Regression scores", cex=2,
     main=paste("Stage 32 Rsq = ", regadjpvalg2[3,2],"P-value = ", round(regadjpvalg2[3,3], digits = 3), sep = " "), bg = g2clr[3])
par(mfrow = c(1, 1))

#non allometric shape
#s23
s23sr <- arrayspecs(rg23$residuals,p=dim(gpag1$coords[,,dfg1$stage==23])[1],k=dim(gpag1$coords[,,dfg1$stage==23])[2])
s23adj <- s23sr + array(mshape(gpag1$coords[,,dfg1$stage==23]), dim(s23sr)) 

morphol.disparity(s23adj~1)
morphol.disparity(gpag1$coords[,,dfg1$stage==23]~1)
v23adj <- abind(s23adj,gpag1$coords[,,dfg1$stage==23])
dimnames(v23adj)[[3]]<- c(1:40)
#comparing variance of shape before and after removal of allometry
morphol.disparity(v23adj~factor(c(rep("A", 20),rep("B",20))), groups = factor(c(rep("A", 20),rep("B",20))) , iter = 1000)

