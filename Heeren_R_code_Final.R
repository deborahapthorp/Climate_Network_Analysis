rm(list=ls())
getwd()
dir()


#------------PACKAGES ----------------
if (!require(huge)) install.packages("huge"); require(huge)
if (!require(networktools)) install.packages("networktools"); require(networktools)
if (!require(qgraph)) install.packages("qgraph"); require(qgraph)
if (!require(parcor)) install.packages("parcor"); require(parcor)
if (!require(Matrix)) install.packages("Matrix"); require(Matrix)
if (!require(psych)) install.packages("psych"); require(psych)
if (!require(dplyr)) install.packages("dplyr"); require(dplyr)
if (!require(ggplot2)) install.packages("ggplot2"); require(ggplot2)
if (!require(igraph)) install.packages("igraph"); require(igraph)
if (!require(car)) install.packages("car"); require(car)
if (!require(ppcor)) install.packages("ppcor"); require(ppcor)
if (!require(bootnet)) install.packages("bootnet"); require(bootnet)
if (!require(mgm)) install.packages("mgm"); require(mgm)
if (!require(corrplot)) install.packages("corrplot"); require(corrplot)
if (!require(digest)) install.packages("digest"); require(digest)
if (!require(devtools)) install.packages("devtools"); require(devtools)
if (!require(xlsx)) install.packages("xlsx"); require(xlsx)
if (!require(qgraph)) install.packages("qgraph"); require(qgraph)
if (!require(bootnet)) install.packages("bootnet"); require(bootnet)
if(!require(foreign)) {install.packages("foreign"); require(foreign)}
if(!require(qgraph)) {install.packages("qgraph"); require(qgraph)}
if(!require(grDevices)) {install.packages("grDevices"); require(grDevices)}
if(!require(ape)) {install.packages("ape"); require(ape)}
if (!require(networktools)) install.packages("networktools"); require(networktools)
if (!require(matrixcalc)) install.packages("matrixcalc"); require(matrixcalc)
if (!require(psy)) install.packages("psy"); require(psy)
if (!require(psych)) install.packages("psy"); require(psy)
if(!require(ppcor)) {install.packages("ppcor"); require(ppcor)}
if(!require(bootnet)) {install.packages("bootnet"); require(bootnet)}
if(!require(bnlearn)) {install.packages("bnlearn"); require(bnlearn)}
if(!require(corrplot)) {install.packages("corrplot"); require(corrplot)}
if(!require(smacof)) {install.packages("smacof"); require(smacof)}
if(!require(colorspace)) {install.packages("colorspace"); require(colorspace)}
if(!require(eigenmodel)) {install.packages("eigenmodel"); require(eigenmodel)}
if(!require(network)) {install.packages("network"); require(network)}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz")

if(!require(tcltk)) {install.packages("tcltk"); require(tcltk)}
if(!require(latentnet)) {install.packages("latentnet"); require(latentnet)}
if(!require(rgl)) {install.packages("rgl"); require(rgl)}

packageVersion("qgraph") #  '1.9.2'
packageVersion("networktools") # '1.5.'
packageVersion("bootnet") # '1.5.'
packageVersion("bnlearn") # '4.7.1'
packageVersion("Rgraphviz") # '2.32.0'




## -------------------- GET the data ------------------- ##
d <- read.csv("data.csv", header = TRUE, sep= ";", dec = ",", fill = TRUE)


summary(d)


data <- na.omit(d)
data<- as.data.frame(apply(data, use='pairwise.complete.obs', 2, as.numeric))
summary(data)
head(data)
tail(data)

###### description score climate anxiety ####### 
cronbach(data[,1:13]) #.89
climate_anxiety <-  rowSums(data[,1:13])
climate_anxiety
summary(climate_anxiety)
describe(climate_anxiety)

# DA: could also do this: 
# library(tidyverse)
# climate_anxiety_qs <- data %>% dplyr::select(CAS_1:CAS_13)
# cronbach(climate_anxiety_qs)
# climate_anxiety <-  rowSums(climate_anxiety_qs)

### option #2: dissociating the features ####


### cognitive_emotional
cronbach(data[,1:8]) #.89
Cognitive_Emo <-  rowSums(data[,1:8])
Cognitive_Emo
summary(Cognitive_Emo)
describe(Cognitive_Emo)

### functional
cronbach(data[,9:13]) #.89
Functional <-  rowSums(data[,9:13])
Functional
summary(Functional)
describe(Functional)

###### description score experience_climate_change####### 
cronbach(data[,14:16]) #.77
Experience <-  rowSums(data[,14:16])
Experience
summary(Experience)
describe(Experience)

###### description score behavioral_engagement####### 
cronbach(data[,17:22]) #.66
Pro_env_behav <-  rowSums(data[,17:22])
Pro_env_behav
summary(Pro_env_behav)
describe(Pro_env_behav)

###### description score Worry####### 
cronbach(data[,23:38]) #.72
Worry <-  rowSums(data[,23:38])
Worry
summary(Worry)
describe(Worry)

# creating the dataframe

DATA_combined <- as.data.frame(cbind(Cognitive_Emo,
                                     Functional,
                                     Pro_env_behav,
                                     Experience, Worry))


DATA_combined
summary(DATA_combined)


#table S1
# DA: This could almost certainly be done more nicely 
table.s1<-data.frame()
if (!require(e1071)) install.packages("e1071"); require(e1071) # for skewness
# step 3 for each variable, computing the indices of interest 
for (i in 1:ncol(DATA_combined)){
  table.s1[i,1]<-colnames(DATA_combined)[i]
  table.s1[i,2]<-round(mean(DATA_combined[,i], na.rm=TRUE),2)
  table.s1[i,3]<-round(sd(DATA_combined[,i], na.rm=TRUE),2)
  table.s1[i,4]<-min(DATA_combined[,i], na.rm=TRUE)
  table.s1[i,5]<-max(DATA_combined[,i], na.rm=TRUE)
  table.s1[i,6]<-round(skewness(DATA_combined[,i], na.rm=TRUE),2)
  table.s1[i,7]<-round(kurtosis(DATA_combined[,i], na.rm=TRUE),2)
}
colnames(table.s1)<-c("variable", "M", "SD","Min", "Max", "Skewness", "Kurtosis")
table.s1 # look at the table.s1 in the console 


#correlation (FIGURE S1) 
library(corrplot)
M <- cor(DATA_combined)
corrplot(M, method="shade") #Figure S1
#check  
is.positive.definite(cor(DATA_combined))#true 

## Use goldbricker function to search for potential "bad pairs"
## This function compares correlations in a psychometric network in order to identify
## nodes which most likely measure the same underlying construct (i.e., are colinear).
## A "bad pairs" includes two nodes that are highly correlated with each other and that 
## correlate with other nodes in the network in a highly similar pattern.
gb <- goldbricker(DATA_combined, p = 0.05, 
                  method = "hittner2003", threshold=0.25, 
                  corMin=.50, progressbar = TRUE)
## Produce list of "bad pairs"
gb

# Nonparanormal transformation
data_npn <- huge.npn(DATA_combined, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
data_npn <- as.data.frame(apply(data_npn, use='pairwise.complete.obs', 2, as.numeric))
summary(data_npn)


##############################################################################################################################
## GGM & NODE PREDICTABILITY 
##############################################################################################################################

group.ega<- c("Climate Anxiety", "Climate Anxiety", 
              "Climate Change", "Climate Change", "Worry") 
p <- ncol(data_npn)

data_npn <- as.matrix(data_npn) # DA: I had to add this line to avoid an error 

fit <- mgm(data_npn, 
           type = rep('g', p),
           lev = rep(1, p),
           rule.reg = 'OR')

pred <- predict(fit, data_npn, 
                errorCon = 'R2')
pred$error

NET <- estimateNetwork(data_npn, default="ggmModSelect", tuning = 0) 
plot(NET,node.width=1.4,layout = 'spring', 
     label.prop = 0.9, borders = TRUE, 
     pie = pred$error[,2], pieColor = rep('#377EB8',p),
     labels = colnames(data_npn),
     edge.labels=FALSE, details=FALSE, groups = group.ega,color= c("grey", "grey", 
                                                                  "yellow", "yellow", "red"), 
     vsize=12, border.width=1.9, label.cex = .9,
     label.scale.equal=TRUE, 
     color='#bbbbbb', border.width=6, legend = FALSE)

centralityPlot(NET, theme_bw=FALSE, scale = "raw0", 
               include = c("ExpectedInfluence"), decreasing= TRUE)

#####----- EBICglASSO-----#####

Network <- EBICglasso(cor(data_npn), 
                      n = nrow(data_npn), 
                      gamma = 0.5, threshold = TRUE) 

y <- qgraph(Network, layout = 'spring', 
            details = FALSE, labels=colnames(data_npn),
            vsize=16, color='#bbbbbb', border.width=2, 
            edge.labels=FALSE, label.cex = 1.2, 
            label.scale.equal=TRUE, legend = FALSE)

centralityPlot(y, theme_bw=TRUE, 
               scale = "raw0", 
               include = c("ExpectedInfluence"), 
               decreasing= TRUE)
centrality_auto(y)



##############################################################################################################################
## Bridge centrality 
##############################################################################################################################
if (!require(networktools)) install.packages("networktools"); require(networktools)
library(igraph)
if (!require(tcltk)) install.packages("tcltk"); require(tcltk)

library(rgl)
library(ape)

adjm <- getWmat(NET)
matrix<- graph_from_adjacency_matrix(adjm, mode="upper", diag=FALSE, weighted=TRUE)


## "walktrap" Community Analysis
walktrap_community <- cluster_walktrap(matrix)
head(walktrap_community)
dendPlot(walktrap_community, mode="phylo")


## bridge estimation using the communites resulting from the spinglass algorithm
library(networktools)
my_communities <- walktrap_community$membership
is.vector(my_communities)
b <-bridge(matrix,communities=my_communities)
BPLOT <- plot(b, 
              include=c("Bridge Expected Influence (1-step)"),
              theme_bw=FALSE, raw0 = TRUE, signed=TRUE, decreasing= TRUE)
BPLOT + theme(axis.text=element_text(size=5.5))




##############################################################################################################################
### STABILITY AND ACCURACCY OF THE GGM
##############################################################################################################################
#edge weights stability 
nonParametricBoot <- bootnet(data_npn, nBoots = 1000, 
                             default = "ggmModSelect", 
                             type="nonparametric", 
                             statistics="all", communities=my_communities)
plot(nonParametricBoot, labels=FALSE, order="sample")
plot(nonParametricBoot, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
summary(nonParametricBoot) %>% ungroup %>% filter(type == "edge") %>% arrange(-sample) %>% View
plot(nonParametricBoot, statistics="ExpectedInfluence", plot="difference",order = "mean")
plot(nonParametricBoot, statistics="bridgeExpectedInfluence", plot="difference",order = "mean")


#stability of the expected influence centrality 

caseDroppingBoot <- bootnet(data_npn, nBoots = 1000, 
                            default = "EBICglasso",  type="case", 
                            statistics="all", communities=my_communities)
corStability(caseDroppingBoot)
plot(caseDroppingBoot, statistics = c("ExpectedInfluence","bridgeExpectedInfluence"))
 
      

### DIFFERENTIAL VARIABILITY OF THE GGM
#  restricted range problem
c<- centrality_auto(matrix)
x <- c$node.centrality$ExpectedInfluence 
SD <- apply(data_npn,2,sd)
cor(SD,x)
cor.test(SD,x) 

######## DIRECTED ACYCLIC GRAPH (DAG) #################################################################################################

## Analyses performed in R

---- # AIM 2: Estimate Bayesian network -----------------------


## Fit a first Bayesian network, based on 50 random re-starts and 100 perturbations for each re-start 
netdata <- as.data.frame(data_npn) # DA: Need to turn it back into a data frame for this analsyis

# AIM 2: Estimate Bayesian network -----------------------
## Fit a first Bayesian network, based on 50 random re-starts and 100 perturbations for each re-start. 
set.seed(123)
fitBN1 <- hc(netdata, restart = 50, perturb = 100)  ## hc gives directed graph
fitBN1

bnlearn::score(fitBN1, data = netdata)         
astr <- arc.strength(fitBN1, netdata, "bic-g")  ## connection strength
astr[order(astr[,3]), ]  ## sorted edge strength from strongest to weakest

strength.plot(fitBN1, astr, shape = "ellipse")

## Now we stabilize the network across multiple samples through bootstrapping:
## Learn 10000 network structures (may take a few min)

# DA note: This was originally 012 as the seed - this produces a completely different result 
## to the paper, with "worry" being outside the complex with no arrows. 
## However, setting the seed as 321 produced an identical result to 123. Go figure. 
set.seed(321)
bootnet <- boot.strength(netdata, R = 10000, algorithm = "hc", 
                         algorithm.args = list(restart = 5, perturb = 10), debug = TRUE)  
head(bootnet)

# Save the bootstrapped networks
#saveRDS(bootnet, "./bootDAG.rds")

# Reload the saved bootstrapped networks
#bootnet <- readRDS("./bootDAG.rds")

## net#1 optimal cutpoint, according to Scurati & Nagarajan (2013)
avgnet <- averaged.network(bootnet)
avgnet      
bnlearn::score(avgnet, data = netdata) # -5126.457
thresh <- avgnet$learning$args[[1]]
thresh                                  ## optimal significance threshold = 0.3319
astr <- arc.strength(avgnet, netdata, "bic-g")   ## compute edge strengths
astr

suppressWarnings(strength.plot(avgnet, astr, shape = "ellipse"))
# FOR interpretation: Edge thickness signifies the probability of prediction is in the direction depicted.


## net2: use net1 threshold, edge strengths are determined by direction probability
boottab <- bootnet[bootnet$strength > thresh & bootnet$direction > 0.5, ]  ## edges in net2
boottab
astr2 <- boottab   ## table with direction probabilities
astr2$strength <- astr2$direction  ## use the direction probabilities for edge width

#strength.plot(avgnet, astr2, shape = "ellipse")
strength.plot(avgnet, astr2, shape = "ellipse", threshold = 0.50)

## thick arrows indicate high directional probabilties, thin arrows low directional probabilities
# FOR interpretation,   Edge thickness signifies the probability of prediction is in the direction depicted.

