#------------PACKAGES ----------------
require(huge)
require(networktools)
require(qgraph)
#if (!require(parcor)) install.packages("parcor"); require(parcor)
require(Matrix)
require(psych)
require(dplyr)
require(ggplot2)
require(igraph)
require(car)
require(ppcor)
require(bootnet)
require(mgm)
require(corrplot)
require(digest)
require(devtools)
require(xlsx)
require(qgraph)
require(bootnet)
require(foreign)
require(qgraph)
require(grDevices)
require(ape)
require(networktools)
require(matrixcalc)
require(psy)
require(psych)
require(ppcor)
require(bootnet)
require(bnlearn)
require(corrplot)
require(smacof)
require(colorspace)
require(eigenmodel)
require(network)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rgraphviz")

require(tcltk)
require(latentnet)
require(rgl)

library(haven)

youth_climate_data <- read_sav("Youth Mental Health and Climate Change_2023_IMPUTED_FOR ANALYSIS_wID_877_01_02_2026.sav")

youth_climate_data<-na.omit(youth_climate_data) # remove NAs

youth_climate_data$DASS_total <- (youth_climate_data$DASS_anxiety+youth_climate_data$DASS_depression+youth_climate_data$DASS_stress)

Cognitive_Emo <-  youth_climate_data$CCAS_CE_impairment
Functional <-  youth_climate_data$CCAS_functional_impair
Experience <-  youth_climate_data$CCAS_CC_experience
Pro_env_behav <-  youth_climate_data$CCAS_behav_engage
Worry <- (youth_climate_data$DASS_total) # Was using DASS anxiety - switched to total on Amy's advice

DATA_combined <- as.data.frame(cbind(Cognitive_Emo,
                                     Functional,
                                     Pro_env_behav,
                                     Experience, Worry))

summary(DATA_combined)

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


M <- cor(DATA_combined)
corrplot(M, method="shade") #Figure S1

check <- is.positive.definite(cor(DATA_combined))#true 

check

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

#transform data with nonparanormal transform

data_npn <- huge.npn(DATA_combined, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
data_npn <- as.data.frame(apply(data_npn, use='pairwise.complete.obs', 2, as.numeric))
summary(data_npn)

#### GGM & NODE PREDICTABILITY 

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

## Estimate the network

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

## Lasso analysis 

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


## Bridge centrality analysis 

adjm <- getWmat(NET)
matrix<- graph_from_adjacency_matrix(adjm, mode="upper", diag=FALSE, weighted=TRUE)


## "walktrap" Community Analysis

walktrap_community <- cluster_walktrap(matrix)
# head(walktrap_community)
plot_dendrogram(walktrap_community, mode="phylo")

# Bridge expected influence of the Gaussian graphic model

my_communities <- walktrap_community$membership
checkVector <- is.vector(my_communities)
b <-bridge(matrix,communities=my_communities)
BPLOT <- plot(b, 
              include=c("Bridge Expected Influence (1-step)"),
              theme_bw=FALSE, raw0 = TRUE, signed=TRUE, decreasing= TRUE)
BPLOT + theme(axis.text=element_text(size=6.5))

## Stability and accuracy of the GGM

# This takes ages to run so leave it out for the quick demo version! 
nonParametricBoot <- bootnet(data_npn, nBoots = 1000, 
                             default = "ggmModSelect", 
                             type="nonparametric", 
                             statistics="all", communities=my_communities)
plot(nonParametricBoot, labels=FALSE, order="sample")
plot(nonParametricBoot, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
summary(nonParametricBoot) %>% ungroup %>% filter(type == "edge") %>% arrange(-sample) %>% View
plot(nonParametricBoot, statistics="ExpectedInfluence", plot="difference",order = "mean")
plot(nonParametricBoot, statistics="bridgeExpectedInfluence", plot="difference",order = "mean")


#stability of the expected influence centrality - this also takes ages

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



# AIM 2: Estimate Bayesian network -----------------------
## Fit a first Bayesian network, based on 50 random re-starts and 100 perturbations for each re-start.
netdata <- as.data.frame(data_npn) # DA: Need to turn it back into a data frame for this analsyis
set.seed(123)
fitBN1 <- hc(netdata, restart = 50, perturb = 100)  ## hc gives directed graph

#fitBN1

bnlearn::score(fitBN1, data = netdata)         
astr <- arc.strength(fitBN1, netdata, "bic-g")  ## connection strength
astr[order(astr[,3]), ]  ## sorted edge strength from strongest to weakest

strength.plot(fitBN1, astr, shape = "ellipse")

## Now we stabilize the network across multiple samples through bootstrapping:
## Learn 10000 network structures (may take a few min)

# DA note: This was originally 012 as the seed - this produces a completely different result 
## to the paper, with "worry" being outside the complex with no arrows. 
## However, setting the seed as 321 produced an identical result to 123. Go figure.


set.seed(123)
bootnet <- boot.strength(netdata, R = 10000, algorithm = "hc", 
                         algorithm.args = list(restart = 5, perturb = 10), debug = TRUE)  
head(bootnet)


# Save the bootstrapped networks
saveRDS(bootnet, "./DASStotal_bootDAG.rds")

# Reload the saved bootstrapped networks
#bootnet <- readRDS("./Anxiety_bootDAG.rds")

## net#1 optimal cutpoint, according to Scurati & Nagarajan (2013)
avgnet <- averaged.network(bootnet) 
# avgnet      #This is currently not working!  Not sure why, the rest of it works! 
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

## thick arrows indicate high directional probabilities, thin arrows low directional probabilities
# FOR interpretation,   Edge thickness signifies the probability of prediction is in the direction depicted.


