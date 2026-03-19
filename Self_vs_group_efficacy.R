#------------PACKAGES ----------------

# To do: Investigate whether all of these packages are really used. 
require(huge)
require(networktools)
require(qgraph)
require(parcor) # This one is kind of a pain because the package is deprecated. Check if it's really required? 
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
require(Rgraphviz)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rgraphviz")

require(tcltk)
require(latentnet)
require(rgl)

library(haven)

## NEXT NETWORK MODEL FOR AMY
## CCAS cog-emotional, CCAS functional, PEB total, Self vs group efficacy and EAI_EMA_total

youth_climate_data <- read_sav("Youth Mental Health and Climate Change_2023_IMPUTED_FOR ANALYSIS_wID_877_01_02_2026.sav")

Cognitive_Emo <-  youth_climate_data$CCAS_CE_impairment
Functional <-  youth_climate_data$CCAS_functional_impair

PersonalBeh <- youth_climate_data$PEB_total
Activism <- youth_climate_data$EAI_EMA_total

SelfEfficacy <- youth_climate_data$Self_efficacy_total
GroupEfficacy <- youth_climate_data$Group_efficacy_total


DATA_combined_social <- as.data.frame(cbind(Cognitive_Emo,Functional,
                                            PersonalBeh, Activism,
                                            SelfEfficacy, GroupEfficacy))

summary(DATA_combined_social)

M <- cor(DATA_combined_social)
corrplot(M, method="shade") #Figure S1

is.positive.definite(cor(DATA_combined_social))#true 

## Use goldbricker function to search for potential "bad pairs"
## This function compares correlations in a psychometric network in order to identify
## nodes which most likely measure the same underlying construct (i.e., are colinear).
## A "bad pairs" includes two nodes that are highly correlated with each other and that 
## correlate with other nodes in the network in a highly similar pattern.
gb <- goldbricker(DATA_combined_social, p = 0.05, 
                  method = "hittner2003", threshold=0.25, 
                  corMin=.50, progressbar = TRUE)
## Produce list of "bad pairs"
gb

#transform data with nonparanormal transform

data_npn_social <- huge.npn(DATA_combined_social, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
data_npn_social <- as.data.frame(apply(data_npn_social, use='pairwise.complete.obs', 2, as.numeric))
summary(data_npn_social)

#### GGM & NODE PREDICTABILITY 

group.ega<- c("Climate Anxiety", "Climate Anxiety", 
              "Behaviour", "Behaviour", "Efficacy", "Efficacy") 
p <- ncol(data_npn_social)

data_npn_social <- as.matrix(data_npn_social) # DA: I had to add this line to avoid an error 

fit <- mgm(data_npn_social, 
           type = rep('g', p),
           lev = rep(1, p),
           rule.reg = 'OR')

pred <- predict(fit, data_npn_social, 
                errorCon = 'R2')
pred$error

## Estimate the network

NET <- estimateNetwork(data_npn_social, default="ggmModSelect", tuning = 0) 
plot(NET,node.width=1.4,layout = 'spring', 
     label.prop = 0.9, borders = TRUE, 
     pie = pred$error[,2], pieColor = rep('#377EB8',p),
     labels = colnames(data_npn_social),
     edge.labels=FALSE, details=FALSE, groups = group.ega,color= c("grey", "grey", 
                                                                   "yellow","yellow",
                                                                   "red", "red"), 
     vsize=12, border.width=1.9, label.cex = .9,
     label.scale.equal=TRUE, 
     color='#bbbbbb', border.width=6, legend = FALSE)

centralityPlot(NET, theme_bw=TRUE, scale = "raw0", 
               include = c("ExpectedInfluence"), decreasing= TRUE)

## Lasso analysis 

Network_social <- EBICglasso(cor(data_npn_social), 
                             n = nrow(data_npn_social), 
                             gamma = 0.5, threshold = TRUE) 

y <- qgraph(Network_social, layout = 'spring', 
            details = FALSE, labels=colnames(data_npn_social),
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
# This currently does not work because of negative weights. 
# Advice from the internet: If your weights represent similarities and negative weights represent dissimilarity,
# the easiest is just to drop the edges with negative weights -- 
# the absence of an edge should be enough indication for walktrap 
# to consider two nodes as dissimilar.

# matrix <- matrix %>%
 # delete_edges("Functional|Activism")

walktrap_community <- cluster_walktrap(matrix)
# head(walktrap_community)
plot_dendrogram(walktrap_community, mode="phylo")

# Bridge expected influence of the Gaussian graphic model

my_communities <- walktrap_community$membership
checkVector <- is.vector(my_communities)
checkVector
b <-bridge(matrix,communities=my_communities)
BPLOT <- plot(b, 
              include=c("Bridge Expected Influence (1-step)"),
              theme_bw=FALSE, raw0 = TRUE, signed=TRUE, decreasing= TRUE)
BPLOT + theme(axis.text=element_text(size=6.5))

# This takes ages to run so leave it out for the quick demo version! 
nonParametricBoot <- bootnet(data_npn_social, nBoots = 1000, 
                             default = "ggmModSelect", 
                             type="nonparametric", 
                             statistics="all", communities=my_communities)
plot(nonParametricBoot, labels=FALSE, order="sample")
plot(nonParametricBoot, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
summary(nonParametricBoot) %>% ungroup %>% filter(type == "edge") %>% arrange(-sample) %>% View
plot(nonParametricBoot, statistics="ExpectedInfluence", plot="difference",order = "mean")
plot(nonParametricBoot, statistics="bridgeExpectedInfluence", plot="difference",order = "mean")

#stability of the expected influence centrality - this also takes ages
# In the initial model this seems to throw a lot of messages about the assumption of sparsity being violated

caseDroppingBoot <- bootnet(data_npn_social, nBoots = 1000, 
                            default = "EBICglasso",  type="case", 
                            statistics="all", communities=my_communities)
corStability(caseDroppingBoot)
plot(caseDroppingBoot, statistics = c("ExpectedInfluence","bridgeExpectedInfluence"))

### DIFFERENTIAL VARIABILITY OF THE GGM
#  restricted range problem
c<- centrality_auto(matrix)
x <- c$node.centrality$ExpectedInfluence 
SD <- apply(data_npn_social,2,sd)
cor(SD,x)
cor.test(SD,x) 


# AIM 2: Estimate Bayesian network -----------------------
## Fit a first Bayesian network, based on 50 random re-starts and 100 perturbations for each re-start.
netdata <- as.data.frame(data_npn_social) # DA: Need to turn it back into a data frame for this analsyis
set.seed(321)
fitBN1 <- hc(netdata, restart = 50, perturb = 100)  ## hc gives directed graph
fitBN1

## Now we stabilize the network across multiple samples through bootstrapping:
## Learn 10000 network structures (may take a few min).


set.seed(42)
bootnet <- boot.strength(netdata, R = 10000, algorithm = "hc", 
                         algorithm.args = list(restart = 5, perturb = 10), debug = TRUE)  
head(bootnet)

# Save the bootstrapped networks
saveRDS(bootnet, "./social_bootDAG.rds")

# Reload the saved bootstrapped networks
#bootnet <- readRDS("./social_bootDAG.rds")

## net#1 optimal cutpoint, according to Scurati & Nagarajan (2013)
avgnet <- averaged.network(bootnet)
#avgnet      
bnlearn::score(avgnet, data = netdata) # -6474.202
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







