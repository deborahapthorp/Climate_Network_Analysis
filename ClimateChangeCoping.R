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

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Rgraphviz")

require(Rgraphviz)

require(tcltk)
require(latentnet)
require(rgl)

library(haven)

## NEXT NETWORK MODEL FOR AMY

## CCAS cog-emotional, CCAS functional, PEB total, Coping - CCCoping_meaning, CCCoping_de_emphasizing, CCCoping_problem


youth_climate_data <- read_sav("Youth Mental Health and Climate Change_2023_IMPUTED_FOR ANALYSIS_wID_877_01_02_2026.sav")

youth_climate_data<-na.omit(youth_climate_data) # remove NAs

Cognitive_Emo <-  youth_climate_data$CCAS_CE_impairment
Functional <-  youth_climate_data$CCAS_functional_impair

PEB <- youth_climate_data$PEB_total
Coping_meaning <- youth_climate_data$CCCoping_meaning
Coping_de_emph <- youth_climate_data$CCCoping_de_emphasizing
Coping_problem <- youth_climate_data$CCCoping_problem

DATA_combined_coping <- as.data.frame(cbind(Cognitive_Emo,
                                            Functional,
                                            PEB,
                                            Coping_meaning, Coping_problem, Coping_de_emph))

summary(DATA_combined_coping)

M <- cor(DATA_combined_coping)
corrplot(M, method="shade") #Figure S1

check <- is.positive.definite(cor(DATA_combined_coping))#true 

check

## Use goldbricker function to search for potential "bad pairs"
## This function compares correlations in a psychometric network in order to identify
## nodes which most likely measure the same underlying construct (i.e., are colinear).
## A "bad pairs" includes two nodes that are highly correlated with each other and that 
## correlate with other nodes in the network in a highly similar pattern.
gb <- goldbricker(DATA_combined_coping, p = 0.05, 
                  method = "hittner2003", threshold=0.25, 
                  corMin=.50, progressbar = TRUE)
## Produce list of "bad pairs"
gb


#transform data with nonparanormal transform

data_npn_coping <- huge.npn(DATA_combined_coping, npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
data_npn_coping <- as.data.frame(apply(data_npn_coping, use='pairwise.complete.obs', 2, as.numeric))
summary(data_npn_coping)



#### GGM & NODE PREDICTABILITY 

group.ega<- c("Climate Anxiety", "Climate Anxiety", 
              "Behaviour", "Coping", "Coping", "Coping") 
p <- ncol(data_npn_coping)

data_npn_coping <- as.matrix(data_npn_coping) # DA: I had to add this line to avoid an error 

fit <- mgm(data_npn_coping, 
           type = rep('g', p),
           lev = rep(1, p),
           rule.reg = 'OR')

pred <- predict(fit, data_npn_coping, 
                errorCon = 'R2')
pred$error

## Estimate the network

NET <- estimateNetwork(data_npn_coping, default="ggmModSelect", tuning = 0) 
plot(NET,node.width=1.4,layout = 'spring', 
     label.prop = 0.9, borders = TRUE, 
     pie = pred$error[,2], pieColor = rep('#377EB8',p),
     labels = colnames(data_npn_coping),
     edge.labels=FALSE, details=FALSE, groups = group.ega,color= c("grey", "grey", 
                                                                   "yellow", "red", "red", "red"), 
     vsize=12, border.width=1.9, label.cex = .9,
     label.scale.equal=TRUE, 
     color='#bbbbbb', border.width=6, legend = FALSE)

centralityPlot(NET, theme_bw=TRUE, scale = "raw0", 
               include = c("ExpectedInfluence"), decreasing= TRUE)

## Lasso analysis 

Network_coping <- EBICglasso(cor(data_npn_coping), 
                             n = nrow(data_npn_coping), 
                             gamma = 0.5, threshold = TRUE) 

y <- qgraph(Network_coping, layout = 'spring', 
            details = FALSE, labels=colnames(data_npn_coping),
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

# need to delete negative edges for the walktrap communuty analysis to work. 
matrix <- delete_edges(matrix, "Coping_problem|Coping_de_emph") 
matrix <- delete_edges(matrix, "PEB|Coping_de_emph")

## "walktrap" Community Analysis

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
nonParametricBoot <- bootnet(data_npn_coping, nBoots = 1000, 
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
SD <- apply(data_npn_coping,2,sd)
cor(SD,x)
cor.test(SD,x) 

# AIM 2: Estimate Bayesian network -----------------------
## Fit a first Bayesian network, based on 50 random re-starts and 100 perturbations for each re-start.
netdata <- as.data.frame(data_npn_coping) # DA: Need to turn it back into a data frame for this analsyis
set.seed(123)
fitBN1 <- hc(netdata, restart = 50, perturb = 100)  ## hc gives directed graph

bnlearn::score(fitBN1, data = netdata)         
astr <- arc.strength(fitBN1, netdata, "bic-g")  ## connection strength
astr[order(astr[,3]), ]  ## sorted edge strength from strongest to weakest

strength.plot(fitBN1, astr, shape = "ellipse")

## Now we stabilize the network across multiple samples through bootstrapping:
## Learn 10000 network structures (may take a few min).


set.seed(321)
bootnet_coping <- boot.strength(netdata, R = 10000, algorithm = "hc", 
                         algorithm.args = list(restart = 5, perturb = 10), debug = TRUE)  
head(bootnet_coping)




# Save the bootstrapped networks
# saveRDS(bootnet_coping, "./coping_bootDAG.rds")

# Reload the saved bootstrapped networks
#bootnet_coping <- readRDS("./coping_bootDAG.rds")

## net#1 optimal cutpoint, according to Scurati & Nagarajan (2013)
avgnet <- averaged.network(bootnet_coping)
avgnet      
bnlearn::score(avgnet, data = netdata) # -5126.457
thresh <- avgnet$learning$args[[1]]
thresh                                  ## optimal significance threshold = 0.3319
astr <- arc.strength(avgnet, netdata, "bic-g")   ## compute edge strengths
astr

suppressWarnings(strength.plot(avgnet, astr, shape = "ellipse"))
# FOR interpretation: Edge thickness signifies the probability of prediction is in the direction depicted.

## CCAS cog-emotional, CCAS functional, PEB total, Self vs group efficacy and EAI_EMA_total