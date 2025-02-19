#!/usr/bin/env Rscript

###############################################################################
## calculate_tree_distance_metrics.R
## v0.1
## 
## Calculate tree distance metrics between two trees, and additional custom
## per-tree metrics.
##
## Author: Isin Altinkaya
## Date Created: 2025-02-17
## Email: isinaltinkaya@gmail.com
###############################################################################
##
## Notes:
##
###############################################################################
# -> Load libraries
library(ape)
library(TreeDist)
library(optparse)
library(stringr)

#library(phangorn)
#library(phytools)
#library(dynamicTreeCut)

###############################################################################
# -> Define generic functions

assert <- function(condition, error_message = "Assertion failed", success_message = NULL) { 
    condition_result <- tryCatch(condition, error = function(e) e)
    
    if (!is.logical(condition_result)) {
        stop("\n\n##########\n# ERROR: The condition does not return a logical value.\n##########\n", call. = FALSE)
    }
    
    if (condition_result != TRUE) {
        calls <- sys.calls()
        call <- calls[1:(length(calls) - 1)]
        error_message <- gsub("\n$", "", error_message)
        stop(paste0("\n\n##########\n# ERROR: ", error_message,
                    "\n# (", call, ") is FALSE",
                    "\n##########\n"), call. = FALSE)
    }
    
    if (!is.null(success_message)) {
        cat(success_message)
    }
}

extract_pop<-function(str){
	sapply(str_split(str, "_"), function(x) x[1])
}

###############################################################################
# -> Read input arguments

option_list <- list(
  make_option(c("-tt", "--truthtree"), type="character", default=NULL, 
			  help="Path to the first tree file"),
  make_option(c("-et", "--esttree"), type="character", default=NULL,
			  help="Path to the second tree file"),
  make_option(c("-o", "--output"), type="character", default=NULL,
			  help="Path to the output file")
)

parser <- OptionParser(option_list=option_list, usage="%prog --truthtree <truthtree> --esttree <tree2> --output <output>")
args <- parse_args(parser)

ttfile <- args$truthtree
etfile <- args$esttree
ofile <- args$output

cat("\n-> Calculating tree distance metrics between the truth tree at ", ttfile, " and the estimated tree at ", etfile, "\n")

###############################################################################
# -> Read trees and perform sanity checks

# test
#ttfile="tools/ngsAMOVA/utils/sim_amova_2412-model1-c1-csrep0.nwk"
#etfile="tools/ngsAMOVA/utils/sim_amova_2412-model1-c1-csrep1.nwk"
##etfile="tools/ngsAMOVA/utils/sim_amova_2412-model1-merged5-csrep0-d100-e0.002-glrep0.truth.nwk"

# Read trees
tt <- read.tree(ttfile)
et <- read.tree(etfile)

# -> Arbitrary nodes are chosen as root for unrooted tree output from the NJ in ngsAMOVA
# So unroot the trees before calculating the metrics
tt <- unroot(tt)
et <- unroot(et)

assert(is.rooted(tt)==is.rooted(et))
is_rooted=is.rooted(tt)

# -> Check if the trees are binary
# As sanity check
assert(is.binary(tt)==is.binary(et))
assert(is.binary(tt)==TRUE, "The trees are not binary")


###############################################################################
# -> Get ratio of misplaced populations in the subtrees
# Calculate how many subtrees have misplaced populations (i.e. populations from different true populations are in the same subtree, if the subtree has less than 10 individuals as each population contains 10 individuals)
get_ratio_misplaced_pops_subtrees <- function(tree){
	l<-subtrees(tree)
	cc<-c()
	for (i in 1:length(l)){
		# bc each population contains 10 individuals
		if(length(unique(l[[i]]$node.label))<=10){
			# check if all individuals in the subtree are from the same population
			cc<-c(cc,(length(unique(extract_pop(l[[i]]$tip.label)))==1))
		}
	}
	return(1-(sum(cc)/length(cc)))
}

###############################################################################
# -> Get ratio of misplaced populations in the clusters from hierarchical clustering
# Calculate how many clusters have misplaced populations (i.e. populations from different true populations are in the same cluster)
get_ratio_misplaced_pops_clusters<-function(tree, n_clusters=4){
	dist_matrix <- cophenetic.phylo(tree)
	# perform hierarchical clustering on the distance matrix
	hc <- hclust(as.dist(dist_matrix))
	# cut the tree into n_clusters clusters 
	clusters <- cutree(hc, k = n_clusters)
	names(clusters)<-sub("_.*","",names(clusters))
	uniq_clusters<-c()
	for(i in unique(clusters)){
		uniq_clusters<-c(uniq_clusters,length(unique(names(clusters)[clusters==i]))==1)
	}
	return(1-(sum(uniq_clusters)/length(uniq_clusters)))
}


###############################################################################
# -> Get ratio of non-monophyletic populations in the tree
# Calculate how many populations are not monophyletic in the tree

get_ratio_non_monophyletic_pops<-function(tree, pops=c("popA","popB","popC","popD")){
	c_pop_is_monophyletic<-c()
	for(pop in pops){
		pop_tips <- tree$tip.label[grep(pop, tree$tip.label)]
		c_pop_is_monophyletic<-c(c_pop_is_monophyletic,is.monophyletic(tree, tip = pop_tips))
	}
	return(1-(sum(c_pop_is_monophyletic)/length(c_pop_is_monophyletic)))
}

###############################################################################

###############################################################################

#clusters <- cutree(hc, k = 4)
#tip<-NULL
#for(g in 1:length(unique(clusters))){
#	individuals<-names(clusters)[clusters==g]
#	populations<-sub("_.*","",individuals)
#	unique_individuals<-individuals[!duplicated(populations)]
#	tip<-c(tip,unique_individuals)
#}

## Create a new tree containing only the representative tips
#new_tree <- keep.tip(tree, tip)

## Plot the new tree of groups
#plot(new_tree, main = "Tree of Groups with Representative Individuals as Tips")


###############################################################################
# -> Others:

#ClusteringInfoDistance(tt,et,reportMatching=TRUE)
#ExpectedVariation(tt,et)

##Nearest Neighbour Interchange (NNI) Distance
#NNIDist(tt,et)

#library(Quartet)
#qd <- as.dist ( Quartet::QuartetDivergence(
#	Quartet::ManyToManyQuartetAgreement(list(tt, et)), similarity = FALSE
#))
###############################################################################



###############################################################################
# -> Calculate tree distance metrics

metrics <- data.frame(
##Information‐Based Generalized Robinson–Foulds Distances
	"tree_distance" = TreeDistance(tt,et),

## -> same
	"diff_phylogenetic_info" = DifferentPhylogeneticInfo(tt,et),
	#"phylogenetic_info_distance" = PhylogeneticInfoDistance(tt,et),
## <-

	"clustering_info_distance" = ClusteringInfoDistance(tt,et,reportMatching=FALSE),
	"matching_split_info_distance" = MatchingSplitInfoDistance(tt,et),

##Nye et al. (2006) Tree Comparison
	"nye_similarity" = NyeSimilarity(tt,et),

##Kendall–Colijn Distance (Incorporates Branch Lengths)
	"kendall_colijn" = KendallColijn(tt,et),

## -> same
	"branch_score_difference" = treedist(tt, et, check.labels = TRUE)["branch.score.difference"][[1]],
	#"KF_dist" = KF.dist(tt, et , check.labels = TRUE, rooted=is_rooted),
## <-

	"quadratic_path_difference" = treedist(tt, et, check.labels = TRUE)["quadratic.path.difference"][[1]],

## -> same
	"symmetric_difference" = treedist(tt, et, check.labels = TRUE)["symmetric.difference"][[1]],
	#"RF_dist" = RF.dist(tt, et , normalize = FALSE, check.labels = TRUE, rooted=is_rooted),
	#"sprdist_rf" = sprdist(tt, et)["rf"][[1]],
## <-

# N.B. is_rooted==TRUE gives "Some trees are not binary. Result may not what you expect!"
	#"wRF_dist" = wRF.dist(tt, et , normalize = FALSE, check.labels = TRUE, rooted=is_rooted),

## -> same
##Path Distance
	#"path_difference" = treedist(tt, et, check.labels = TRUE)["path.difference"][[1]],
	#"pathDist" = PathDist(tt,et),
	"path_dist" = path.dist(tt, et , check.labels = TRUE, use.weight = FALSE),
## <-

##Matching Split Distance
	"matching_split_distance" = MatchingSplitDistance(tt,et),

##Jaccard–Robinson–Foulds Metric
	"jaccard_robinson_foulds" = JaccardRobinsonFoulds(tt,et),

##Maximum Agreement Subtree (MAST) Size
# N.B. is_rooted==TRUE gives "Error in Func(tree1, tree2, tipLabels = labels1, nTip = nTip, ...) : Both trees must be rooted if rooted = TRUE"
	#"MAST_size" = MASTSize(tt,et,rooted=is_rooted),

# N.B. is_rooted==TRUE gives "Error in Func(tree1, tree2, tipLabels = labels1, nTip = nTip, ...) : Both trees must be rooted if rooted = TRUE"
	#"MAST_info" = MASTInfo(tt,et,rooted=is_rooted),

##Robinson–Foulds Distances (with Phylogenetic Information Adjustments)
	"info_robinson_foulds" = InfoRobinsonFoulds(tt,et),

##Subtree Prune and Regraft (SPR) Distance
	"SPR_dist" = SPRDist(tt,et),
	#"sprdist" = sprdist(tt, et),
	"sprdist_spr" = sprdist(tt, et)["spr"][[1]],
	"sprdist_spr_extra" = sprdist(tt, et)["spr_extra"][[1]],
	"sprdist_hdist" = sprdist(tt, et)["hdist"][[1]],
	"quartet_divergence" = as.dist(Quartet::QuartetDivergence(Quartet::ManyToManyQuartetAgreement(list(tt, et)), similarity = FALSE))[[1]],

##Custom Metrics
	"truth_ratio_misplaced_pops_subtrees" = get_ratio_misplaced_pops_subtrees(tt),
	"estimated_ratio_misplaced_pops_subtrees" = get_ratio_misplaced_pops_subtrees(et),
	"truth_ratio_misplaced_pops_clusters" = get_ratio_misplaced_pops_clusters(tt),
	"estimated_ratio_misplaced_pops_clusters" = get_ratio_misplaced_pops_clusters(et),
	"truth_ratio_non_monophyletic_pops" = get_ratio_non_monophyletic_pops(tt),
	"estimated_ratio_non_monophyletic_pops" = get_ratio_non_monophyletic_pops(et)
)

# -> Duplicate metrics to exclude
#metrics$diff_phylogenetic_info==metrics$phylogenetic_info_distance
#metrics$KF_dist==metrics$branch_score_difference
#metrics$symmetric_difference==metrics$RF_dist
#metrics$RF_dist==metrics$sprdist_rf
#metrics$path_difference==metrics$path_dist

###############################################################################


# -> Write output to file
cat("\n-> Writing output to file: ", ofile, "\n")
write.table(metrics, file=ofile, sep="\t", row.names=FALSE, col.names=FALSE)
