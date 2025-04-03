## code to prepare `DATASET` dataset goes here

library(readxl)
library(ape)
library(phangorn)
library(iNEXT.3D)
library(iNEXT.beta3D)

otu_file <- "OTU_tables.xlsx"


sheet_names <- excel_sheets(otu_file)

otu_data <- read_excel(otu_file, sheet = sheet_names[2])

head(otu_data)

otu_ids <- unique(as.character(otu_data$OTU))


tree_file <- "tree_data_Figure2.tre"
phylo_tree <- read.tree(tree_file)


tip_labels <- phylo_tree$tip.label


extract_middle_number <- function(label) {
  matches <- regmatches(label, regexpr("[0-9]{5}", label))  # 提取 5 位數字
  return(ifelse(length(matches) > 0, paste0("OTU_", matches), NA))  # 加上 "OTU_"
}


tip_numbers <- sapply(tip_labels, extract_middle_number, USE.NAMES = FALSE)




matching_tips <- tip_labels[tip_numbers %in% otu_ids]


cat("匹配的 Tip Labels 數量：", length(matching_tips), "\n")


pruned_tree <- keep.tip(phylo_tree, matching_tips)


duplicated_tips <- pruned_tree$tip.label[duplicated(pruned_tree$tip.label)]
if (length(duplicated_tips) > 0) {
  pruned_tree$tip.label <- make.unique(pruned_tree$tip.label)  # 使 Tip Labels 唯一
}




if (!is.null(pruned_tree$node.label)) {
  pruned_tree$node.label <- paste0("NODE_", seq_along(pruned_tree$node.label))
} else {

  n_internal <- pruned_tree$Nnode
  pruned_tree$node.label <- paste0("NODE_", seq_len(n_internal))
}


# pruned_tree <- multi2di(pruned_tree, random=TRUE)  # 轉換為二分樹
# pruned_tree <- midpoint.root(pruned_tree)  # 進行 Midpoint rooting
# pruned_tree <- compute.brlen(pruned_tree, method="Grafen")  # 計算分支長度
# pruned_tree <-compute.brlen(pruned_tree,power =0.5)
tree_clean <- di2multi(pruned_tree, tol=1e-8)
pruned_tree <- chronos(tree_clean,lambda = 10)
plot(pruned_tree, main="Pruned Phylogenetic Tree", edge.width=2, show.tip.label=FALSE)

# T03_tree = pruned_tree
# T09_tree = pruned_tree


data <- otu_data
data <- as.data.frame(data)


rownames(data) <- data$OTU
data <- data[,-1]
data


iNEXT.3D::DataInfo3D(data)
# iNEXT.beta3D::DataInfobeta3D(data,diversity = "PD",PDtree = pruned_tree)


# T03T04
T03T04 <- data[, c(1,7)]
colnames(T03T04) <- c("T03", "T04")

# T09T10
T09T10 <- data[, c(6, 8)]
colnames(T09T10) <- c("T09", "T10")



extract_middle_number <- function(label) {
  matches <- regmatches(label, regexpr("[0-9]{5}", label))  # 提取 5 位數字
  return(ifelse(length(matches) > 0, paste0("OTU_", matches), label))  # 如果沒匹配，保留原標籤
}


tree_labels <- sapply(pruned_tree$tip.label, extract_middle_number, USE.NAMES = FALSE)

pruned_tree$tip.label <- tree_labels


pruned_tree$tip.label <- make.unique(pruned_tree$tip.label)

head(pruned_tree$tip.label)
plot(pruned_tree, main="Updated Phylogenetic Tree", edge.width=2, show.tip.label=FALSE)

otu_ids <- rownames(T03T04)
tree_labels <- pruned_tree$tip.label
T03T04 <- T03T04[otu_ids %in% tree_labels, , drop = FALSE]

otu_ids <- rownames(T09T10)
tree_labels <- pruned_tree$tip.label
T09T10 <- T09T10[otu_ids %in% tree_labels, , drop = FALSE]

fungi = list(T03T04 = T03T04, T09T10 = T09T10)

# save(T03T04, T09T10, pruned_tree, file = "data.RData")
fungi_tree <- as.phylo(pruned_tree)
class(fungi_tree) <-"phylo"

usethis::use_data(fungi,fungi_tree, overwrite = TRUE)




library(poppr)
library(tidyverse)
library(ape)
library(phytools)
library(phyclust)
library(iNEXT.seq)
library(iNEXT.beta3D)
library(ade4)
library(phyloseq)

data(GlobalPatterns)

otu_table <- otu_table(GlobalPatterns)
iNEXT.3D::DataInfo3D(otu_table)
iNEXT.beta3D::DataInfobeta3D(otu_table)
sample_data <- sample_data(GlobalPatterns)
phylo_tree_ori <- phy_tree(GlobalPatterns)
# plot(phylo_tree_ori, direction = "downwards")
# tax_table <- tax_table(GlobalPatterns)["951",]
fresh_ocean_samples <- phyloseq::subset_samples(GlobalPatterns, SampleType == "Freshwater" | SampleType == "Ocean" | SampleType == "Sediment (estuary)")
# ocean_samples <- subset_samples(GlobalPatterns, SampleType == "Ocean")
otu_table <- otu_table(fresh_ocean_samples)
sample_table <- sample_data(fresh_ocean_samples)
tax_table <- tax_table(fresh_ocean_samples)
phylo_tree_ori <- phy_tree(fresh_ocean_samples)
phylo_tree_ori$node.label <- NULL


mat <- matrix(data = 0, nrow = 3, ncol = 7)
mat[1,] = "All"
mat[2,1:3] = "Sediment"
mat[2,4:5] = "Ocean"
mat[2,6:7] = "Freshwater"
mat[3,] = c("TRRsed1", "TRRsed2","TRRsed3" ,"NP2", "NP3 ","SLEpi20M","LMEpi24M")

#TRR sed1 河口沉積物 AQC1cm 淡水SLEpi20M 淡水 LMEpi24M 淡水
dat_Type <- list(TRRsed1 = otu_table[,sample_table$X.SampleID == "TRRsed1"],
                 TRRsed2 = otu_table[,sample_table$X.SampleID == "TRRsed2"],
                 TRRsed3 = otu_table[,sample_table$X.SampleID == "TRRsed3"],
                 NP2 = otu_table[,sample_table$X.SampleID == "NP2"],
                 NP3 = otu_table[,sample_table$X.SampleID == "NP3"],
                 # NP5 = otu_table[,sample_table$X.SampleID == "NP5"],
                 SLEpi20M = otu_table[,sample_table$X.SampleID == "SLEpi20M"],
                 LMEpi24M = otu_table[,sample_table$X.SampleID == "LMEpi24M"])



dat <- sapply(dat_Type, function(x) rowSums(x))
dat_check <- dat[rowSums(dat)>0,]
iNEXT.3D::DataInfo3D(dat_check)
iNEXT.beta3D::DataInfobeta3D(dat_check)




#第一種資料縮減模式
# 去除所有行的總和為 0 的行
dat <- sapply(dat_Type, function(x) rowSums(x))

dat_nonzero <- dat[rowSums(dat)>0,]

#只取5%
set.seed(112024510)
keep_rows <- sample(1:nrow(dat_nonzero), size = round(nrow(dat_nonzero) / 20))
dat_nonzero <- dat_nonzero[keep_rows, ]

rownames(dat_nonzero) <- paste0("OTU_", rownames(dat_nonzero))
iNEXT.3D::DataInfo3D(dat_nonzero)
iNEXT.beta3D::DataInfobeta3D(dat_nonzero)
phylo_tree_ori$tip.label <- paste0("OTU_", phylo_tree_ori$tip.label)
absent <- phylo_tree_ori$tip.label[!(phylo_tree_ori$tip.label %in% rownames(dat_nonzero))]
# absent <- phylo_tree_ori$tip.label[!(phylo_tree_ori$tip.label %in% rownames(dat))]
if (length(absent) != 0) {
  phylo_tree <- drop.tip(phylo_tree_ori, absent)
}

plot(phylo_tree, direction = "downwards", show.tip.label = F)
phylo_tree <- multi2di(phylo_tree, control = list(tol = 1e-10))
# 將 phylo_tree 計算成 ultrametric tree
phylo_tree <- compute.brlen(phylo_tree, method = "Grafen")
# tree = phylo_tree %>% chronos() #%>% multi2di() %>% compute.brlen()
tree =phylo_tree
# tree = phylo_tree 
class(tree) = "phylo"
# write.tree(tree, file = "global_tree_I.txt") #存到 Unifrac那個資料夾
# tree = read.tree("global_tree_I.txt")  #且這個樹是用chrono用的
# plot(phylo_tree, direction = "downwards", show.tip.label = F)
plot(tree,direction = "downwards", show.tip.label = F)

global_tree = tree
global_mat = mat
global <- as.data.frame(dat_nonzero)
# save(global, global_tree,global_mat, file = "Global_data.RData")
usethis::use_data(global, global_tree,global_mat, overwrite = TRUE)



