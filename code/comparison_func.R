
library(ggplot2)

# colors
pal <- c("#ffb6db", "#33A02C", "#b66dff", "#FEC44F", "#41B6C4", "#8E0152", "#0868AC", "#807DBA", "#E7298A", "#810F7C", 
         "#00441B", "#525252", "#4D9221", "#8B5742", "#D8DAEB", "#7cdd2d", "#980043", "#8C96C6", "#EC7014", "#92C5DE", 
         "#FDAE61", "#1D91C0", "#A6DBA0", "#4292C6", "#BF812D", "#01665E", "#41AB5D", "#FE9929", "#252525", "#A6761D")
names(pal) <- 1:30

pal2 <- c("#7CAE00", "#F8766D", "#00BFC4")
names(pal2) <- c("Common", "A Enriched", "B Enriched")

# function to visualize clusters for all cells
plotcluster_func <- function(gene_umap, cluster, method, filename) {
  df <- data.frame(UMAP1 = gene_umap[,1], UMAP2 = gene_umap[,2], 
                   Cluster = as.character(cluster))
  df_naomit <- na.omit(df)
  
  ggplot(df_naomit, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 0.5) +
    scale_color_manual(values=c(pal)) +
    guides(colour=guide_legend(override.aes=list(size = 5))) +
    ggtitle(method)
  ggsave(filename,height = 3, width = 3.5)
}

# # function to visualize cluster types for all cells
# plotclustertype_func <- function(gene_umap, cluster_type, method, filename) {
#   cluster_type_temp <- rep(NA, length(cluster_type))
#   cluster_type_temp[which(cluster_type == "C")] <- rep("Common", length(which(cluster_type == "C")))
#   cluster_type_temp[which(cluster_type == "A")] <- rep("A Enriched", length(which(cluster_type == "A")))
#   cluster_type_temp[which(cluster_type == "B")] <- rep("B Enriched", length(which(cluster_type == "B")))
#   cluster_type_temp <- factor(cluster_type_temp, levels=c("Common", "A Enriched", "B Enriched"))
#   df <- data.frame(UMAP1 = gene_umap[,1], UMAP2 = gene_umap[,2], 
#                    Cluster_type = as.character(cluster_type_temp))
#   df_naomit <- na.omit(df)
#   
#   ggplot(df_naomit, aes(x = UMAP1, y = UMAP2, color = Cluster_type)) +
#     geom_point(size = 0.3) +
#     scale_color_manual(values=c(pal2)) +
#     guides(colour=guide_legend(override.aes=list(size = 5))) +
#     ggtitle(method)
#   ggsave(filename,height = 3, width = 3.5)
# }

# function to calculate num of cells from two conditions in each cluster
cluster_num_func <- function(cluster, condition) {
  ncluster_temp <- length(unique(cluster))
  cluster_num <- matrix(0, ncol=ncluster_temp, nrow=2)
  cluster_num[1, ] <- table(cluster[condition == "A"])
  cluster_num[2, ] <- table(cluster[condition == "B"])
  colnames(cluster_num) <- levels(cluster)
  rownames(cluster_num) <- c("A", "B")
  return(cluster_num)
}

# function to calculate cluster type
cluster_type_func <- function(cluster, cluster_num, wdiff_tol, alpha_ci = 0.05) {
  z <- qnorm(alpha_ci/2)
  th_star <- apply(cluster_num, 1, sum)
  th_star <- th_star/sum(th_star)
  th_star <- th_star[1]
  cluster_total <- apply(cluster_num, 2, sum)
  th_hat <- t(t(cluster_num) / cluster_total)[1, ]
  ci_increa <- z*sqrt(th_hat*(1-th_hat)/cluster_total)
  weight_ci <- rbind(th_hat - ci_increa, th_hat + ci_increa)
  Cind <- which(weight_ci[1, ] >= th_star & weight_ci[2, ] <= th_star)
  Aind <- which(weight_ci[2, ] > th_star)
  Bind <- which(weight_ci[1, ] < th_star)
  cluster_type <- rep(NA, length(cluster))
  cluster_type[which(cluster %in% Cind)] <- "C"
  cluster_type[which(cluster %in% Bind)] <- "B"
  cluster_type[which(cluster %in% Aind)] <- "A"
  cluster_type <- factor(cluster_type, levels = c("C", "A", "B"))
  return(cluster_type)
}

# function to calculate clustering acc
acc_func <- function(predicted, truth) {
  confusion_mtx <- table(predicted, truth)
  
  acc <-  sum(diag(confusion_mtx))/sum(confusion_mtx)
  
  # 'TP', 'FP' and 'FN' - 'True Positives', 'False Positives' and 'False Negatives'
  TP <- diag(confusion_mtx); FP <- rowSums(confusion_mtx)-TP; FN <- colSums(confusion_mtx)-TP
  # 'TN' - 'True Negatives'
  TN <- sapply(1:length(TP), function(y, TP){ sum(TP[-y], na.rm = TRUE) }, TP)
  # 'specificity' = TN/RN = TN/(TN+FP)
  recall <- TP/(TP+FN)
  # 'precision' ('PPV') = TP/PP = TP/(TP+FP)
  specificity <- TN/(TN+FP)
  balanced_acc <- (recall+specificity)/2
  
  return(list(confusion_mtx = confusion_mtx, 
              acc = acc, balanced_acc = balanced_acc))
}


getmarker_func <- function(sim_data, cluster, num_marker = 100) {
  sim_data_DE <- SetIdent(sim_data, value = cluster)
  markers = FindAllMarkers(sim_data_DE, only.pos = T)
  
  marker.sortedByPvaladj = markers[order(markers$p_val_adj),]
  filtered_markers_merged_screening <- marker.sortedByPvaladj %>%
    group_by(cluster) %>%
    filter(p_val_adj < 0.05) %>% 
    filter(avg_log2FC > 0.1)
  top_marker_temp <- filtered_markers_merged_screening %>%
    group_by(cluster) %>%
    slice_max(n = round(num_marker/length(unique(cluster))), order_by = avg_log2FC)
  top_marker <- top_marker_temp
  return(top_marker)
}

get_TPR_FDR <- function(truth, predict) {
  truth <- na.omit(truth)
  predict <- na.omit(predict)
  TPR <- length(which(table(c(truth, predict)) > 1)) / length(truth)
  FDR <- 1 - length(which(table(c(truth, predict)) > 1)) / length(predict)
  TPR_FDR <- c(TPR, FDR)
  return(TPR_FDR)
}














