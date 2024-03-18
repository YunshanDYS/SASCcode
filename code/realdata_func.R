preprocess <- function (data_seurat, resolution = 1, num_markers_per = 2, 
                        fdr_thhd = 1e-70) {
  cat("Start pre-processing ... \n")
  input <- data_seurat
  ## find DE as biomarkers to include in the features
  # get seurat clusters with lower resolution for further analysis of DE
  input <- NormalizeData(object = input, normalization.method = "LogNormalize", scale.factor = 10000)
  input <- FindVariableFeatures(object = input, selection.method = "vst", nfeatures = 2000)
  all.genes = rownames(input)
  input <- ScaleData(object = input, features = all.genes)
  input <- RunPCA(object = input, npcs = dim(input@reductions$pca)[2]) # dim
  input <- FindNeighbors(object = input, dim = 1:dim(input@reductions$pca)[2])
  input <- FindClusters(object = input, resolution = resolution)
  input <- SetIdent(input, value = input@meta.data$seurat_clusters)
  # find markers
  markers <-  FindAllMarkers(input, min.pct = 0.5)
  nclus <- length(unique(input@meta.data$seurat_clusters))
  
  ## calculate family wise p value for each gene and then contorl FDR
  # gene significant or not in one of the clusters
  markers_fw <- markers %>%
    group_by(gene) %>%
    summarize(fwp = min(min(p_val_adj)*nclus, 1) )
  # calculate FDR
  fdr_out <- p.adjust(markers_fw$fwp, "BH")
  cat("Control FDR", fdr_thhd, "adding", length(which(fdr_out < fdr_thhd)) , "biomarkers \n")
  markers_sel <- markers_fw$gene[which(fdr_out < fdr_thhd)]
  
  # get the normalized gene expressions of the biomarkers selected
  normed_gene <- input@assays$RNA@data
  gene_markers <- t(as.matrix(normed_gene[intersect(markers_sel, rownames(normed_gene)), ]))
  gene_PCA <- input@reductions$pca@cell.embeddings
  # rescale the markers to have the same scale of PCA1
  gene_markers <- gene_markers / sqrt(apply(gene_markers^2, 2, sum)) * sqrt(sum((gene_PCA[,1])^2))
  
  features <- cbind(gene_PCA, gene_markers)
  cat("Finish pre-processing! \n")
  return(list(input = input,
              features = features))
}

run888model <- function(data_seurat, H, resolution = 1, num_markers_per = 2, 
                        seed = 1999, fdr_thhd = 1e-70) {
  set.seed(seed)
  preprocess_out <- preprocess(data_seurat, resolution, num_markers_per, fdr_thhd)
  data <- preprocess_out$input
  features <- preprocess_out$features
  
  cell_names <- data@assays$RNA@counts@Dimnames[[2]]
  onset <- data@meta.data$Onset
  LOCRC_cells <- cell_names[which(onset == "LOCRC")]
  y0 <- as.matrix(features[LOCRC_cells, ])
  YOCRC_cells <- cell_names[which(onset == "YOCRC")]
  y1 <- as.matrix(features[YOCRC_cells, ])
  
  # some parameters (might need to tune)
  # params for the priors of the weights
  beta0_low <- 1
  beta0_high <- length(cell_names)/H/3/2
  # control the threshold of merging cluster in post processing
  # clustered with weights less than weight_tol are merged
  weight_tol_merge <- 1/H/15 
  
  out888 <- func888(features, y0, y1, H, seed, cell_names, 
                    LOCRC_cells, YOCRC_cells, beta0_low=beta0_low, beta0_high=beta0_high)
  
  z0_post <- out888$z0_post
  z1_post <- out888$z1_post
  weights <- out888$weights
  pi0_htemp_list <- out888$pi0_htemp_list
  pi1_htemp_list <- out888$pi1_htemp_list
  z0_mat_fixxi <- out888$z0_mat_fixxi
  z1_mat_fixxi <- out888$z1_mat_fixxi
  weights_mat <- out888$weights_mat
  xi_est <- out888$xi_est
  cluster_type_CLE <- out888$cluster_type_CLE
  
  cluster <- rep(NA, nrow(features))
  names(cluster) <- rownames(features)
  cluster[LOCRC_cells] <- as.character(z0_post)
  cluster[YOCRC_cells] <- as.character(z1_post)
  
  Cind <- which(cluster_type_CLE == "C")
  Lind <- which(cluster_type_CLE == "L")
  Eind <- which(cluster_type_CLE == "E")
  
  cluster_type <- rep(NA, length(cluster))
  cluster_type[which(cluster %in% Cind)] <- "C"
  cluster_type[which(cluster %in% Lind)] <- "L"
  cluster_type[which(cluster %in% Eind)] <- "E"
  cluster_type <- factor(cluster_type, levels = c("C", "E", "L"))
  
  out <- out888
  out$cluster_raw <- cluster
  out$xi_est_raw <- xi_est
  out$weights_raw <- weights
  out$weights_mat_raw <- weights_mat
  out$cluster_type <- cluster_type
  
  ## Uncertainty quantification
  prob_in_h_mat <- matrix(NA, ncol = nrow(features), nrow = 3*H)
  
  for (h_temp in 1:(3*H)){
    pi0_htemp <- pi0_htemp_list[[h_temp]]
    pi1_htemp <- pi1_htemp_list[[h_temp]]
    
    prob_in_h <- rep(NA, nrow(features))
    names(prob_in_h) <- rownames(features)
    prob_in_h[LOCRC_cells] <- apply(pi0_htemp, 2, mean)
    prob_in_h[YOCRC_cells] <- apply(pi1_htemp, 2, mean)
    prob_in_h_mat[h_temp, ] <- prob_in_h
  }
  out$prob_in_h_mat_raw <- prob_in_h_mat
  
  # post processing after 888 model
  # combine the clusters that have too few number of cells to other clusters
  dis_xi <- matrix(1e+10, nrow = 3*H, ncol = 3*H)
  for (i in 1:(3*H-1)) {
    for (j in (i+1):(3*H)) {
      dis_xi[i,j] <- dis_xi[j,i] <- sqrt(sum((xi_est[, i] - xi_est[, j])^2))
    }
  }
  
  colnames(weights) <- (1:(3*H))
  for (i in 1:2) {
    colnames(weights_mat[[i]]) <- (1:(3*H))
  }
  prob_in_h_mat_out <- prob_in_h_mat[(1:(3*H)), ]
  rownames(prob_in_h_mat_out) <- (1:(3*H))
  for (h in 1:(3*H)) {
    if (length(which(cluster == h)) < weight_tol_merge * length(cell_names)) {
      h_left <- which(apply(weights, 2, sum) > 0)
      h_replace <- h_left[which.min(dis_xi[h, h_left])]
      cluster[which(cluster == h)] <- h_replace
      weights[, h_replace] <- weights[, h_replace] + weights[, h]
      weights[, h] <- 0
      for (i in 1:2) {
        weights_mat[[i]][, h_replace] <- weights_mat[[i]][, h_replace] + weights_mat[[i]][, h]
        weights_mat[[i]][, h] <- 0
      }
      prob_in_h_mat_out[h_replace, ] <- log(exp(prob_in_h_mat_out[h_replace, ]) + exp(prob_in_h_mat_out[h, ]))
      prob_in_h_mat_out[h, ] <- 0
    }
  }
  cluster_old <- cluster
  cluster_table <- table(cluster_old)
  weights <- weights[, names(cluster_table)]
  names(cluster_type_CLE) <- (1:(3*H))
  cluster_type_CLE <- cluster_type_CLE[names(cluster_table)]
  for (i in 1:length(cluster_table)) {
    cluster[which(cluster_old == as.numeric(names(cluster_table)[i]))] <- i 
  }
  colnames(weights) <-  (1:length(cluster_table))
  for (i in 1:2) {
    weights_mat[[i]] <- weights_mat[[i]][, names(cluster_table)]
    colnames(weights_mat[[i]]) <- (1:length(cluster_table))
  }
  prob_in_h_mat_out <- prob_in_h_mat_out[names(cluster_table), ]
  rownames(prob_in_h_mat_out) <- (1:length(cluster_table))
  
  colnames(xi_est) <- 1:(3*H) 
  xi_est <- xi_est[, names(cluster_table)]
  colnames(xi_est) <- (1:length(cluster_table))
  
  out$cluster <- cluster
  out$xi_est <- xi_est
  out$weights <- weights
  out$weights_mat <- weights_mat
  out$prob_in_h_mat_out <- prob_in_h_mat_out
  out$cluster_type_CLE <- cluster_type_CLE
  
  return(out)
}

