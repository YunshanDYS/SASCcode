# Function to calculate y bar
func_mean <- function(y, z, h) {
  if (length(which(z == h)) > 1) {
    mean <- apply(y[which(z == h), ], 2, mean)
  } else if (length(which(z == h)) == 1) {
    mean <- y[which(z == h), ]
  }
  return(mean)
}

# Function for the posterior inference of the 888 model
func888 <- function(features, y0, y1, H, seed, cell_names, LOCRC_cells, YOCRC_cells, condition_names = c("L", "E"),
                    beta0_low = 1, beta0_high = 10, niter = 2000, niter_extra = 500, iden_eps = TRUE,
                    tmodel = FALSE, t0 = NULL, t1 = NULL, annot = NULL, c2 = 0, c3 = 0, alpha_t = 0.5, alpha_ci = 0.05) {
  set.seed(seed)
  ################
  # initialization using k means
  if (ncol(y0) >= 25) {
    kmeans_out <- kmeans(rbind(y0, y1)[,1:25], centers = 3*H, nstart = 3*H, iter.max = 100)
    xi <- c()
    Sig <- list()
    for (h in 1:(3*H)) {
      hh <- h 
      Sig[[h]] <- diag(apply(features[which(kmeans_out$cluster == hh),], 2, var) + 1e-5)
      xi <- cbind(xi, apply(features[which(kmeans_out$cluster == hh),], 2, mean))
    }
  } else {
    kmeans_out <- kmeans(rbind(y0, y1), centers = 3*H, nstart = 3*H, iter.max = 100)
    xi <- t(kmeans_out$centers) 
    Sig <- list()
    for (h in 1:(3*H)) {
      hh <- h %% H
      if (hh == 0){
        hh <- H
      }
      Sig[[h]] <- diag(apply(gene[which(kmeans_out$cluster == hh),], 2, var)) * 10
    }
  }
  ######################
  n_feat <- dim(y0)[2]
  # hyper parameters
  alpha <- rep(1, H)
  a0 <- H*beta0_high 
  b0 <- H*beta0_high + H*beta0_low 
  mu0 <- apply(features, 2, mean)
  Sig0 <- diag(apply(features, 2, var)) 
  Sig0_inv <- diag(1/(apply(features, 2, var))) 
  nu0 <- n_feat + 2 
  Phi0 <- Sig0 

  ####################
  ## MCMC 
  num0 <- dim(y0)[1]
  num1 <- dim(y1)[1]
  ## initialization
  u00 <- u01 <- u10 <- u11 <- rdirichlet(1, alpha)
  w <- rdirichlet(1, alpha)
  eta0 <- rbeta(1, beta0_high, beta0_low)
  eta1 <- rbeta(1, beta0_low, beta0_high)
  eps0 <- rbeta(1, a0, b0) #1/3
  eps1 <- rbeta(1, a0, b0) #1/3

  prob0_temp <- c(eps0*w, (1 - eps0)*eta0*u00, (1 - eps0)*(1 - eta0)*u01)
  prob1_temp <- c(eps1*w, (1 - eps1)*eta0*u10, (1 - eps1)*(1 - eta1)*u11)
  z0 <- sample(1:(3*H), size = num0, replace = T, prob = prob0_temp)
  z1 <- sample(1:(3*H), size = num1, replace = T, prob = prob1_temp)

  y_mean0 <- apply(y0, 2, mean)
  y_mean1 <- apply(y1, 2, mean)
  
  if (tmodel) {
    J <- length(table(annot))
    pkj <- rdirichlet((3*H), rep(alpha_t, J))
  } else {
    J = 2
    pkj <- matrix(1, nrow = (3*H), ncol = J)
    t0 <- rep(1, length(LOCRC_cells))
    t1 <- rep(1, length(YOCRC_cells))
  }
  
  ## prior params
  Sig_tilde_h <- Sig[[1]]
  xi_tilde_h <- xi[, 1]
  Sig0_inv_mu0 <- Sig0_inv %*% mu0

  pi0 <- matrix(0, nrow = num0, ncol = 3*H)
  pi1 <- matrix(0, nrow = num1, ncol = 3*H)
  m0h <- m1h <- rep(0, 3*H)
  y_bar <- list()

  ## store mcmc samples
  u0_mat <- u1_mat <- matrix(0, nrow = niter, ncol = 2*H)
  w_mat <- matrix(0, nrow = niter, ncol = H)
  z0_mat <- matrix(0, nrow = niter, ncol = num0)
  z1_mat <- matrix(0, nrow = niter, ncol = num1)
  xi_list <- list()
  Sig_list <- list()
  eps0_mat <- eps1_mat <- rep(NA, niter)
  eta0_mat <- eta1_mat <- rep(NA, niter)
  cat("MCMC iterations: \n")
  ####
  for (iter in 1:niter) {
    # update z0, z1
    for (h in 1:H) {
      pi0[ ,h] <- eps0 * w[h] * dmvnorm(y0, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t0]
      pi1[ ,h] <- eps1 * w[h] * dmvnorm(y1, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t1]
    }
    for (h in (H + 1): (2*H)) {
      pi0[ ,h] <- (1 - eps0) * eta0 * u00[h - H] * dmvnorm(y0, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t0]
      pi1[ ,h] <- (1 - eps1) * eta1 * u10[h - H] * dmvnorm(y1, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t1]
    }
    for (h in (2*H + 1): (3*H)) {
      pi0[ ,h] <- (1 - eps0) * (1 - eta0) * u01[h - 2*H] * dmvnorm(y0, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t0]
      pi1[ ,h] <- (1 - eps1) * (1 - eta1) * u11[h - 2*H] * dmvnorm(y1, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t1]
    }
    for (i in 1:num0) {
      z0[i] <- sample(1:(3*H), size = 1, prob = pi0[i, ])
    }
    for (i in 1:num1) {
      z1[i] <- sample(1:(3*H), size = 1, prob = pi1[i, ])
    }
    z0_mat[iter, ] <- z0
    z1_mat[iter, ] <- z1
    
    # update xi
    # mh
    for (h in 1:(3*H)) {
      m0h[h] <- sum(z0 == h)
      m1h[h] <- sum(z1 == h)
    }
    mh <- m0h + m1h
    # y bar
    for (h in 1:(3*H)) {
      if (min(m0h[h], m1h[h]) > 0) {
        y_bar[[h]] <- func_mean(y0, z0, h) + func_mean(y1, z1, h)
      } else if (m0h[h] > 0 & m1h[h] == 0) {
        y_bar[[h]] <- func_mean(y0, z0, h)
      } else if (m0h[h] == 0 & m1h[h] > 0) {
        y_bar[[h]] <- func_mean(y1, z1, h)
      } else {
        y_bar[[h]] <- rep(0, n_feat)
      }
    }
    # xi
    for (h in 1:(3*H)) {
      Sig_h_inv <- solve(Sig[[h]])
      Sig_tilde_h <- solve(Sig0_inv + mh[h] * Sig_h_inv)
      xi_tilde_h <- Sig_tilde_h %*% (Sig0_inv_mu0 + mh[h] * Sig_h_inv %*% y_bar[[h]])
      xi[, h] <- rmvnorm(1, mean = xi_tilde_h, sigma = Sig_tilde_h)
    }
    xi_list[[iter]] <- xi
    
    # update Sig
    for (h in 1:(3*H)) {
      nu_tilde_h <- nu0 + mh[h]
      Vh_temp <- matrix(0, nrow = n_feat, ncol = n_feat)
      if (mh[h] != 0) {
        if (m0h[h] != 0) {
          for (i in which(z0 == h)) {
            diff_temp <- matrix(y0[i, ] - xi[, h], ncol = 1)
            Vh_temp <- Vh_temp + diff_temp %*% t(diff_temp)
          }
        }
        if (m1h[h] != 0) {
          for (i in which(z1 == h)) {
            diff_temp <- matrix(y1[i, ] - xi[, h], ncol = 1)
            Vh_temp <- Vh_temp + diff_temp %*% t(diff_temp)
          }
        }
        diff_temp <- matrix(xi[, h] - y_bar[[h]], ncol = 1)
        Phi_tilde_h <- Phi0 + Vh_temp + mh[h] * diff_temp %*% t(diff_temp)
      } else {
        Phi_tilde_h <- Phi0 + Vh_temp
      }
      
      # Sig
      Sig[[h]] <- riwish(nu_tilde_h, Phi_tilde_h)
    }
    Sig_list[[iter]] <- Sig
    
    # update w, u0, u1, eta0, eta1
    w <- rdirichlet(1, alpha + mh[1:H])
    w_mat[iter, ] <- w
    
    eta0 <- rbeta(1, beta0_high + sum(m0h[(H+1):(2*H)]), beta0_low + sum(m0h[(2*H+1):(3*H)]))
    eta1 <- rbeta(1, beta0_low + sum(m1h[(H+1):(2*H)]), beta0_high + sum(m1h[(2*H+1):(3*H)]))
    eta0_mat[iter] <- eta0
    eta1_mat[iter] <- eta1
    
    u00 <- rdirichlet(1, alpha + m0h[(H+1):(2*H)])
    u01 <- rdirichlet(1, alpha + m0h[(2*H+1):(3*H)])
    u10 <- rdirichlet(1, alpha + m1h[(H+1):(2*H)])
    u11 <- rdirichlet(1, alpha + m1h[(2*H+1):(3*H)])
    u0_mat[iter, ] <- c(u00, u01)
    u1_mat[iter, ] <- c(u10, u11)
    
    # update eps0, eps1
    if (iden_eps) {
      eps0 <- eps1 <- rbeta(1, a0 + sum(mh[1:H]), b0 + sum(mh[(H+1):(3*H)]))
    } else {
      eps0 <- rbeta(1, a0 + sum(m0h[1:H]), b0 + sum(m0h[(H+1):(3*H)]))
      eps1 <- rbeta(1, a0 + sum(m1h[1:H]), b0 + sum(m1h[(H+1):(3*H)]))
    }
    eps0_mat[iter] <- eps0
    eps1_mat[iter] <- eps1
    
    # update pkj
    if (tmodel) {
      for (k in 1:(3*H)) {
        nt0k <- nt1k <- rep(NA, J)
        for (j in 1:J) {
          nt0k[j] <- sum((z0 == k & t0 == j))
          nt1k[j] <- sum((z1 == k & t1 == j))
        }
        ntk <- nt0k + nt1k
        pkj[k, ] <- rdirichlet(1, rep(alpha_t, J) + ntk)
      }
    }
    
    if (iter %% 50 == 0) {
      cat(iter, "\t")
    }
  }

  iters_left <- seq(niter/2, niter, by = 2)
  z0_mat_post <- z0_mat[iters_left, ]
  z1_mat_post <- z1_mat[iters_left, ]

  # # take a look at the convergence of the mc chain
  # plot(w_mat[,1], type = "l")
  # plot(w_mat[,2], col = 'brown', type = "l")

  ##################################
  ## point estimate of the partition
  cat("\n Calculating the loss of the partitions... \n")
  z_mat_post <- cbind(z0_mat_post, z1_mat_post)
  ncells_total <- ncol(z_mat_post)
  nmc_left <- nrow(z_mat_post)
  loss <- rep(NA,nmc_left)
  for (j in 1:nmc_left) {     
    loss[j] <- partition.loss(truth = z_mat_post, estimate = z_mat_post[j, ], loss = VI())
  }
  # calculate added loss
  if (c2 + c3 == 0) {
    z_post <- z_mat_post[which.min(loss)[1], ]
    z0_post <- z_post[1:ncol(z0_mat_post)]
    z1_post <- z_post[(ncol(z0_mat_post)+1):ncells_total]
  } else {
    names(annot) <- cell_names
    annot_reorder <- c(annot[LOCRC_cells], annot[YOCRC_cells])
    # top 50 partitions based on VI loss in salso
    top_index <- order(loss, decreasing = F)[1:50]
    top_loss <- loss[top_index]
    added_loss <- top_loss
    ncells_avg <- ncells_total/H/3
    for (ind in 1:100) {
      par_temp <- z_mat_post[top_index[ind], ]
      nk <- rep(0, (3*H))
      for (kk in 1:(3*H)) {
        nk[kk] <- length(which(par_temp == kk))
      }
      nk[which(nk > ncells_avg/4)] <- ncells_avg
      loss_ncell <- sum(ncells_avg - nk)/ncells_avg
      loss_homo <- 0
      for (kk in 1:(3*H)) {
        if (length(which(par_temp == kk)) > 1) {
          lk <- table(annot_reorder[which(par_temp == kk)])
          lk <- c(lk/sum(lk))
          loss_homo <- loss_homo + (1 - max(lk))
        }
      }
      loss_homo <- loss_homo
      loss_temp <- c2*loss_ncell + c3*loss_homo
      added_loss[ind] <- added_loss[ind] + loss_temp
    }
    z_post <- z_mat_post[top_index[which.min(added_loss)[1]], ]
    z0_post <- z_post[1:ncol(z0_mat_post)]
    z1_post <- z_post[(ncol(z0_mat_post)+1):ncells_total]
  }
  
  #####################################
  ## calculate the cluster specific parameters fixing the partition as the point estimate
  z0 <- z0_post
  z1 <- z1_post
  # calculate estimated xi
  # mh
  for (h in 1:(3*H)) {
    m0h[h] <- sum(z0 == h)
    m1h[h] <- sum(z1 == h)
  }
  mh <- m0h + m1h
  # y bar
  for (h in 1:(3*H)) {
    if (min(m0h[h], m1h[h]) > 0){
      y_bar[[h]] <- func_mean(y0, z0, h) + func_mean(y1, z1, h)
    } else if (m0h[h] > 0 & m1h[h] == 0){
      y_bar[[h]] <- func_mean(y0, z0, h)
    } else if (m0h[h] == 0 & m1h[h] > 0){
      y_bar[[h]] <- func_mean(y1, z1, h)
    } else {
      y_bar[[h]] <- rep(0, n_feat)
    }
  }
  # xi
  for (h in 1:(3*H)) {
    Sig_h_inv <- solve(Sig[[h]])
    Sig_tilde_h <- solve(Sig0_inv + mh[h] * Sig_h_inv)
    xi_tilde_h <- Sig_tilde_h %*% (Sig0_inv_mu0 + mh[h] * Sig_h_inv %*% y_bar[[h]])
    xi[, h] <- rmvnorm(1, mean = xi_tilde_h, sigma = Sig_tilde_h)
  }
  xi_est <- xi
  # calculate estimated Sig
  for (h in 1:(3*H)) {
    nu_tilde_h <- nu0 + mh[h]
    Vh_temp <- matrix(0, nrow = n_feat, ncol = n_feat)
    if (mh[h] != 0){
      if (m0h[h] != 0) {
        for (i in which(z0 == h)) {
          diff_temp <- matrix(y0[i, ] - xi[, h], ncol = 1)
          Vh_temp <- Vh_temp + diff_temp %*% t(diff_temp)
        }
      }
      if (m1h[h] != 0) {
        for (i in which(z1 == h)) {
          diff_temp <- matrix(y1[i, ] - xi[, h], ncol = 1)
          Vh_temp <- Vh_temp + diff_temp %*% t(diff_temp)
        }
      }
      diff_temp <- matrix(xi[, h] - y_bar[[h]], ncol = 1)
      Phi_tilde_h <- Phi0 + Vh_temp + mh[h] * diff_temp %*% t(diff_temp)
    } else {
      Phi_tilde_h <- Phi0 + Vh_temp
    }
    
    
    # Sig
    Sig[[h]] <- riwish(nu_tilde_h, Phi_tilde_h)
  }
  Sig_est <- Sig
  cat("\n Got the point estimate of the partition! \n")

  ####################################################
  ## Uncertainty quantification

  ## more mcmc samples to calculate the posterior mean weights
  cat('\n Extra mcmc for estimation of the weights. \n')
  ## store mcmc samples
  u0_mat <- u1_mat <- matrix(0, nrow = niter_extra, ncol = 2*H)
  w_mat <- matrix(0, nrow = niter_extra, ncol = H)
  z0_mat <- matrix(0, nrow = niter_extra, ncol = num0)
  z1_mat <- matrix(0, nrow = niter_extra, ncol = num1)
  eps0_mat <- eps1_mat <- eta0_mat <- eta1_mat <- rep(NA, niter_extra)

  pi0_htemp_list <- pi1_htemp_list <- list()
  for (h in 1:(3*H)) {
    pi0_htemp_list[[h]] <- matrix(0, nrow = niter_extra, ncol = num0)
    pi1_htemp_list[[h]] <- matrix(0, nrow = niter_extra, ncol = num1)
  }

  for (iter in 1:niter_extra) {
    # update z0, z1
    for (h in 1:H) {
      pi0[ ,h] <- eps0 * w[h] * dmvnorm(y0, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t0]
      pi1[ ,h] <- eps1 * w[h] * dmvnorm(y1, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t1]
    }
    for (h in (H + 1): (2*H)) {
      pi0[ ,h] <- (1 - eps0) * eta0 * u00[h - H] * dmvnorm(y0, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t0]
      pi1[ ,h] <- (1 - eps1) * eta1 * u10[h - H] * dmvnorm(y1, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t1]
    }
    for (h in (2*H + 1): (3*H)) {
      pi0[ ,h] <- (1 - eps0) * (1 - eta0) * u01[h - 2*H] * dmvnorm(y0, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t0]
      pi1[ ,h] <- (1 - eps1) * (1 - eta1) * u11[h - 2*H] * dmvnorm(y1, mean = xi[, h], sigma = Sig[[h]]) * pkj[h, t1]
    }
    for (i in 1:num0) {
      z0[i] <- sample(1:(3*H), size = 1, prob = pi0[i, ])
    }
    for (i in 1:num1) {
      z1[i] <- sample(1:(3*H), size = 1, prob = pi1[i, ])
    }
    z0_mat[iter, ] <- z0
    z1_mat[iter, ] <- z1
    
    for (h_temp in 1:(3*H)) {
      if (h_temp <= H) {
        pi0_htemp_list[[h_temp]][iter, ] <- log(eps0 * w[h_temp]) + dmvnorm(y0, mean = xi[, h_temp], sigma = Sig[[h_temp]], log = T) + log(pkj[h_temp, t0])
        pi1_htemp_list[[h_temp]][iter, ] <- log(eps1 * w[h_temp]) + dmvnorm(y1, mean = xi[, h_temp], sigma = Sig[[h_temp]], log = T) + log(pkj[h_temp, t1])
      } else if (h_temp %in% (H+1):(2*H)) {
        pi0_htemp_list[[h_temp]][iter, ] <- log((1 - eps0) * eta0 * u00[h_temp - H]) + dmvnorm(y0, mean = xi[, h_temp], sigma = Sig[[h_temp]], log = T) + log(pkj[h_temp, t0])
        pi1_htemp_list[[h_temp]][iter, ] <- log((1 - eps1) * eta1 * u10[h_temp - H]) + dmvnorm(y1, mean = xi[, h_temp], sigma = Sig[[h_temp]], log = T) + log(pkj[h_temp, t1])
      } else {
        pi0_htemp_list[[h_temp]][iter, ] <- log((1 - eps0) * (1 - eta0) * u01[h_temp - 2*H]) + dmvnorm(y0, mean = xi[, h_temp], sigma = Sig[[h_temp]], log = T) + log(pkj[h_temp, t0])
        pi1_htemp_list[[h_temp]][iter, ] <- log((1 - eps1) * (1 - eta1) * u11[h_temp - 2*H]) + dmvnorm(y1, mean = xi[, h_temp], sigma = Sig[[h_temp]], log = T) + log(pkj[h_temp, t1])
      }
    }
    
    # mh
    for (h in 1:(3*H)) {
      m0h[h] <- sum(z0 == h)
      m1h[h] <- sum(z1 == h)
    }
    mh <- m0h + m1h
    # update w, u0, u1
    w <- rdirichlet(1, alpha + mh[1:H])
    w_mat[iter, ] <- w
    
    eta0 <- rbeta(1, beta0_high + sum(m0h[(H+1):(2*H)]), beta0_low + sum(m0h[(2*H+1):(3*H)]))
    eta1 <- rbeta(1, beta0_low + sum(m1h[(H+1):(2*H)]), beta0_high + sum(m1h[(2*H+1):(3*H)]))
    eta0_mat[iter] <- eta0
    eta1_mat[iter] <- eta1
    
    u00 <- rdirichlet(1, alpha + m0h[(H+1):(2*H)])
    u01 <- rdirichlet(1, alpha + m0h[(2*H+1):(3*H)])
    u10 <- rdirichlet(1, alpha + m1h[(H+1):(2*H)])
    u11 <- rdirichlet(1, alpha + m1h[(2*H+1):(3*H)])
    u0_mat[iter, ] <- c(u00, u01)
    u1_mat[iter, ] <- c(u10, u11)
    
    # update eps0, eps1
    if (iden_eps) {
      eps0 <- eps1 <- rbeta(1, a0 + sum(mh[1:H]), b0 + sum(mh[(H+1):(3*H)]))
    } else {
      eps0 <- rbeta(1, a0 + sum(m0h[1:H]), b0 + sum(m0h[(H+1):(3*H)]))
      eps1 <- rbeta(1, a0 + sum(m1h[1:H]), b0 + sum(m1h[(H+1):(3*H)]))
    }
    eps0_mat[iter] <- eps0
    eps1_mat[iter] <- eps1
    
    # update pkj
    if (tmodel) {
      for (k in 1:(3*H)) {
        nt0k <- nt1k <- rep(NA, J)
        for (j in 1:J) {
          nt0k[j] <- sum((z0 == k & t0 == j))
          nt1k[j] <- sum((z1 == k & t1 == j))
        }
        ntk <- nt0k + nt1k
        pkj[k, ] <- rdirichlet(1, rep(alpha_t, J) + ntk)
      }
    }
    
    if (iter %% 50 == 0) {
      cat(iter, "\t")
    }
  }
  
  weights_mat <- list()
  weights_mat[[1]] <- cbind(eps0_mat*w_mat, (1-eps0_mat)*eta0_mat*u0_mat[, 1:H], (1-eps0_mat)*(1-eta0_mat)*u0_mat[, (H+1):(2*H)] )
  weights_mat[[2]] <- cbind(eps1_mat*w_mat, (1-eps1_mat)*eta1_mat*u1_mat[, 1:H], (1-eps1_mat)*(1-eta1_mat)*u1_mat[, (H+1):(2*H)] )
  
  weights <- rbind(apply(weights_mat[[1]], 2, mean), apply(weights_mat[[2]], 2, mean) )
  
  # identify cluster type
  weights_ci <- rbind(apply(weights_mat[[1]] - weights_mat[[2]], 2, quantile, alpha_ci/2),
                      apply(weights_mat[[1]] - weights_mat[[2]], 2, quantile, 1 - alpha_ci/2))
  cluster_type_CLE <- rep(NA, ncol(weights_ci))
  for (i in 1:ncol(weights_ci)) {
    if (weights_ci[1,i] <= 0 & weights_ci[2,i] >= 0) {
      cluster_type_CLE[i] <- "C"
    } else if (weights_ci[2,i] < 0) {
      cluster_type_CLE[i] <- condition_names[2]
    } else {
      cluster_type_CLE[i] <- condition_names[1]
    }
  }

  return(out = list(z0_post = z0_post, z1_post = z1_post, weights = weights,
  pi0_htemp_list = pi0_htemp_list, pi1_htemp_list = pi1_htemp_list,
  xi_est = xi_est, Sig_est = Sig_est,
  z0_mat_fixxi = z0_mat, z1_mat_fixxi = z1_mat,
  weights_mat = weights_mat, cluster_type_CLE = cluster_type_CLE))
}




