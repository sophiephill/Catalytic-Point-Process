
# setwd("~/Desktop/smrf/code/")
# source("Cmp_functions.R")
# source("simetas2022.r")

# set.seed(1234)
# unif_t_length = 2;
# theta_beta = 5 # exp rate
# a3 = 0; b3 = 200; # gauss
# c = 2; p = 2.5 # pareto
# theta_K = 0.5; mu = 1;
# X1 = 1; Y1 = 1
# T=100
# 
# z0 = simhawk(x1=X1,y1=Y1, T=T, gmi=pointprod, gt=expgt, gxy = unifxy)
# n0 <- z0$n


misd <- function(tpp, grange, tot_time, num_iter, bw = 3, normN = T, verbose = T, plott = F) {
  supDist <- function (x, y) return (max (abs (x - y)))
  if (min(grange) > 0) grange <- c(0, grange)
  N <- length(tpp)
  delta_t <- diff(grange)  # bin width of histogram estimator g based on input grange
  g <- double(length(delta_t))  # intialize triggering function
  tdiffsmtx <- outer(tpp, tpp, FUN = "-")  # N by N matrix of times differences for tpp 
  tdiffs <- tdiffsmtx[lower.tri(tdiffsmtx, diag = F)]	
  
  A <- list(0)  
  bin_ind <- rep(length(A) + 1, length(tdiffs))
  converge = 0
  for(i in 2:length(grange)) {
    inds <- which(tdiffs > grange[i-1] & tdiffs <= grange[i]) # Find the indices in each bin
    if ((length(inds) == 1) & is.na(inds[1])) inds <- 0
    A[[i-1]] <- inds
    bin_ind[inds] <- i - 1
  }  
  
  # intitialize p_{ij} = 1/i for j <= i and 0 otherwise
  Pmtx_new <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    Pmtx_new[i,] <- rep(1/i, N)
  }
  P_new <- Pmtx_new[lower.tri(Pmtx_new)] #note: diag = F is the default
  
  # start convergence algorithm 
  for(k in 1:num_iter) {
    Pmtx_old <- Pmtx_new
    P_old <- P_new		
    
    for(j in 1:length(delta_t)) {
      Z <- delta_t[j] * ifelse(normN, N, sum(P_old))
      g[j] <- sum(P_old[A[[j]]]) / Z   
    }
    
    
    mu <- sum(diag(Pmtx_old)) / tot_time
    
    # take the probabilities from the M-step and put them back in their right place in the g matrix using the bin indicator 
    fix <- g[bin_ind]
    fix[is.na(fix)] <- 0 
    
    gmtx <- matrix(0, nrow = N, ncol = N)
    gmtx[lower.tri(gmtx)] <- fix 
    
    # set the diagonal elements of the g matrix to the new background rate from the E-step 
    diag(gmtx) <- mu
    
    # E-step
    Pmtx_new <- matrix(0, nrow = N, ncol = N)  
    for(i in 1:N) {
      Pmtx_new[i,] <- gmtx[i,] / sum(gmtx[i,])
    }
    
    P_new <- Pmtx_new[lower.tri(Pmtx_new)]
    
    if(k > 2 & k %%10 == 0) {
      sup <- supDist(Pmtx_old, Pmtx_new)
      if(verbose) {
        cat (
          "Iteration: ", k,
          "SupDist: ", formatC(sup, digits = 8, width = 12, forma = "f"),
          "\n")
      }
      if( sup < 1e-5 ) {
        converge = 1
        break
      }	
    }
  }
  if (plott) {
    plot(grange[-1], g)
    lines(grange[-1], gauss.filt(g, bins = grange[-1], bw = bw)) 
  }
  return(list(mu = mu, g = g, grange = grange, delta_t = delta_t, converge = converge))	
  
}


misdInit <- function(tpp, grange, tot_time, g0,
                     num_iter = 1e5, bw = 3, normN = T, verbose = T, plott = F) {
  
  supDist <- function (x, y) return (max (abs (x - y)))
  if (min(grange) > 0) grange <- c(0, grange)
  N <- length(tpp)
  delta_t <- diff(grange)  # bin width of histogram estimator g based on input grange
  g <- double(length(delta_t))  # intialize triggering function
  tdiffsmtx <- outer(tpp, tpp, FUN = "-")  # N by N matrix of times differences for tpp 
  tdiffs <- tdiffsmtx[lower.tri(tdiffsmtx, diag = F)]	
  
  tdiffIds <-  matrix(1:N, nrow = N, ncol = N, byrow = TRUE)
  tdiffIds <- tdiffIds[lower.tri(tdiffIds, diag = F)]	
    
  A <- list(0)  
  bin_ind <- rep(length(A) + 1, length(tdiffs))
  Pmtx_new <- matrix(0, nrow = N, ncol = N)
  for (i in 2:length(grange)) {
    inds <- which(tdiffs > grange[i-1] & tdiffs <= grange[i]) # Find the indices in each bin
    if ((length(inds) == 1) & is.na(inds[1])) inds <- 0
    A[[i-1]] <- inds
    bin_ind[inds] <- i - 1
    
    colInds <- sort(tdiffIds[inds])
    cnts <- as.vector(table(colInds))
    Pmtx_new[i-1, unique(colInds)] <- g0[i-1] * cnts
  }  
  P_new <- Pmtx_new[lower.tri(Pmtx_new)] #note: diag = F is the default
  
  # start convergence algorithm 
  converge = 0
  for(k in 1:num_iter) {
    Pmtx_old <- Pmtx_new
    P_old <- P_new		
    
    for(j in 1:length(delta_t)) {
      Z <- delta_t[j] * ifelse(normN, N, sum(P_old))
      g[j] <- sum(P_old[A[[j]]]) / Z   
    }
  
    mu <- sum(diag(Pmtx_old)) / tot_time
    
    # take the probabilities from the M-step and put them back in their right place in the g matrix using the bin indicator 
    fix <- g[bin_ind]
    fix[is.na(fix)] <- 0 
    
    gmtx <- matrix(0, nrow = N, ncol = N)
    gmtx[lower.tri(gmtx)] <- fix 
    
    # set the diagonal elements of the g matrix to the new background rate from the E-step 
    diag(gmtx) <- mu
    
    # E-step
    Pmtx_new <- matrix(0, nrow = N, ncol = N)  
    for(i in 1:N) {
      Pmtx_new[i,] <- gmtx[i,] / sum(gmtx[i,])
    }
    
    P_new <- Pmtx_new[lower.tri(Pmtx_new)]
    
    if(k > 2 & k %%10 == 0) {
      sup <- supDist(Pmtx_old, Pmtx_new)
      if(verbose) {
        cat (
          "Iteration: ", k,
          "SupDist: ", formatC(sup, digits = 8, width = 12, forma = "f"),
          "\n")
      }
      if( sup < 1e-5 ) {
        converge = 1
        break
      }	
    }
  }
  if (plott) {
    plot(grange[-1], g)
    lines(grange[-1], gauss.filt(g, bins = grange[-1], bw = bw)) 
  }
  return(list(mu = mu, g = g, grange = grange, delta_t = delta_t, converge = converge))	
  
}

# grange <- seq(0, 100, length = n0)
# start <- Sys.time()
# out <- misd(z0$t, grange = grange, num_iter = 1e3, tot_time = T)
# print(Sys.time() - start)
# 
# library(gplm)
# m3 = 100 ## I used 1100 for normal and exponential with a normal kernel, 880 for uniform, 500par.
# u = c(1:m3)/m3
# v = kernel.function(u, kernel = "gaussian") # returns apply(dnorm(u),1,prod)
# gf <- myfilter(out$g, v)
# 
# plot(grange, dexp(grange, rate = 1/theta_beta), col = "red", lty = 2,
#      type = "l", xlim = c(0, tstop))
# lines(grange[-1], gauss.filt(out$g, bins = grange[-1], bw = 2), type = "l")
# lines(grange[-1], gf, type = "l")


