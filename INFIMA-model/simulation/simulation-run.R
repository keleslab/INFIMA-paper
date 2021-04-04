## ----------------------------------------------------------
##           run simulations for the INFIMA model
## ----------------------------------------------------------


library(data.table)
library(gtools)
library(parallel)
rm(list = ls())


for(scheme in c('NI', 'MI', 'HI')){
  outdir <- paste0('INFIMA-model/simulation/', scheme)

  set.seed(2022)

  setwd(outdir)

  #### get the input data and simulate 
  load('INFIMA-model/RData/EM_input_main.RData')
  load('pseudo_counts.RData')
  load('INFIMA-model/RData/EM_results.RData')


  E.flat <- unlist(E.g)
  p <- as.numeric(table(E.flat))/length(E.flat) 

  ## simulate Y
  x <- params$alpha[1]
  # y <- params$alpha[2]/(p[2]+1)
  y <- params$alpha[2]
  # z <- params$alpha[3]/(p[1]+p[3])
  z <- params$alpha[3]
  tran.mat.a <- matrix(c(x, y, z,
                        y, x, y,
                        z, y, x), nrow = 3, ncol = 3, byrow = T)
  x <- params$alpha.0[1]
  y <- params$alpha.0[2]
  z <- params$alpha.0[3]
  tran.mat.a0 <- matrix(c(x, y, z,
                        y, x, y,
                        z, y, x), nrow = 3, ncol = 3, byrow = T)

  # transition matrix for wild strains
  x <- params$beta[1]
  y <- params$beta[2]
  z <- params$beta[3]
  tran.mat.b <- matrix(c(x, y, z,
                        y, x, y,
                        z, y, x), nrow = 3, ncol = 3, byrow = T)
  x <- params$beta.0[1]
  y <- params$beta.0[2]
  z <- params$beta.0[3]
  tran.mat.b0 <- matrix(c(x, y, z,
                        y, x, y,
                        z, y, x), nrow = 3, ncol = 3, byrow = T)

  sim.Y <- function(R, V){
    if(V == 1){
      tm.a <- tran.mat.a
      tm.b <- tran.mat.b
    }
    else { # null generative
      tm.a <- tran.mat.a0
      tm.b <- tran.mat.b0
    }
    Y <- rep(0,8)
    for(i in c(1,2,3,5,6)){
      Y[i] <- sample(c(-1,0,1), 1, prob = tm.a[R[i]+2,])
    }
    for(i in c(4,7,8)){
      Y[i] <- sample(c(-1,0,1), 1, prob = tm.b[R[i]+2,])
    }
    return(Y)
  }

  #### Full EM model implementation ####
  # res <- ModelFitting(D.g, p.g, B.g.t, Y.g, PI)
  ModelFitting <- function(.D.g, .p.g, .B.g.t, .Y.g, .PI){
    G <- length(.D.g)
    
    ####### The counts for the causal generative model
    # alpha <=> a1, beta <=> b1
    # preprocess the count data
    n0 <- vector('list', G)
    n1 <- vector('list', G)
    n2 <- vector('list', G)
    m0 <- vector('list', G)
    m1 <- vector('list', G)
    m2 <- vector('list', G)
    
    # the following are weighted sum counts with Z.g
    N0 <- rep(0,G)
    N1 <- rep(0,G)
    N2 <- rep(0,G)
    M0 <- rep(0,G)
    M1 <- rep(0,G)
    M2 <- rep(0,G)
    
    for(g in 1:G){
      if(.p.g[g] > 0){
        D1.g <- .D.g[[g]][,c(1,2,3,5,6)]
        D2.g <- .D.g[[g]][,c(4,7,8)]
        if(.p.g[g] == 1){
          D1.g <- t(as.matrix(D1.g))
          D2.g <- t(as.matrix(D2.g))
        }
        n0[[g]] <- apply(D1.g, 1, function(x){sum(x==0)})
        n1[[g]] <- apply(D1.g, 1, function(x){sum(x==1)})
        n2[[g]] <- apply(D1.g, 1, function(x){sum(x==2)})
        m0[[g]] <- apply(D2.g, 1, function(x){sum(x==0)})
        m1[[g]] <- apply(D2.g, 1, function(x){sum(x==1)})
        m2[[g]] <- apply(D2.g, 1, function(x){sum(x==2)})
      }
    }
    
    
    
    ####### The counts for the null generative model
    # alpha.0 <=> a0, beta.0 <=> b0
    N00 <- rep(0,G)
    N01 <- rep(0,G)
    N02 <- rep(0,G)
    M00 <- rep(0,G)
    M01 <- rep(0,G)
    M02 <- rep(0,G)
    
    for(g in 1:G){
      if(.p.g[g] > 0){
        tmp <- abs(.B.g.t[[g]] - .Y.g[[g]])
        tmp1 <- tmp[c(1,2,3,5,6)]
        tmp2 <- tmp[c(4,7,8)]
        
        N00[g] <- sum(tmp1 == 0)
        N01[g] <- sum(tmp1 == 1)
        N02[g] <- sum(tmp1 == 2)
        M00[g] <- sum(tmp2 == 0)
        M01[g] <- sum(tmp2 == 1)
        M02[g] <- sum(tmp2 == 2)
      }
    }
    
    
    ##### the reward terms for enforcing concordance
    lambda_0 <- 0.1
    lambda_1 <- 0.01
    lambda_2 <- 0
    
    ### initial values
    alpha.init <- c(0.6,0.3,0.1) # cannot contain any probability equal to zero
    beta.init <- c(0.5,0.4,0.1)
    alpha <- alpha.init
    beta <- beta.init
    alpha.0.init <- c(0.6,0.3,0.1) # cannot contain any probability equal to zero
    beta.0.init <- c(0.5,0.4,0.1)
    alpha.0 <- alpha.0.init
    beta.0 <- beta.0.init
    
    gamma <- 0.4
    
    Z.g <- vector('list', G)
    for(g in 1:G){
      if(.p.g[g] == 1){
        Z.g[[g]] <- 1
      }
      if(.p.g[g] > 1){
        Z.g[[g]] <- rep(1/.p.g[g], .p.g[g])
      }
    }
    
    theta.g <- Z.g
    V.g <- rep(0, G)
    
    alpha.prev <- alpha
    beta.prev <- beta
    for(em.iter in 1:100){
      # E-step
      for(g in 1:G){
        if(.p.g[g] > 0){
          for(k in 1:.p.g[g]){
            Z.g[[g]][k] <- gamma*theta.g[[g]][k]*prod(c(alpha,beta)^c(n0[[g]][k],n1[[g]][k],n2[[g]][k],m0[[g]][k],m1[[g]][k],m2[[g]][k]))
          }
          
          # null probability
          prob.null <- (1-gamma)*prod(c(alpha.0,beta.0)^c(N00[g],N01[g],N02[g],M00[g],M01[g],M02[g]))
          Z.g[[g]] <- Z.g[[g]]/(sum(Z.g[[g]]) + prob.null)
        }
      }
      
      # M-step
      for(g in 1:G){
        if(.p.g[g] > 0){
          V.g[g] <- sum(Z.g[[g]])
          theta.g[[g]] <- .PI[[g]] + Z.g[[g]]*V.g[g]
          theta.g[[g]] <- theta.g[[g]]/sum(theta.g[[g]])
          N0[g] <- sum(Z.g[[g]]*n0[[g]]) + lambda_0*.p.g[g]
          N1[g] <- sum(Z.g[[g]]*n1[[g]]) + lambda_1*.p.g[g]
          N2[g] <- sum(Z.g[[g]]*n2[[g]]) + lambda_2*.p.g[g]
          M0[g] <- sum(Z.g[[g]]*m0[[g]]) + lambda_0*.p.g[g]
          M1[g] <- sum(Z.g[[g]]*m1[[g]]) + lambda_1*.p.g[g]
          M2[g] <- sum(Z.g[[g]]*m2[[g]]) + lambda_2*.p.g[g]
        }
      }
      alpha <- c(sum(N0*V.g), sum(N1*V.g), sum(N2*V.g))
      beta <- c(sum(M0*V.g), sum(M1*V.g), sum(M2*V.g))
      alpha <- alpha/sum(alpha)
      beta <- beta/sum(beta)
      alpha.0 <- c(sum(N00*(1-V.g)), sum(N01*(1-V.g)), sum(N02*(1-V.g)))
      beta.0 <- c(sum(M00*(1-V.g)), sum(M01*(1-V.g)), sum(M02*(1-V.g)))
      alpha.0 <- alpha.0/sum(alpha.0)
      beta.0 <- beta.0/sum(beta.0)
      
      gamma <- mean(V.g)
      
      # print(em.iter)
      # print(sprintf('alpha = %.2f, %.2f, %.2f', alpha[1], alpha[2], alpha[3]))
      # print(sprintf('beta = %.2f, %.2f, %.2f', beta[1], beta[2], beta[3]))
      # print(sprintf('alpha.0 = %.2f, %.2f, %.2f', alpha.0[1], alpha.0[2], alpha.0[3]))
      # print(sprintf('beta.0 = %.2f, %.2f, %.2f', beta.0[1], beta.0[2], beta.0[3]))
      # print(sprintf('gamma = %.2f', gamma))
      
      if(sum((alpha.prev - alpha)^2) + sum((beta.prev - beta)^2) < 1e-6){
        break
      }
      
      alpha.prev <- alpha
      beta.prev <- beta
    }
    
    params <- list(alpha = alpha, beta = beta,
                  alpha.0 = alpha.0, beta.0 = beta.0,
                  gamma = gamma, em.iter = em.iter,
                  V.g = V.g, Z.g = Z.g)
    return(params)
  }


  # params <- list(alpha = alpha, beta = beta,
  #                alpha.0 = alpha.0, beta.0 = beta.0,
  #                gamma = gamma, em.iter = em.iter)
  sim.once <- function(params, PI, E.g, B.g.t, p.g, filename){

    ######## simulation of data ########
    G <- length(p.g)
    Z.g <- vector('list', G)
    theta.g <- vector('list', G)
    D.g <- vector('list', G)
    Y.g <- vector('list', G)
    V.g <- rbinom(G, 1, prob = params$gamma)
    for(g in 1:G){
      if(p.g[g] > 0){
        ### Theta.g follows Dirichlet(PI.g)
        theta.g[[g]] <- as.vector(rdirichlet(1, PI[[g]]))
        ### Z.g follows mixture of Multinomial and zero point mass
        Z.g[[g]] <- V.g[g]*as.vector(rmultinom(1, 1, theta.g[[g]]))
        ### simulate Y.g
        if(V.g[g] == 0){
          Y.g[[g]] <- sim.Y(B.g.t[[g]], V.g[[g]])
        }
        else{
          R.g <- E.g[[g]][which(Z.g[[g]] == 1),]
          Y.g[[g]] <- sim.Y(R.g, V.g[[g]])
        }
        
        # compute the absolute distance
        D.g[[g]] <- matrix(0, nrow(E.g[[g]]), ncol(E.g[[g]]))
        for(i in 1:nrow(E.g[[g]])){
          D.g[[g]][i,] <- abs(E.g[[g]][i,] - Y.g[[g]])
        }
      }
    }
    
    ######### fit the model ########
    params.sim <- ModelFitting(D.g, p.g, B.g.t, Y.g, PI)
    params$Z.g <- Z.g
    params$V.g <- V.g
    params.real <- params
    save(params.sim, params.real, file = filename)
  }


  res <- mclapply(1:1000, function(x){
    print(x)
    sim.once(params, PI, E.g, B.g.t, p.g, filename = paste0('sim-',x,'.RData'))
  }, mc.cores = 20)



}