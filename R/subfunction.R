phy.H.rel <- function(dat, tmp, q, rtreephy, wk, formula){
  n <- sum(dat)
  N <- ncol(dat)
  ga <- tmp$Gamma[,1]; gL <- tmp$Gamma[,2]
  names(ga) <- rownames(tmp$Gamma)
  names(gL) <- rownames(tmp$Gamma)
  gp = ga/n; TT = sum(gp*gL);
  aa <- sapply(tmp$Alpha, function(X){
    abun <- X[,1]
    names(abun) = rownames(X)
    return(abun)
  })
  B <- length(ga)
  pik <- sapply(1:N, function(k) {
    abun <- aa[,k]
    nk <- sum(dat[,k])
    return(abun/nk)
  }) #zik|z+k
  weight.p <- pik*t(replicate(B,wk)) #(zik|z+k)*wk => zik|z++
  if(formula == "mle"){
    qDk <- apply(pik, 2, get("mle.phy.q"), LL = gL, TT = TT, q) ##\sum{p^q} or -\sum{p*log(p)})
    qDg <- mle.phy.q(gp, LL = gL, TT = TT, q)
  } else if(formula == "est"){
    qDk <- sapply(seq_len(ncol(dat)), function(k){
      nk <- sum(dat[,k])
      get("est.phy.q")(xx = aa[,k], LL = gL, TT = TT, n = nk, q, rtreephy)
    })
    qDg <- est.phy.q(ga, LL = gL, TT = TT, n = n, q, rtreephy)
  }
  if(q!=1){
    alpha.R <- (wk%*%qDk)^(1/(1-q))
    joint.R <- (wk^q%*%qDk)^(1/(1-q))
    gamma.R <- qDg^(1/(1-q))
  } else{
    alpha.R <- exp(wk%*%(qDk/TT+log(TT)))
    joint.R <- exp(-wk%*%log(wk) + wk%*%(qDk/TT+log(TT)))
    gamma.R <- exp(qDg+log(TT))
  }
  beta.R <- gamma.R/alpha.R
  betamax.R <- joint.R/alpha.R
  
  ifelse(q==1, C_1 <- log(beta.R)/log(betamax.R), C_1 <- (beta.R^(1-q)-1)/(betamax.R^(1-q)-1))
  ifelse(q==1, U_1 <- log(beta.R)/log(betamax.R), U_1 <- (beta.R^(q-1)-1)/(betamax.R^(q-1)-1))
  U_1/C_1
  V_1 <- (beta.R-1)/(betamax.R-1)
  S_1 <- (beta.R^-1-1)/(betamax.R^-1-1)
  diff <- c("1-CqN" = C_1, "1-UqN" = U_1, "1-VqN" = V_1, "1-SqN" = S_1)
  out <- c("Gamma" = gamma.R, "Alpha" = alpha.R, "Beta" = beta.R, diff)
  
  return(out)
}


phy.H.abs <- function(dat, tmp, q, rtreephy, wk, formula){
  n <- sum(dat)
  N <- ncol(dat)
  ga <- tmp$Gamma[,1]; gL <- tmp$Gamma[,2]
  names(ga) <- rownames(tmp$Gamma)
  names(gL) <- rownames(tmp$Gamma)
  gp = ga/n; TT = sum(gp*gL);
  aa <- sapply(tmp$Alpha, function(X){
    abun <- X[,1]
    names(abun) = rownames(X)
    return(abun)
  })
  B <- length(ga)
  pik <- aa/n #zik|z++
  if(formula == "mle"){
    qDa <- mle.phy.q(c(pik), LL = rep(gL, N), TT = TT, q) ##\sum{p^q} or -\sum{p*log(p)})
    qDg <- mle.phy.q(gp, LL = gL, TT = TT, q)
  } else if(formula == "est"){
    aa_joint <- c(aa) %>% `names<-`(rep(rownames(aa), N))
    qDa <- est.phy.q(aa_joint, LL = rep(gL, N), TT = TT, n = n, q, rtreephy)
    qDg <- est.phy.q(ga, LL = gL, TT = TT, n = n, q, rtreephy)
  }
  if(q!=1){
    alpha.C <- qDa^(1/(1-q))/N
    gamma.C <- qDg^(1/(1-q))
  } else{
    alpha.C <- exp(qDa/TT+log(TT))/N
    gamma.C <- exp(qDg+log(TT))
  }
  beta.C <- gamma.C/alpha.C
  
  ifelse(q==1, C_1 <- log(beta.C)/log(N), C_1 <- (beta.C^(1-q)-1)/(N^(1-q)-1))
  ifelse(q==1, U_1 <- log(beta.C)/log(N), U_1 <- (beta.C^(q-1)-1)/(N^(q-1)-1))
  
  V_1 <- (beta.C-1)/(N-1)
  S_1 <- (beta.C^-1-1)/(N^-1-1)
  diff <- c("1-CqN" = C_1, "1-UqN" = U_1, "1-VqN" = V_1, "1-SqN" = S_1)
  out <- c("Gamma" = gamma.C, "Alpha" = alpha.C, "Beta" = beta.C, diff)
  
  return(out)
}


mle.phy.q <- function(pp, LL, TT, q){
  LL <- LL[pp>0]
  pp <- pp[pp>0]
  if(q!=1){
    LL%*%((pp/TT)^q)
    #LL*((pp/T)^q)
  }else{
    -LL%*%(pp*log(pp))
    #-LL*((pp/T)*log(pp/T))
  }
}


est.phy.q <- function(xx, LL, TT, n, q, rtreephy){ # proposed
  LL <- LL[names(xx)]
  tmp <- data.frame(abun = xx, length = LL)
  PD_obs <- sum(tmp[tmp[,1]>0,2])
  #f1 <- ifelse(datatype=='incidence', sum(rowSums(data)==1), sum(data==1) )
  f1 = sum(tmp[, 1]==1)
  f2 = sum(tmp[, 1]==2)
  g1 = sum(tmp[tmp[, 1]==1, 2])
  g2 = sum(tmp[tmp[, 1]==2, 2])
  node <- names(rtreephy$parts)
  node2 = names(xx)[xx==2][names(xx)[xx==2] %in% node]
  if(length(node2)>0){
    rep2 <- sapply(node2, function(tx){
      sum(xx[rtreephy$parts[[tx]]] %in% c(2,0)) == length(rtreephy$parts[[tx]])
    })
    f2 <- f2 - sum(rep2)
  }
  node1 = names(xx)[xx==1][names(xx)[xx==1] %in% node]
  if(length(node1)>0){
    rep1 <- sapply(node1, function(tx){
      sum(xx[rtreephy$parts[[tx]]] %in% c(1,0)) == length(rtreephy$parts[[tx]])
    })
    f1 <- f1 - sum(rep1)
  }
  if(f2 > 0){
    A = 2*f2/((n-1)*f1+2*f2)
  }else if(f2 == 0 & f1 > 0){
    A = 2/((n-1)*(f1-1)+2)
  }else{
    A = 1
  }
  t_bar <- sum(xx*LL/n)
  # position0 <- which(q==0)
  # position1 <- which(q==1)
  # position_else = c(1:length(q))[-c(position1)]
  #position_else = c(1:length(q))[-c(position0,position1)]
  #ans = rep(0, length(q))
  if(q==0){
    ans = PD_obs+ifelse(g2>0, (n-1)/n*g1^2/(2*g2), (n-1)/n*g1*(f1-1)/2*(f2-1))
  }else if(q==1){
    q1 = sum(sapply(1:(n-1), function(r) {(1-A)^r/r} ))
    if(A < 1) h2 = (g1/n)*((1-A)^(-n+1))*(-log(A)-q1)
    if(A == 1) h2 = 0
    tmp2 = subset(tmp, tmp[,1]>=1 & tmp[,1]<=(n-1) )
    tmp2 = cbind(tmp2, sapply(tmp2[,1], function(x) sum( 1/(x:(n-1)))))
    h1 = sum(apply(tmp2, 1, prod))/n
    h = h1+h2
    #ans[position1] = t_bar*exp(h/t_bar)
    ans = h
  } else{
    r = 0 : (n-1)
    de = delta(tmp, r , n , A)
    a = sum( choose(q-1, r)*(-1)^r*de )/((t_bar)^q)
    if(A < 1) b = (g1*((1-A)^(1-n))/n)*(A^(q-1)-sum(choose(q-1, r)*(A-1)^r))/((t_bar)^q)
    if(A == 1) b = 0
    ans = a+b
  }
  return( ans )
}


bootstrap.q.Beta <- function(data, mat, rtree, tmp, q, nboot, wij, formula, method){
  H <- nrow(mat)
  out = array(0, dim=c(7*H, length(q), nboot))
  pool <- rowSums(data)
  rtreephy <- newick2phylog(convToNewick(rtree))
  #if(datatype == "abundance"){
  n = colSums(data) ; N = ncol(data)
  pop = Boots.pop(data, rtree, tmp$Alpha)
  S = nrow(data)
  #B = length(c(phytree$leaves,phytree$parts))
  if(pop$unseen == 0) p = pop$p[1:S,]
  if(pop$unseen != 0) p = pop$p[c(1:S,tail(1:nrow(pop$p), pop$unseen)),]
  boot.data = array(0, dim = dim(p))
  rownames(boot.data) <- rownames(p)
  S = nrow(data)
  B = S+rtree$Nnode
  for(i in 1:nboot){
    L = pop$L
    if(pop$unseen == 0) p = pop$p[1:S,]
    if(pop$unseen != 0) p = pop$p[c(1:S,(B+1):nrow(pop$p)),]
    boot.data = array(0, dim = dim(p))
    for(j in 1:ncol(p)) boot.data[,j] = rmultinom(1,n[j],p[,j])
    rownames(boot.data) <- rownames(p)
    unseen = boot.data[-(1:S),]
    boot.data.obs <- boot.data[1:S,]
    #boot.data.obs <- boot.data.obs[rowSums(boot.data.obs)>0, ]
    #tip.boot <- names(rtreephy$leaves)[!names(rtreephy$leaves)%in%rownames(boot.data.obs)]
    #rtree.boot <- drop.tip(rtree, tip.boot)
    #rtreephy.boot <- newick2phylog(convToNewick(rtree.boot))
    boot.datatmp = apply(boot.data.obs, 2, function(x){
      #names(x) = names(rtreephy$leaves)
      tmp <- choose_data(x, rtreephy)
      abun <- tmp[ ,1]
      names(abun) <- rownames(tmp)
      return(abun)
    })
    boot.datatmp = rbind(boot.datatmp, unseen)
    boot.gamma = rowSums(boot.datatmp)
    #boot.gamma = boot.gamma[boot.gamma]
    #boot.datatmp <- boot.datatmp[names(boot.gamma), ]
    boot.gamma <- boot.gamma[boot.gamma>0]
    boot.datatmp <- boot.datatmp[names(boot.gamma), ]
    L <- L[names(boot.gamma), ]
    L.gamma = rowSums(boot.datatmp[names(boot.gamma), ] * L) / boot.gamma
    boot.alpha <- apply(boot.datatmp, 2, function(s) data.frame(branch_abun = s, branch_length = L.gamma))
    boot.gamma <- data.frame(branch_abun = boot.gamma, branch_length = L.gamma)
    boot.tmp <- list(Gamma = boot.gamma, Alpha = boot.alpha)
    #sapply(q, function(qq) method(dat = boot.data, mat, boot.tmp, qq, rtreephy = rtreephy, wij, formula))
    out[,,i] = sapply(q, function(qq) method(dat = boot.data, mat, boot.tmp, qq, rtreephy = rtreephy, wij, formula))
  }
  #print(sum(is.infinite(apply(out, 3, sum))))
  #out[ , ,!is.infinite(apply(out, 3, sum))]
  #}
  return(out)
}


