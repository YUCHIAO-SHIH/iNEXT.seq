choose_data = function(data, tree){
  
  tmp <- data[names(tree$leaves)]  
  for(i in 1:length(tree$parts)){
    tmp[1+length(tmp)] <- sum(tmp[tree$parts[[i]]])
    names(tmp)[length(tmp)] <- names(tree$parts)[i]
  }
  tmp <- data.frame('branch_abun' = tmp,"branch_length" = c(tree$leaves,tree$nodes))
  
  return(tmp)
}


convToNewick <- function(tree){
  tree<-reorder.phylo(tree,"cladewise")
  n<-length(tree$tip)
  string<-vector(); string[1]<-"("; j<-2
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,2]<=n){
      string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
      if(!is.null(tree$edge.length)){
        string[j]<-paste(c(":",round(tree$edge.length[i],10)), collapse="")
        j<-j+1
      }
      v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
      while(length(v)>0&&k==v[length(v)]){
        string[j]<-")"; j<-j+1
        w<-which(tree$edge[,2]==tree$edge[k,1])
        if(!is.null(tree$edge.length)){
          string[j]<-paste(c(":",round(tree$edge.length[w],10)), collapse="")
          j<-j+1
        }
        v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w
      }
      string[j]<-","; j<-j+1
    } else if(tree$edge[i,2]>=n){
      string[j]<-"("; j<-j+1
    }
  }
  if(is.null(tree$edge.length)) string<-c(string[1:(length(string)-1)], ";")
  else string<-c(string[1:(length(string)-2)],";")
  string<-paste(string,collapse="")
  return(string)
}


TranMul <- function(data, tree){
  rtree = newick2phylog(convToNewick(tree))
  data = cbind(data[names(rtree$leaves), ])
  rdata = apply(data, 1, sum)
  AlphaTmp = list()
  GammaTmp = choose_data(rdata, rtree)
  GammaTbar = sum(GammaTmp[, 1]*GammaTmp[, 2]) / sum(rdata)
  for(i in 1:ncol(data)){
    adata = aadata = data[, i]
    names(aadata) = rownames(data)
    names(adata) = tree$tip.label
    if(length(rtree$leaves)<=2){
      abun <- c(adata[names(adata)%in%rtree$tip.label], root = sum(adata))
      tmp = data.frame('branch_abun'=abun, "branch_length" = c(rtree$edge.length,0))
    } else{
      tmp = choose_data(aadata, rtree)
    }
    AlphaTmp[[i]] = tmp
  }
  output = list(Alpha=AlphaTmp, Gamma=GammaTmp)
  return(output)
}


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
  }) #zik|z+k = pi|k
  wk.pik <- pik*t(replicate(B,wk)) #(zik|z+k)*wk => zik|z++ = pik
  if(formula == "mle"){
    qDk <- apply(pik, 2, get("mle.phy.q"), LL = gL, TT = TT, q) ##\sum{p^q} or -\sum{p*log(p)})
    qDg <- mle.phy.q(rowSums(wk.pik), LL = gL, TT = TT, q)
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


bootstrap.q.Beta <- function(data, rtree, tmp, q, nboot, wk, formula, method){
  out = array(0, dim = c(7, length(q), nboot))
  pool <- rowSums(data)
  rtreephy <- newick2phylog(convToNewick(rtree))
  #if(datatype == "abundance"){
  n = colSums(data) ; N = ncol(data)
  pop = Boots.pop(data, rtree, tmp$Alpha) #boots population
  S = nrow(data)
  #B = length(c(phytree$leaves,phytree$parts))
  if (pop$unseen == 0) p = pop$p[1:S,]
  if (pop$unseen != 0) p = pop$p[c(1:S, tail(1:nrow(pop$p), pop$unseen)),]
  boot.data = array(0, dim = dim(p))
  rownames(boot.data) <- rownames(p)
  
  S = nrow(data)
  B = S+rtree$Nnode
  for(i in 1:nboot){
    L = pop$L
    if (pop$unseen == 0) p = pop$p[1:S,]
    if (pop$unseen != 0) p = pop$p[c(1:S, (B+1):nrow(pop$p)),]
    boot.data = array(0, dim = dim(p))
    for (j in 1:ncol(p)) boot.data[,j] = rmultinom(1, n[j], p[,j])
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
    
    out[,,i] = sapply(q, function(qq) method(dat = boot.data, boot.tmp, qq, rtreephy = rtreephy, wk, formula))
  }
  #print(sum(is.infinite(apply(out, 3, sum))))
  #out[ , ,!is.infinite(apply(out, 3, sum))]
  #}
  return(out)
}


Boots.pop <- function(data, rtree, tmp){
  # if(datatype == "abundance"){
  N = ncol(data); n = colSums(data) #z+k
  pool = rowSums(data); OBS = length(pool) #S_pool_obs
  rtreephy <- newick2phylog(convToNewick(rtree))
  OBS_B <- dim(tmp[[1]])[1] #B_pool_obs
  obs <- colSums(data>0) #S_k_obs
  TT <- sum(tmp[[1]][1]/n[1]*tmp[[1]][2])
  
  F1 = sum(pool == 1); F2 = sum(pool == 2)
  F0 = ifelse(F2==0, F1*(F1-1)/2, F1^2/(2*F2))*(sum(n)-1)/sum(n) #pool assemblage f0 estimate
  F0_N <- round(F0)
  
  f1 = sapply(1:N,function(k) sum(data[,k]==1))
  f2 = sapply(1:N,function(k) sum(data[,k]==2))
  g1 = unlist(lapply(tmp, function(tmp_k) sum(tmp_k[tmp_k[,1]==1,2]) ))
  g2 = unlist(lapply(tmp, function(tmp_k) sum(tmp_k[tmp_k[,1]==2,2]) ))
  C <- ifelse(f2 == 0, C <- 1 - f1/n*(n-1)*(f1-1)/((n-1)*(f1-1) + 2), C <- 1 - f1/n*(n-1)*f1/((n-1)*f1 + 2*f2)) #C_k_hat
  f0 <- ifelse(f2 == 0, f0 <- f1*(f1-1)/2, f0 <- f1^2/(2*f2))*(n-1)/n #each assemblage f0 estimate
  f0_N <- round(f0)
  
  r.data = sapply(1:N, function(k) data[,k]/n[k])
  W = sapply(1:N,function(k) (1-C[k])/sum(r.data[,k]*(1-r.data[,k])^n[k])) #lambda_k_hat
  #Chao1-PD
  g0 = sapply(1:N, function(k)
    if ((2*g2[k]*f1[k])>g1[k]*f2[k]) (n[k]-1)/n[k]*g1[k]^2/2/g2[k]
    else (n[k]-1)/n[k]*g1[k]*(f1[k]-1)/2/(f2[k]+1) )
  
  if (F0 > 0){
    boots.pop = rbind(r.data, matrix(0, ncol=N, nrow=F0_N))
  } #obs species abundance + undetected
  else {boots.pop = r.data}
  L = matrix(0, nrow=(OBS_B+F0_N), ncol=N)
  boots.pop2 = matrix(0, nrow=(OBS_B+F0_N), ncol=N) #obs branch abundance + undetected
  for (i in 1:N)
  {
    if (f0_N[i]>0)
    {
      f0_N[i] = ifelse(f0_N[i]+obs[i]>OBS+F0_N, OBS+F0_N-obs[i], f0_N[i])
      boots.pop[,i][1:OBS] <- r.data[ ,i]*(1-W[i]*(1-r.data[ ,i])^n[i]) #species, p_ik_hat, i<S_k_obs
      I = which(boots.pop[,i]==0) #ai+bi
      II = sample(I, f0_N[i])
      u.p <- (1-C[i])/f0_N[i] #p_ik_hat, i>S_k_obs
      boots.pop[II, i] = rep(u.p, f0_N[i])
      da = boots.pop[1:OBS, i] #corrected observed relative species abundance
      names(da) = rownames(data)
      mat = choose_data(da, rtreephy) #corrected observed relative node abundance
      boots.pop2[,i] = c(mat[,1], boots.pop[,i][-(1:OBS)]) #add undetected species
      F00 = sum(II > OBS) #not detect in pool assemblage
      L[1:nrow(mat), i] = mat[,2]
      if (F00>0){
        index = which(boots.pop2[,i] > 0)[which(boots.pop2[,i] > 0) > nrow(mat)]
        #un.sp <- rownames(data)[II[II<OBS]]
        #g0r <- g0[i]- sum(mat[un.sp,2])
        L[index, i] = g0[i]/F00
      }
    }else{
      L[seq_len(OBS_B), i] = tmp[[i]][ ,2]
      boots.pop2[seq_len(OBS_B), i] = tmp[[i]][,1]
    }
  }
  if (F0_N==0){
    rownames(L) <- rownames(mat)
    rownames(boots.pop2) <- rownames(mat)
  } else{
    rownames(L) <- c(rownames(mat), paste0("u", seq_len(F0_N)))
    rownames(boots.pop2) <- c(rownames(mat), paste0("u", seq_len(F0_N)))
  }
  L[L>TT] <- TT
  return(list(p=boots.pop2, L=L, unseen=F0_N))
  # }
}


transconf = function(Bresult, est, conf){
  est.btse = apply(Bresult, 1, sd)
  est.LCL = est - qnorm(1-(1-conf)/2) * est.btse 
  est.UCL = est + qnorm(1-(1-conf)/2) * est.btse
  # if(any(est.LCL<0)) est.LCL[est.LCL<0] <- 0
  # if(any(est.UCL>1)) est.UCL[est.UCL>1] <- 1
  cbind(est = est, btse=est.btse, LCL=est.LCL, UCL = est.UCL)
}

