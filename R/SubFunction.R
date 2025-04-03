choose_data = function(data, tree){
  
  tmp <- data[names(tree$leaves)]  
  for(i in 1:length(tree$parts)){
    tmp[1+length(tmp)] <- sum(tmp[tree$parts[[i]]])
    names(tmp)[length(tmp)] <- names(tree$parts)[i]
  }
  tmp <- data.frame('branch_abun' = tmp,"branch_length" = c(tree$leaves,tree$nodes))
  
  return(tmp)
}

TranMul = function(data, tree){
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


delta = function(data, k, n, A){
  ans = sapply(1:length(k), function(i){
    if(k[i]<n){      
      data1 = data[data[,1]<=(n-k[i]) &  data[,1]>=1,]
      if( class(data1) == "numeric" ) data1 = t(as.matrix(data1))
      sum( data1[,2]*(data1[,1]/n)*exp(lchoose(n-data1[,1], k[i])-lchoose(n-1, k[i])) ) 
    }else{
      g1 = sum(data[data==1,2])
      g1*(1-A)^(k[i]-n+1)/n 
    }
  })
  return( ans )
}



Boots.pop = function(data, rtree, tmp){
  choose_data = function(data, tree){
    
    tmp <- data[names(tree$leaves)]  
    for(i in 1:length(tree$parts)){
      tmp[1+length(tmp)] <- sum(tmp[tree$parts[[i]]])
      names(tmp)[length(tmp)] <- names(tree$parts)[i]
    }
    tmp <- data.frame('branch_abun' = tmp,"branch_length" = c(tree$leaves,tree$nodes))
    
    return(tmp)
  }
  # if(datatype == "abundance"){
  N = ncol(data); n = colSums(data)
  pool=rowSums(data) ; OBS=length(pool)
  rtreephy <- newick2phylog(convToNewick(rtree))
  OBS_B <- dim(tmp[[1]])[1]  #這邊的 tmp 是傳tmp$alpha 近來
  obs <- colSums(data>0)
  TT <- sum(tmp[[1]][,1]/n[1]*tmp[[1]][,2])   #樹根的總長度
  F1=sum(pool==1);F2=sum(pool==2)
  F0=ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2))*(sum(n)-1)/sum(n)  #pool assemblage f0 estimate
  F0_N <- round(F0)
  f1=sapply(1:N,function(k) sum(data[,k]==1))
  f2=sapply(1:N,function(k) sum(data[,k]==2))
  g1=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==1,2]) ))
  g2=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==2,2]) ))
  C <- ifelse(f2 == 0, C <- 1 - f1/n*(n-1)*(f1-1)/((n-1)*(f1-1) + 2), C <- 1 - f1/n*(n-1)*f1/((n-1)*f1 + 2*f2))
  f0 <- ifelse(f2 == 0, f0 <- f1*(f1-1)/2, f0 <- f1^2/(2*f2))*(n-1)/n
  f0_N <- round(f0)  #跟著第一層的跑
  r.data=sapply(1:N,function(k) data[,k]/n[k]) #relative species abundance
  W=sapply(1:N,function(k) (1-C[k])/sum(r.data[,k]*(1-r.data[,k])^n[k])) #N個lambda
  g0 = sapply(1:N, function(k)  
    if((2*g2[k]*f1[k])>g1[k]*f2[k]) (n[k]-1)/n[k]*g1[k]^2/2/g2[k]
    else (n[k]-1)/n[k]*g1[k]*(f1[k]-1)/2/(f2[k]+1) )
  if(F0>0){boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=F0_N))  
  # 如果混合群落有未被觀察到的物種，就新增多少row上去
  }else{boots.pop=r.data}
  L = matrix(0, nrow=(OBS_B+F0_N), ncol=N) # 包含觀察和未觀察物種的L
  boots.pop2 = matrix(0, nrow=(OBS_B+F0_N), ncol=N) 
  #obs branch abundance + undetected豐度矩陣
  
  #針對每個樣本，調整相對豐度，生成自助重抽樣數據。
  # 如果樣本中有未觀察物種，則進行相對豐度的均分。
  # L：進化樹各節點的分支長度
  for(i in 1:N)
  { 
    if(f0_N[i]>0)
    {
      f0_N[i]=ifelse(f0_N[i]+obs[i]>OBS+F0_N, OBS+F0_N-obs[i],f0_N[i])
      boots.pop[,i][1:OBS] <- r.data[ ,i]*(1-W[i]*(1-r.data[ ,i])^n[i])   #dectected species 修正
      I=which(boots.pop[,i]==0) #ai+bi
      II=sample(I,f0_N[i])
      u.p <- (1-C[i])/f0_N[i]
      boots.pop[II,i]=rep(u.p,f0_N[i])
      da = boots.pop[1:OBS,i] 
      #corrected observed relative species abundance
      names(da) = rownames(data)
      mat = choose_data(da, rtreephy)  #mat[,1]:corrected observed branch_abun
      #corrected observed relative node abundance
      boots.pop2[,i]=c(mat[,1], boots.pop[,i][-(1:OBS)]) #add undetected species 
      F00 = sum(II > OBS) #not detect in pool assemblage 
      L[1:nrow(mat),i] = mat[,2] #mat[,2]:corrected observed branch_length
      if(F00>0 ){
        index = which(boots.pop2[,i] > 0)[which(boots.pop2[,i] > 0) > nrow(mat)]
        #un.sp <- rownames(data)[II[II<OBS]]
        #g0r <- g0[i]- sum(mat[un.sp,2])
        L[index, i] = g0[i]/F00  #多層次論文中的公式p.44 (2.60)
      }
    }else{
      
      ##yayun revise##
      da = boots.pop[1:OBS,i]
      names(da) = rownames(data)
      # mat = hiDIP:::choose_data(da, rtreephy)
      mat = choose_data(da, rtreephy)
      ###############
      
      L[seq_len(OBS_B), i] = tmp[[i]][ ,2] 
      boots.pop2[seq_len(OBS_B), i] = tmp[[i]][ ,1]
    }
  }
  if(F0_N==0){
    rownames(L) <- rownames(mat)
    rownames(boots.pop2) <- rownames(mat)
  } else{
    rownames(L) <- c(rownames(mat), paste0("u", seq_len(F0_N)))
    rownames(boots.pop2) <- c(rownames(mat), paste0("u", seq_len(F0_N)))
  }
  L[L>TT] <- TT
  return(list(p=boots.pop2,L=L,unseen=F0_N))
  # }
}

bootstrap.q.Beta = function(data, mat, rtree, tmp, q, nboot, wij, type, method){
  H <- nrow(mat)
  ##chunyu revise## number of value is 10*H-6
  out = array(0, dim=c((10*H-6), length(q), nboot))
  pool <- rowSums(data)
  rtreephy <- newick2phylog(convToNewick(rtree))
  #if(datatype == "abundance"){
  n = colSums(data) ; N = ncol(data)
  pop = Boots.pop(data, rtree, tmp$Alpha)  #回傳 L,p,unseen
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
    ##chunyu revise##
    if(is.vector(unseen)) unseen = matrix(unseen, nrow = 1) %>% `row.names<-`("u1")
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
    #sapply(q, function(qq) method(dat = boot.data, mat, boot.tmp, qq, rtreephy = rtreephy, wij, type))
    out[,,i] = sapply(q, function(qq) method(dat = boot.data, mat, boot.tmp, qq, rtreephy = rtreephy, wij, type))
  }
  #print(sum(is.infinite(apply(out, 3, sum))))
  #out[ , ,!is.infinite(apply(out, 3, sum))]
  #}
  return(out)
}

transconf = function(Bresult, est, conf){
  est.btse = sd(Bresult)
  est.LCL = est - qnorm(1-(1-conf)/2) * est.btse 
  est.UCL = est + qnorm(1-(1-conf)/2) * est.btse
  # if(any(est.LCL<0)) est.LCL[est.LCL<0] <- 0
  # if(any(est.UCL>1)) est.UCL[est.UCL>1] <- 1
  cbind(est = est, btse=est.btse, LCL=est.LCL, UCL = est.UCL)
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



# TranMul <- function(data, tree){
#   rtree = newick2phylog(convToNewick(tree))
#   data = cbind(data[names(rtree$leaves), ])
#   rdata = apply(data, 1, sum)
#   AlphaTmp = list()
#   GammaTmp = choose_data(rdata, rtree)
#   GammaTbar = sum(GammaTmp[, 1]*GammaTmp[, 2]) / sum(rdata)
#   for(i in 1:ncol(data)){
#     adata = aadata = data[, i]
#     names(aadata) = rownames(data)
#     names(adata) = tree$tip.label
#     if(length(rtree$leaves)<=2){
#       abun <- c(adata[names(adata)%in%rtree$tip.label], root = sum(adata))
#       tmp = data.frame('branch_abun'=abun, "branch_length" = c(rtree$edge.length,0))
#     } else{
#       tmp = choose_data(aadata, rtree)
#     }
#     AlphaTmp[[i]] = tmp
#   }
#   output = list(Alpha=AlphaTmp, Gamma=GammaTmp)
#   return(output)
# }

phy.rel <- function(dat, tmp, q, rtreephy, wk, formula){
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
  
  ##chunyu revise##
  if (formula == "spader"){
    nk <- colSums(dat)
    diff <- est.spader(xa = aa, xg = ga, LL = gL, TT = TT, nk = nk, wk = wk, q, rtreephy)
    out <- c("1-CqN" = diff$C_1, "1-UqN" = diff$U_1)
  }
  
  else {
    if (formula == "mle"){
      qDk <- apply(pik, 2, get("mle.phy.q"), LL = gL, TT = TT, q) ##sum{p^q} or -sum{p*log(p)})
      qDg <- mle.phy.q(rowSums(wk.pik), LL = gL, TT = TT, q)
    }
    else if(formula == "est"){
      qDk <- sapply(seq_len(ncol(dat)), function(k){
        nk <- sum(dat[,k])
        get("est.phy.q")(xx = aa[,k], LL = gL, TT = TT, n = nk, q, rtreephy)
      })
      qDg <- est.phy.q(xx = ga, LL = gL, TT = TT, n = n, q, rtreephy)
    }
    
    if (q!=1) {
      alpha.R <- (wk%*%qDk)^(1/(1-q))
      joint.R <- (wk^q%*%qDk)^(1/(1-q))
      gamma.R <- qDg^(1/(1-q))
    } else {
      alpha.R <- exp(wk%*%(qDk/TT+log(TT)))
      joint.R <- exp(-wk%*%log(wk) + wk%*%(qDk/TT+log(TT)))
      gamma.R <- exp(qDg/TT+log(TT))
    }
    beta.R <- gamma.R/alpha.R
    betamax.R <- joint.R/alpha.R
    
    ifelse(q==1, C_1 <- log(beta.R)/log(betamax.R), C_1 <- (beta.R^(1-q)-1)/(betamax.R^(1-q)-1))
    ifelse(q==1, U_1 <- log(beta.R)/log(betamax.R), U_1 <- (beta.R^(q-1)-1)/(betamax.R^(q-1)-1))
    
    V_1 <- (beta.R-1)/(betamax.R-1)
    S_1 <- (beta.R^-1-1)/(betamax.R^-1-1)
    diff <- c("1-CqN*" = C_1, "1-UqN*" = U_1, "1-VqN*" = V_1, "1-SqN*" = S_1)
    out <- c("Gamma" = gamma.R, "Alpha" = alpha.R, "Beta" = beta.R, diff)
  }
  
  return(out)
}


phy.abs <- function(dat, tmp, q, rtreephy, wk, formula){
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
  }
  else if(formula == "est"){
    #aa_joint <- c(aa) %>% `names<-`(rep(rownames(aa), N))
    qDa <- est.phy.q.joint(xx = aa, LL = gL, TT = TT, n = n, q, rtreephy)
    qDg <- est.phy.q(xx = ga, LL = gL, TT = TT, n = n, q, rtreephy)
  }
  if(q!=1) {
    alpha.C <- qDa^(1/(1-q))/N
    gamma.C <- qDg^(1/(1-q))
  } else {
    alpha.C <- exp(qDa/TT+log(TT))/N
    gamma.C <- exp(qDg/TT+log(TT))
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

    ##chunyu revise##
    ans = PD_obs + ifelse((2*f1)*g2>(g1*f2), (n-1)/n*g1^2/(2*g2), (n-1)/n*g1*(f1-1)/(2*(f2+1)))
    #ans = PD_obs + ifelse(g2>0, (n-1)/n*g1^2/(2*g2), (n-1)/n*g1*(f1-1)/2*(f2-1))

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


est.phy.q.joint <- function(xx, LL, TT, n, q, rtreephy){ # proposed
  N <- ncol(xx)

  #each assemblage calculate once f1, f2
  f12 <- sapply(1:N, function(k){
    xxk <- xx[,k]
    LLk <- LL[names(xxk)]
    tmp <- data.frame(abun = xxk, length = LLk)
    #f1 <- ifelse(datatype=='incidence', sum(rowSums(data)==1), sum(data==1) )
    f1 = sum(tmp[, 1]==1)
    f2 = sum(tmp[, 1]==2)
    node <- names(rtreephy$parts)
    node2 = names(xxk)[xxk==2][names(xxk)[xxk==2] %in% node]
    if(length(node2)>0){
      rep2 <- sapply(node2, function(tx){
        sum(xxk[rtreephy$parts[[tx]]] %in% c(2,0)) == length(rtreephy$parts[[tx]])
      })
      f2 <- f2 - sum(rep2)
    }
    node1 = names(xxk)[xxk==1][names(xxk)[xxk==1] %in% node]
    if(length(node1)>0){
      rep1 <- sapply(node1, function(tx){
        sum(xxk[rtreephy$parts[[tx]]] %in% c(1,0)) == length(rtreephy$parts[[tx]])
      })
      f1 <- f1 - sum(rep1)
    }
    return(c(f1,f2))
  }) %>% rowSums()

  xx <- c(xx) %>% `names<-`(rep(rownames(xx), N))
  LL <- rep(LL, N)
  LL <- LL[names(xx)]
  tmp <- data.frame(abun = xx, length = LL)
  PD_obs <- sum(tmp[tmp[,1]>0,2])
  g1 = sum(tmp[tmp[, 1]==1, 2])
  g2 = sum(tmp[tmp[, 1]==2, 2])
  f1 <- f12[1]; f2 <- f12[2]

  if(f2 > 0){
    A = 2*f2/((n-1)*f1+2*f2)
  }else if(f2 == 0 & f1 > 0){
    A = 2/((n-1)*(f1-1)+2)
  }else{
    A = 1
  }
  t_bar <- sum(xx*LL/n)

  if(q==0){

    ##chunyu revise##
    ans = PD_obs + ifelse((2*f1)*g2>(g1*f2), (n-1)/n*g1^2/(2*g2), (n-1)/n*g1*(f1-1)/(2*(f2+1)))
    #ans = PD_obs + ifelse(g2>0, (n-1)/n*g1^2/(2*g2), (n-1)/n*g1*(f1-1)/2*(f2-1))

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


est.spader <- function(xa, xg, LL, TT, nk, wk, q, rtreephy){
  Xip = xg
  Xik = xa
  Li = LL
  N = ncol(Xik)

  if (q == 1) {
    W = sum(-wk * log(wk))
    r.data = sapply(1:N, function(k) Xik[, k]/nk[k]) #pik
    r.pool = c(r.data %*% wk)
    U = numeric(N)
    K = numeric(N)
    for (k in 1:N) {
      I = which(Xik[, k] * (rowSums(Xik) - Xik[, k]) > 0)
      is = Xik[, k][I]
      L.is = Li[I]
      pools = rowSums(Xik)[I] - is
      r.is = is/nk[k]
      r.pools = r.pool[I]
      U1 = sum(L.is * r.is)
      sf1 = sum(pools == 1)
      sf2 = sum(pools == 2)
      sf2 = ifelse(sf2 == 0, 1, sf2)
      U2 = sum(L.is[pools == 1] * r.is[pools == 1]) * (sf1/(2 * sf2))
      U[k] = max(0, 1 - U1 - U2) * (-wk[k] * log(wk[k]))
      K[k] = -sum(L.is * wk[k] * r.is * log(r.pools/r.is))
    }
    est = (sum(U) + sum(K))/W
    C_1 = est
    U_1 = est
  }
  else if (q == 0){
    PDk <- sapply(1:N, function(k){
      get("est.phy.q")(xx = Xik[,k], LL = Li, TT = TT, n = nk[k], q, rtreephy)
    })
    PD <- est.phy.q(xx = Xip, LL = Li, TT = TT, n = sum(nk), q, rtreephy)
    C_1 = (PD-sum(wk*PDk))/sum((1-wk)*PDk)
    U_1 = (PD-sum(wk*PDk))/(1-sum(wk*PDk)/sum(PDk))/PD
  }
  else if (q == 2){
    qDk = sapply(1:N, function(k)
      sum(Li * Xik[, k] * (Xik[, k] - 1)/(nk[k] * (nk[k] - 1))))
    pik = sapply(1:N, function(k) Xik[,k]/nk[k])
    pik.1 = sapply(1:N, function(k) (Xik[,k]-1)/(nk[k]-1))
    temp = sapply(1:nrow(pik), function(i)
      (sum((pik[i, ] %*% t(pik[i, ]))*(wk %*% t(wk))) - sum((pik[i, ]*wk)^2)) + sum(pik[i, ]*pik.1[i, ]*wk^2))

    alpha.R <- (wk%*%qDk/TT^2)^(-1)
    joint.R <- (wk^2%*%qDk/TT^2)^(-1)
    gamma.R <- (sum(Li * temp)/TT^2)^(-1)

    beta.R <- gamma.R/alpha.R
    betamax.R <- joint.R/alpha.R

    C_1 = (beta.R^(-1)-1)/(betamax.R^(-1)-1)
    U_1 = (beta.R-1)/(betamax.R-1)
  }
  else {C_1 = NA; U_1 = NA}

  out = list(C_1 = C_1, U_1 = U_1)
  return(out)
}

#delta 有改
delta <- function(data, k, n, A){
  ans = sapply(1:length(k), function(i){
    if(k[i]<n){
      data1 = data[data[,1]<=(n-k[i]) &  data[,1]>=1,]
      if(inherits(data1, "numeric")) data1 = t(as.matrix(data1))
      sum( data1[,2]*(data1[,1]/n)*exp(lchoose(n-data1[,1], k[i])-lchoose(n-1, k[i])) )
    }else{
      g1 = sum(data[data==1,2])
      g1*(1-A)^(k[i]-n+1)/n
    }
  })
  return( ans )
}

# bootstrap.q.Beta = function(data, mat, rtree, tmp, q, nboot, wij, type, method){
#   H <- nrow(mat)
#   ##chunyu revise## number of value is 10*H-6
#   out = array(0, dim=c((10*H-6), length(q), nboot))
#   pool <- rowSums(data)
#   rtreephy <- newick2phylog(convToNewick(rtree))
#   #if(datatype == "abundance"){
#   n = colSums(data) ; N = ncol(data)
#   pop = Boots.pop(data, rtree, tmp$Alpha)  #回傳 L,p,unseen
#   S = nrow(data)
#   #B = length(c(phytree$leaves,phytree$parts))
#   if(pop$unseen == 0) p = pop$p[1:S,]
#   if(pop$unseen != 0) p = pop$p[c(1:S,tail(1:nrow(pop$p), pop$unseen)),]
#   boot.data = array(0, dim = dim(p))
#   rownames(boot.data) <- rownames(p)
#   S = nrow(data)
#   B = S+rtree$Nnode
#   for(i in 1:nboot){
#     L = pop$L
#     if(pop$unseen == 0) p = pop$p[1:S,]
#     if(pop$unseen != 0) p = pop$p[c(1:S,(B+1):nrow(pop$p)),]
#     boot.data = array(0, dim = dim(p))
#     for(j in 1:ncol(p)) boot.data[,j] = rmultinom(1,n[j],p[,j]) 
#     rownames(boot.data) <- rownames(p)
#     unseen = boot.data[-(1:S),]
#     ##chunyu revise##
#     if(is.vector(unseen)) unseen = matrix(unseen, nrow = 1) %>% `row.names<-`("u1")
#     boot.data.obs <- boot.data[1:S,]
#     #boot.data.obs <- boot.data.obs[rowSums(boot.data.obs)>0, ]
#     #tip.boot <- names(rtreephy$leaves)[!names(rtreephy$leaves)%in%rownames(boot.data.obs)]
#     #rtree.boot <- drop.tip(rtree, tip.boot)
#     #rtreephy.boot <- newick2phylog(convToNewick(rtree.boot))
#     boot.datatmp = apply(boot.data.obs, 2, function(x){
#       #names(x) = names(rtreephy$leaves)
#       tmp <- choose_data(x, rtreephy)
#       abun <- tmp[ ,1]
#       names(abun) <- rownames(tmp)
#       return(abun)
#     })  
#     boot.datatmp = rbind(boot.datatmp, unseen)
#     boot.gamma = rowSums(boot.datatmp)
#     #boot.gamma = boot.gamma[boot.gamma]
#     #boot.datatmp <- boot.datatmp[names(boot.gamma), ]
#     boot.gamma <- boot.gamma[boot.gamma>0]
#     boot.datatmp <- boot.datatmp[names(boot.gamma), ]
#     L <- L[names(boot.gamma), ]
#     L.gamma = rowSums(boot.datatmp[names(boot.gamma), ] * L) / boot.gamma
#     boot.alpha <- apply(boot.datatmp, 2, function(s) data.frame(branch_abun = s, branch_length = L.gamma))
#     boot.gamma <- data.frame(branch_abun = boot.gamma, branch_length = L.gamma)
#     boot.tmp <- list(Gamma = boot.gamma, Alpha = boot.alpha)
#     #sapply(q, function(qq) method(dat = boot.data, mat, boot.tmp, qq, rtreephy = rtreephy, wij, type))
#     out[,,i] = sapply(q, function(qq) method(dat = boot.data, mat, boot.tmp, qq, rtreephy = rtreephy, wij, type))
#   }
#   #print(sum(is.infinite(apply(out, 3, sum))))
#   #out[ , ,!is.infinite(apply(out, 3, sum))]
#   #}
#   return(out)
# }

bootstrap.q.Beta_OBS <- function(data, rtree, tmp, q, nboot, wk, formula, method){
  if (formula == "spader"){
    out = array(0, dim = c(2, length(q), nboot))
  }else out = array(0, dim = c(7, length(q), nboot))
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
    ##chunyu revise##
    if(is.vector(unseen)) unseen = matrix(unseen, nrow = 1) %>% `row.names<-`("u1")
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

# Boots.pop = function(data, rtree, tmp){
#   # if(datatype == "abundance"){
#   N = ncol(data); n = colSums(data)
#   pool=rowSums(data) ; OBS=length(pool)
#   rtreephy <- newick2phylog(convToNewick(rtree))
#   OBS_B <- dim(tmp[[1]])[1]  #這邊的 tmp 是傳tmp$alpha 近來
#   obs <- colSums(data>0)
#   TT <- sum(tmp[[1]][,1]/n[1]*tmp[[1]][,2])   #樹根的總長度
#   F1=sum(pool==1);F2=sum(pool==2)
#   F0=ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2))*(sum(n)-1)/sum(n)  #pool assemblage f0 estimate
#   F0_N <- round(F0)
#   f1=sapply(1:N,function(k) sum(data[,k]==1))
#   f2=sapply(1:N,function(k) sum(data[,k]==2))
#   g1=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==1,2]) ))
#   g2=unlist(lapply(tmp, function(tmp) sum(tmp[tmp[,1]==2,2]) ))
#   C <- ifelse(f2 == 0, C <- 1 - f1/n*(n-1)*(f1-1)/((n-1)*(f1-1) + 2), C <- 1 - f1/n*(n-1)*f1/((n-1)*f1 + 2*f2))
#   f0 <- ifelse(f2 == 0, f0 <- f1*(f1-1)/2, f0 <- f1^2/(2*f2))*(n-1)/n
#   f0_N <- round(f0)  #跟著第一層的跑
#   r.data=sapply(1:N,function(k) data[,k]/n[k]) #relative species abundance
#   W=sapply(1:N,function(k) (1-C[k])/sum(r.data[,k]*(1-r.data[,k])^n[k])) #N個lambda
#   g0 = sapply(1:N, function(k)  
#     if((2*g2[k]*f1[k])>g1[k]*f2[k]) (n[k]-1)/n[k]*g1[k]^2/2/g2[k]
#     else (n[k]-1)/n[k]*g1[k]*(f1[k]-1)/2/(f2[k]+1) )
#   if(F0>0){boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=F0_N))  
#   # 如果混合群落有未被觀察到的物種，就新增多少row上去
#   }else{boots.pop=r.data}
#   L = matrix(0, nrow=(OBS_B+F0_N), ncol=N) # 包含觀察和未觀察物種的L
#   boots.pop2 = matrix(0, nrow=(OBS_B+F0_N), ncol=N) 
#   #obs branch abundance + undetected豐度矩陣
#   
#   #針對每個樣本，調整相對豐度，生成自助重抽樣數據。
#   # 如果樣本中有未觀察物種，則進行相對豐度的均分。
#   # L：進化樹各節點的分支長度
#   for(i in 1:N)
#   { 
#     if(f0_N[i]>0)
#     {
#       f0_N[i]=ifelse(f0_N[i]+obs[i]>OBS+F0_N, OBS+F0_N-obs[i],f0_N[i])
#       boots.pop[,i][1:OBS] <- r.data[ ,i]*(1-W[i]*(1-r.data[ ,i])^n[i])   #dectected species 修正
#       I=which(boots.pop[,i]==0) #ai+bi
#       II=sample(I,f0_N[i])
#       u.p <- (1-C[i])/f0_N[i]
#       boots.pop[II,i]=rep(u.p,f0_N[i])
#       da = boots.pop[1:OBS,i] 
#       #corrected observed relative species abundance
#       names(da) = rownames(data)
#       mat = choose_data(da, rtreephy)  #mat[,1]:corrected observed branch_abun
#       #corrected observed relative node abundance
#       boots.pop2[,i]=c(mat[,1], boots.pop[,i][-(1:OBS)]) #add undetected species 
#       F00 = sum(II > OBS) #not detect in pool assemblage 
#       L[1:nrow(mat),i] = mat[,2] #mat[,2]:corrected observed branch_length
#       if(F00>0 ){
#         index = which(boots.pop2[,i] > 0)[which(boots.pop2[,i] > 0) > nrow(mat)]
#         #un.sp <- rownames(data)[II[II<OBS]]
#         #g0r <- g0[i]- sum(mat[un.sp,2])
#         L[index, i] = g0[i]/F00  #多層次論文中的公式p.44 (2.60)
#       }
#     }else{
#       
#       ##yayun revise##
#       da = boots.pop[1:OBS,i]
#       names(da) = rownames(data)
#       mat = hiDIP:::choose_data(da, rtreephy)
#       ###############
#       
#       L[seq_len(OBS_B), i] = tmp[[i]][ ,2] 
#       boots.pop2[seq_len(OBS_B), i] = tmp[[i]][ ,1]
#     }
#   }
#   if(F0_N==0){
#     rownames(L) <- rownames(mat)
#     rownames(boots.pop2) <- rownames(mat)
#   } else{
#     rownames(L) <- c(rownames(mat), paste0("u", seq_len(F0_N)))
#     rownames(boots.pop2) <- c(rownames(mat), paste0("u", seq_len(F0_N)))
#   }
#   L[L>TT] <- TT
#   return(list(p=boots.pop2,L=L,unseen=F0_N))
#   # }
# }


# transconf = function(Bresult, est, conf){
#   est.btse = apply(Bresult, 1, sd)
#   est.LCL = est - qnorm(1-(1-conf)/2) * est.btse
#   est.UCL = est + qnorm(1-(1-conf)/2) * est.btse
#   # if(any(est.LCL<0)) est.LCL[est.LCL<0] <- 0
#   # if(any(est.UCL>1)) est.UCL[est.UCL>1] <- 1
#   cbind(est = est, btse=est.btse, LCL=est.LCL, UCL = est.UCL)
# }



hier.phylogeny <- function(data, mat, tree, q = seq(0, 2, 0.2), weight = "size", nboot = 20,
                           conf = 0.95, type = "mle", decomposition = "relative"){
  dat <- data[rowSums(data)>0, ]
  if(decomposition == "relative"){
    method = phy.H.rel
  }else{
    method = phy.H.abs
  }
  H <- nrow(mat)
  rtip <- tree$tip.label[!tree$tip.label %in% rownames(dat)]
  #tree$tip.label[!rtree$tip.label %in% rownames(dat)]
  rtree <- drop.tip(tree, rtip)
  tmp <- TranMul(dat, rtree)
  rtreephy <- newick2phylog(convToNewick(rtree))
  if(inherits(weight, "numeric")){
    wij <- weight
  } else if (weight == "size"){
    wij <- colSums(dat)/sum(dat)
  } else  {
    wij <- rep(1/ncol(dat), ncol(dat))
  } 
  est <- sapply(q, function(i) method(dat, mat, tmp, i, rtreephy, wij, type))
  if(nboot!=0){
    H <- nrow(mat)
    boot.est <- bootstrap.q.Beta(data = dat, mat, rtree = rtree, tmp = tmp, q = q, nboot = nboot, wij = wij, type = type, method)
    test <- boot.est[seq_len(H+1), , ]
    #is.infinite(sum(boot.est))
    #test <- boot.est[head(seq_len(dim(boot.est)[1]),H),1,]
    id <- apply(test, 1:2, function(x) {
      bb <- x
      q1 <- quantile(bb,0.25)
      q3 <- quantile(bb,0.75)
      q1 <- q1-1.5*(q3-q1)
      q3 <- q1+1.5*(q3-q1)
      which(bb >= q1 & bb <= q3)
    })
    #執行 nboot 次 bootstrap 後，計算每個元素的四分位範圍。
    #然後通過四分位範圍篩選掉偏離範圍的數據，只保留那些在四分位範圍內的數據
    index <- Reduce(function(x,y) {intersect(x,y)}, id)
    boot.est <- boot.est[ ,,index]
    #boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ][boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ]<0] <- 0
    #boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ][boot.est[tail(seq_len(dim(boot.est)[1]),-2*H), ]>1] <- 1
    #dim(boot.est)
    #diff.boot.est <- as.data.frame(apply(boot.est,3, tail, n = H))
    out = lapply(seq_along(q), function(i){
      
      #boot.est[j,i,] 是去選假設qI_alpha1 有四次bootstrap的結果(transconf會去算 sd(Bresult))
      x = sapply(seq_len(nrow(boot.est)), function(j) transconf(Bresult = boot.est[j,i,], est = est[j,i], conf)) %>% t()
      colnames(x) = c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
      rownames(x) = rownames(est)
      return(x)}) %>% do.call(rbind,.)
    
    ##chunyu revise## number of value is nrow(est) = 10*H-6
    Order.q = rep(q, each = nrow(est))
    
    Method = rep(rownames(out)[1:nrow(est)], length(q))
    
    rownames(out) = NULL
    
    out = cbind(Method, Order.q, as.data.frame(out), Decomposition = decomposition)
    
    out
    # 
    # CL = t(sapply(seq_len(nrow(boot.est)), function(j) transconf(matrix(boot.est[j,,], nrow=length(q)), est[j,], conf)))
    # rownames(CL) <- rownames(est)
    # colnames(CL) <- c(sapply(paste0("q=",q), function(k) paste(k,c("est", "bt.sd", "LCL", "UCL"), sep = "_")))
    # 
  }else{
    out <- lapply(seq_along(q), function(i){x = data.frame("Estimator" = est[,i], Bootstraps.e = NA, LB = NA, UB = NA) 
    colnames(x) = c("Estimator", "Bootstrap S.E.", "LCL", "UCL")
    return(x)}) %>% do.call(rbind,.)
    
    Order.q = rep(q, each = nrow(est))
    
    ##chunyu revise## number of value is nrow(est) = 10*H-6
    Method = rep(rownames(out)[1:nrow(est)], length(q))
    
    rownames(out) = NULL
    
    out = cbind(Method, Order.q, as.data.frame(out), Decomposition = decomposition)
    out$'Bootstrap S.E.' = as.numeric(out$'Bootstrap S.E.')
    out$LCL = as.numeric(out$LCL)
    out$UCL = as.numeric(out$UCL)
    
    out
  }
  return(out)
}

phy.H.rel <- function(dat, mat, tmp, q, rtreephy, wij, type){
  n <- sum(dat)
  nsite <- ncol(mat)
  H <- nrow(mat)
  index <- lapply(rev(seq_len(H)), function(i) lapply(as.list(unique(mat[i,])), function(x) which(mat[i,] == as.character(x))))
  w <- lapply(index, function(x) sapply(x, function(y) sum(wij[y])))
  M <- lapply(index, function(x) sapply(x,length))
  ga <- tmp$Gamma[,1];gB <- tmp$Gamma[,2]
  names(gB) <- rownames(tmp$Gamma)
  gp=ga/n;TT=sum(gp*gB);
  aa <- sapply(tmp$Alpha, function(X){
    abun <- X[,1]
    names(abun) = rownames(X)
    return(abun)
  })
  B <-length(ga)
  pij <- sapply(1:nsite, function(x) {
    abun <- aa[,x]
    n <- sum(dat[,x])
    return(abun/n)
  }) #pk|im
  weight.p <- pij*t(replicate(B,wij))
  p <- lapply(index, function(h) sapply(h, function(x) rowSums(as.matrix(weight.p[,x])/sum(wij[x]))))
  x <- lapply(index, function(h) sapply(h, function(x) rowSums(as.matrix(dat[ ,x]))))
  xtmp <- lapply(index, function(h) sapply(h, function(x) rowSums(cbind(aa[ ,x]))))
  if(type == "mle"){
    qD <- sapply(1:H, function(i) apply(p[[i]], 2, get("mle.phy.q"), LL = gB, TT = TT, q)) ##\sum{p^q} or -\sum{p*log(p)})
  } else if(type == "est"){
    qD <- sapply(1:H, function(i) sapply(seq_len(ncol(x[[i]])), function(j){
      n <- sum(x[[i]][,j])
      get("est.phy.q")(xx = xtmp[[i]][,j], LL = gB, TT = TT, n = n, q = q, rtreephy)
    }))
  }
  alpha.relative <- sapply(1:H, function(i){
    qD_est <- qD[[i]]
    if(q!=1){
      (w[[i]]%*%qD_est)^(1/(1-q))
    } else{
      exp(w[[i]]%*%(qD_est/TT+log(TT)))
    }
  })
  names(alpha.relative) <- c(paste0("qPD_alpha",seq_len(H-1)),"qPD_gamma")
  all.pair <- cbind(sapply(1:(H-1), function(i) cbind(i,i+1)), c(1,H))
  Ialpha.relative <- sapply(alpha.relative ,function(x){
    if(q!=1){
      alpha.relative.H.j <- (TT-x^(1-q)*TT^q)/(q-1)
    } else{
      alpha.relative.H.j <- (log(x)-log(TT))*TT
    } 
    alpha.relative.H.j}
  )
  names(Ialpha.relative) <- c(paste0("qI_alpha", seq_len(H-1)), "qI_gamma")
  Ibeta.relative <- apply(all.pair, 2, function(j){
    alpha <- alpha.relative
    if(q!=1){
      alpha.relative.H.j <- (TT-alpha[j]^(1-q)*TT^q)/(q-1)
    } else{
      alpha.relative.H.j <- (log(alpha[j])-log(TT))*TT
    }
    diff(alpha.relative.H.j)
  })
  names(Ibeta.relative) <- c(paste0("qI_Beta ", "level(",apply(all.pair, 2,paste, collapse = "|"), ")"))
  differential.relative <-  apply(all.pair, 2, function(j){
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    num <- sapply(index[[j[1]]], function(a) seq_along(w.down)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    w.down.new <-w.down[num]
    ratio<- w.up/w.down.new
    est <- qD[[j[1]]]*TT^q
    Max <- ifelse(q!=1, sum(w.down.new*ratio*(1-ratio^(q-1))*est)/(q-1),
                  -TT*sum(w.up*log(ratio)))
    alpha <- alpha.relative
    if(q!=1){
      alpha.relative.H.j <- (TT-alpha[j]^(1-q)*TT^q)/(q-1)
    } else{
      alpha.relative.H.j <- (log(alpha[j])-log(TT))*TT
    }
    diff(alpha.relative.H.j)/Max
    #sum(w.down*sapply(index[[j[2]]], function(g) sum(ratio[g]*(1-ratio[g]^(q-1))*group.p[g])))
  })
  names(differential.relative ) <- paste0("Delta ","level(",apply(all.pair, 2,paste, collapse = "|"), ")")
  beta <- apply(all.pair, 2, function(j) alpha.relative[j[2]]/alpha.relative[j[1]])
  betamax <- apply(all.pair, 2, function(j){
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    num <- sapply(index[[j[1]]], function(a) seq_along(w.down)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    w.down.new <-w.down[num]
    ratio<- w.up/w.down.new
    if(q==1){
      out <- exp(-sum(w.up*log(ratio)))
    } else{
      out <- (sum(w.down.new*(ratio)^q*qD[[j[[1]]]])/sum(w.down.new*(ratio)*qD[[j[[1]]]]))^(1/(1-q))
    }
    out
  })
  # gammamax <- betamax*apply(all.pair, 2, function(j){
  #   alpha.relative[j[1]]
  # })
  # gamma <- apply(all.pair, 2, function(j){''
  #   alpha.relative[j[2]]
  # })
  beta <-  beta[1:(H-1)]
  betamax <- betamax[1:(H-1)] 
  # gamma <-  gamma[1:(H-1)]
  # gammamax <- gammamax[1:(H-1)] 
  ifelse(q==1, C_1 <- log(beta)/log(betamax), C_1 <- (beta^(1-q)-1)/(betamax^(1-q)-1))
  ifelse(q==1, U_1 <- log(beta)/log(betamax), U_1 <- (beta^(q-1)-1)/(betamax^(q-1)-1))
  V_1 <- (beta-1)/(betamax-1)
  S_1 <- (beta^-1-1)/(betamax^-1-1)
  diff <- c("1-CqN" = C_1, "1-UqN" = U_1, "1-VqN" = V_1, "1-SqN" = S_1)
  names(diff) <- c(sapply(c("1-CqN", "1-UqN", "1-VqN", "1-SqN"), function(y) paste0(y,"(",apply(all.pair[,1:(H-1)], 2,paste, collapse = "|"), ")")))
  names(beta) <- paste0("qPD_Beta ",seq_len(H-1))
  names(betamax) <- paste0("qPD_Beta_max ",seq_len(H-1))
  out <-  c(Ialpha.relative ,Ibeta.relative, differential.relative, alpha.relative, beta, betamax, diff)
  #out <- data.frame(out)
  return(out)
}

phy.H.abs <- function(dat, mat, tmp, q, rtreephy, wij, type){
  n <- sum(dat)
  nsite <- ncol(mat)
  H <- nrow(mat)
  index <- lapply(rev(seq_len(H)), function(i) lapply(as.list(unique(mat[i,])), function(x) which(mat[i,] == as.character(x))))
  
  #rownames(aa) <- rownames(tmp$Gamma)
  
  # alpha.relative <- numeric(H)
  # alpha.absolute <- numeric(H)
  #wij <- colSums(dat)/n
  w <- lapply(index, function(x) sapply(x, function(y) sum(wij[y])))
  M <- lapply(index, function(x) sapply(x, length))
  ga <- tmp$Gamma[,1]; gB <- tmp$Gamma[,2]
  names(gB) <- rownames(tmp$Gamma)
  gp=ga/n #pi|++
  TT=sum(gp*gB)
  aa <- sapply(tmp$Alpha, function(X){
    abun <- X[,1]
    names(abun) = rownames(X)
    return(abun)
  })
  B <- length(ga)
  pijk <- sapply(1:nsite, function(x) {
    abun <- aa[,x]
    n <- sum(dat[,x])
    return(abun/n)
  }) #pi|jk
  weight.p <- pijk*t(replicate(B,wij)) #wjk*pi|jk
  p <- lapply(index, function(h) sapply(h, function(x) rowSums(as.matrix(weight.p[,x])/sum(wij[x]))))
  #x <- lapply(index, function(h) sapply(h, function(x) rowSums(as.matrix(dat[ ,x]))))
  nlist <- lapply(index, function(h) sapply(h, function(x) sum(dat[ ,x])))
  xtmp <- lapply(index, function(h) sapply(h, function(x) rowSums(cbind(aa[ ,x]))))
  if(type == "mle"){
    qD <- sapply(1:H, function(i) apply(p[[i]], 2, get("mle.phy.q"), LL = gB, TT = TT, q)) ##\sum{p^q} or -\sum{p*log(p)})
  } else if(type == "est"){
    qD <- sapply(1:H, function(i) sapply(seq_along(M[[i]]), function(j){
      n <- nlist[[i]][j]
      get("est.phy.q")(xtmp[[i]][,j], LL = gB, TT = TT, n = n, q = q, rtreephy)
    }))
  }
  alpha.relative <- sapply(1:H, function(i){
    qD_est <- qD[[i]]
    if(q!=1){
      (sum( (nlist[[i]]/n)^q * M[[i]]^(1-q)*qD_est )^(1/(1-q)))/nsite #revised
      #(w[[i]]%*%qD_est)^(1/(1-q))
    } else{
      (exp( -sum(nlist[[i]]/n*(log(nlist[[i]]/n/TT/M[[i]])-qD_est/TT)) ))/nsite #revised
      #exp(w[[i]]%*%(qD_est/TT+log(TT)))
    }
  })
  #names(alpha.absolute) <- c(paste0("alpha",seq_len(H-1)),"gamma")
  names(alpha.relative) <- c(paste0("qPD_alpha",seq_len(H-1)),"qPD_gamma")
  all.pair <- cbind(sapply(1:(H-1), function(i) cbind(i,i+1)), c(1,H))
  Ialpha.relative <- sapply(alpha.relative ,function(x){
    if(q!=1){
      alpha.relative.H.j <- (TT-x^(1-q)*TT^q*(nsite^(1-q)))/(q-1)
    } else{
      alpha.relative.H.j <- (log(x)-log(TT))*TT
    }
    alpha.relative.H.j})
  names(Ialpha.relative) <- c(paste0("qI_alpha", seq_len(H-1)), "qI_gamma")
  Ibeta.relative <- apply(all.pair, 2, function(j){
    alpha <- alpha.relative
    if(q!=1){
      alpha.relative.H.j <- (TT-alpha[j]^(1-q)*TT^q*(nsite^(1-q)))/(q-1)
    } else{
      alpha.relative.H.j <- (log(alpha[j])-log(TT))*TT
    }
    diff(alpha.relative.H.j)
  })
  names(Ibeta.relative) <- c(paste0("qI_Beta ", "level(",apply(all.pair, 2,paste, collapse = "|"), ")"))
  differential.relative <-  apply(all.pair, 2, function(j){
    #qD = get(paste0("qD.", Q))
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    M1 <- M[[j[[1]]]]
    M2 <- M[[j[[2]]]]
    num <- sapply(index[[j[1]]], function(a) seq_along(M2)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    Max <- ifelse(q!=1, sum((M1^(1-q)-M2[num]^(1-q))*(nlist[[j[[1]]]]/n)^q*TT^q*qD[[j[[1]]]])/(q-1),
                  TT*(sum(w.down*log(M2))-sum(w.up*log(M1))))
    alpha <- alpha.relative
    if(q!=1){
      alpha.relative.H.j <- (TT-alpha[j]^(1-q)*TT^q*(nsite^(1-q)))/(q-1)
    } else{
      alpha.relative.H.j <- (log(alpha[j])-log(TT))*TT
    }
    diff(alpha.relative.H.j)/Max
    #sum(w.down*sapply(index[[j[2]]], function(g) sum(ratio[g]*(1-ratio[g]^(q-1))*group.p[g])))
  })
  names(differential.relative ) <- paste0("Delta ","level(",apply(all.pair, 2,paste, collapse = "|"), ")")
  beta <- apply(all.pair, 2, function(j) alpha.relative[j[2]]/alpha.relative[j[1]])
  beta.max <- apply(all.pair, 2, function(j){
    w.up <- unlist(w[[j[1]]])
    w.down <- unlist(w[[j[2]]])
    M1 <- M[[j[[1]]]]
    M2 <- M[[j[[2]]]]
    num <- sapply(index[[j[1]]], function(a) seq_along(M2)[sapply(index[[j[2]]], function(b) all(a%in%b))])
    ifelse(q!=1,
           (sum(M2[num]^(1-q)*((nlist[[j[[1]]]]/n))^q*qD[[j[[1]]]])/sum(M1^(1-q)*((nlist[[j[[1]]]]/n))^q*qD[[j[[1]]]]))^(1/(1-q)),
           exp(sum(w.down*log(M2)))/exp(sum(w.up*log(M1))))
  })
  beta <-  beta[1:(H-1)]
  beta.max <- beta.max[1:(H-1)]
  ifelse(q==1, C_1 <- log(beta)/log(beta.max), C_1 <- (beta^(1-q)-1)/(beta.max^(1-q)-1))
  ifelse(q==1, U_1 <- log(beta)/log(beta.max), U_1 <- (beta^(q-1)-1)/(beta.max^(q-1)-1))
  V_1 <- (beta-1)/(beta.max-1)
  S_1 <- (beta^-1-1)/(beta.max^-1-1)
  diff <- c("1-C_qN" = C_1, "1-U_qN" = U_1, "1-V_qN" = V_1, "1-S_qN" = S_1)
  names(diff) <- c(sapply(c("1-CqN", "1-UqN", "1-VqN", "1-SqN"), function(y) paste0(y,"(",apply(all.pair[,1:(H-1)], 2,paste, collapse = "|"), ")")))
  names(beta) <- paste0("qPD_Beta ",seq_len(H-1))
  names(beta.max) <- paste0("qPD_Beta_max ",seq_len(H-1))
  #names(diff) <- c(sapply(c("1-C", "1-U", "1-V", "1-S"), function(y) withMathJax(paste0("$$",y,"_{",apply(all.pair, 2,paste, collapse = "|"), "}$$"))))
  #names(diff) <- c(sapply(c("1-C", "1-U", "1-V", "1-S"), function(y) paste0("$",y,"{",apply(all.pair, 2,paste, collapse = "|"), "}$")))
  #names(out) <- quote(names(out))
  #names(out) <- c(paste0("alpha of level ",seq_len(H-1)),"gamma", paste0("beta ",beta.pair[1, ]), paste0("differentiation ","level",apply(all.pair, 2,paste, collapse = "|")))
  out <-  c(Ialpha.relative ,Ibeta.relative, differential.relative, alpha.relative, beta, beta.max, diff)
  return(out)
}
gghier_phylogeny = function(outcome, method = 1){
  if(method == 1){
    outcome = outcome[c(grep(c("qI_alpha"), outcome$Method),
                        grep(c("qI_gamma"), outcome$Method)) %>% sort(),]
  }else if(method == 2){
    outcome = outcome[grep("qI_Beta level", outcome$Method),]
  }else if(method == 3){
    outcome = outcome[grep("Delta level", outcome$Method),]
  }else if(method == 4){
    outcome = outcome[c(grep(c("qPD_alpha"), outcome$Method),
                        grep(c("qPD_gamma"), outcome$Method)) %>% sort(),]
  }else if(method == 5){
    outcome = outcome[grep("qPD_Beta ", outcome$Method),]
  }else if(method == 6){
    outcome = outcome[grep("1-", outcome$Method),]
  }
  if(method != 6){
    out = ggplot(outcome, aes(x = Order.q, y = Estimator, colour = Method, 
                              fill = Method))+
      geom_line(size = 1.5) + 
      geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Method), linetype = 0, 
                  alpha = 0.2)
  }else{
    outcome$group = "1-C"
    outcome$group[grep("1-UqN", outcome$Method)] = "1-U"
    outcome$group[grep("1-SqN", outcome$Method)] = "1-S"
    outcome$group[grep("1-VqN", outcome$Method)] = "1-V"
    out = ggplot(outcome, aes(x = Order.q, y = Estimator, colour = Method, 
                              fill = Method))+
      facet_grid(group~.)+
      geom_line(size = 1.5) + geom_ribbon(aes(ymin = LCL, 
                                              ymax = UCL, fill = Method), linetype = 0, 
                                          alpha = 0.2)
  }
  
  out = out  + theme_bw() + 
    theme(legend.position = "bottom",legend.box = "vertical", 
          legend.key.width = unit(1.2, "cm"), 
          legend.title = element_blank(), 
          legend.margin = margin(0, 0, 0, 0), 
          legend.box.margin = margin(-10, -10, -5, -10), 
          text = element_text(size = 16), 
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) + 
    guides(linetype = guide_legend(keywidth = 2.5))
  return(out)
}