library("limSolve")
library("XML")
library("rlist")


rn.xmlextract <- function(f, verbose=F) {
  p <- xmlRoot(xmlParse(f))[["model"]]
  sp.id <- NULL
  sp.name <- NULL
  if ("listOfSpecies" %in% names(p)) {
    sp.l <- xmlApply(p[["listOfSpecies"]],xmlAttrs)
    sp.id <- as.vector(sapply(sp.l,function(e) e["id"]))
    sp.name <- as.vector(sapply(sp.l,function(e) e["name"]))
  }
  sp.idn <- sp.id
  i <- which(grepl("^S[0-9]+$",sp.idn))
  sp.idn[i] <- sp.name[i]
  p <- p[["listOfReactions"]]
  n <- which(xmlSApply(p, xmlName) == "reaction")
  a <- xmlApply(p, xmlAttrs)
  l <- list()
  n <- which(xmlApply(p, xmlName) == "reaction")
  R <- list()
  
  for (i in n) {
    reversible <- a[[i]]["reversible"]
    reversible <- if (!is.na(reversible) && reversible=="false") F else T
    reactants <- if("listOfReactants" %in% names(p[[i]])) xmlApply(p[[i]][["listOfReactants"]], xmlAttrs) else NULL
    products <- if("listOfProducts" %in% names(p[[i]])) xmlApply(p[[i]][["listOfProducts"]], xmlAttrs) else NULL
    r.stoichiometry <- sapply(reactants,function(e) { s <- e["stoichiometry"]; if (is.na(s)) 1 else s } )
    p.stoichiometry <- sapply(products,function(e) { s <- e["stoichiometry"]; if (is.na(s)) 1 else s } )
    names(r.stoichiometry) <- NULL; names(p.stoichiometry) <- NULL
    reactants <- if(is.null(reactants)) NULL else sapply(reactants,function(e) e["species"])
    products <- if(is.null(products)) NULL else sapply(products,function(e) e["species"])
    names(reactants) <- NULL; names(products) <- NULL
    R <- c(R,list(list( reactants=reactants,
                        products=products,
                        r.stoichiometry=as.numeric(r.stoichiometry),
                        p.stoichiometry=as.numeric(p.stoichiometry) )))
    if (reversible)
      R <- c(R,list(list( reactants=products,
                          products=reactants,
                          r.stoichiometry=as.numeric(p.stoichiometry),
                          p.stoichiometry=as.numeric(r.stoichiometry) )))
    if (verbose) {
      print("----------------------------")
      print(reversible)
      print(c(reactants,"-->",products))
      print(c(r.stoichiometry,"-->",p.stoichiometry))
    }
  }
  list(sp.id=sp.id,sp.name=sp.name,sp.idn=sp.idn,reac=R)
}

rn.proc <- function(f=rn.test,tot_org=T) {
  if (is.numeric(f)) {
    f <- paste0("./ReacNet/",dir("ReacNet","*.xml")[f])
  } 
  cat("f =",f,"\n")
  rn <- rn.xmlextract(f)
  scod <- 1:length(rn$sp.id); names(scod) <- rn$sp.id
  mr <- mp <- matrix(0,length(scod),length(rn$reac))
  rownames(mr) <- rn$sp.id; rownames(mp) <- rn$sp.id
  
  
  i <- 1
  for (r in rn$reac) {
    l <- length(r$reactants); lp <- length(r$products)
    j <- 1; for (s in r$reactants) { mr[s,i] <- r$r.stoichiometry[j]; j <- j + 1 } 
    j <- 1; for (s in r$products) { mp[s,i] <- r$p.stoichiometry[j]; j <- j + 1 } 
    i <- i + 1
  }
  cat(nrow(mr),"especies,",ncol(mr), "reacciones\n")
  
  
  if(!tot_org){
    rn <- c(rn,list(spc=scod,mr=mr,mp=mp))
    nsp <- rn.linp_org(rn,rn$sp.id,F)
    rn$spc=scod
    rn$mr=mr
    rn$mp=mp
    rn$nsp=nsp
    rn <<- rn
    return(rn)
  }
  
  m=mp-mr
  
  rsp <- which(sapply(1:length(m[,1]),function(i) any(m[i,]!=0))) 
  csp <- which(sapply(1:length(m[,1]),function(i) all(m[i,]<=0))) 
  rcsp <- intersect(rsp,csp)
  if(length(rcsp!=0)){
    cm <- matrix(0,length(scod),length(rcsp))
    mr <- cbind(mr,cm)
    
    for(i in 1:length(rcsp)){
      cm[rcsp[i],i]=1
      rn$reac <- list.append(rn$reac,list(reactants=NULL,products=rcsp[i],r.stoichiometry=numeric(0),p.stoichiometry=1))
    } 
    mp <- cbind(mp,cm)
    
    
  }
  
  rn$spc=scod
  rn$mr=mr
  rn$mp=mp
  nsp <- rn.linp_org(rn,rn$sp.id,F)
  
  if(length(nsp)==0) {
    rn <<- rn
    return(rn)
  }
  
  cm <- matrix(0,length(scod),length(nsp))
  rn$mr <- cbind(rn$mr,cm)
  for(i in 1:length(nsp)){ 
    cm[nsp[i],i]=1
    rn$reac <- list.append(rn$reac,list(reactants=NULL,products=nsp[i],r.stoichiometry=numeric(0),p.stoichiometry=1))
  }
  rn$mp <- cbind(rn$mp,cm)
  
  rn <<- rn
  return(rn)
}


rn.linp_org <- function(crn,id=crn$sp.id,verbose=T) {
  Ns <- nrow(crn$mr); Nr <- ncol(crn$mr)
  E <- cbind(diag(Ns),-diag(Ns),crn$mp-crn$mr)
  f <- rep(0,Ns)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)); h[1:(2*Ns)] <- 0
  Cost <- rep(0,ncol(E)); Cost[1:Ns] <- 1
  rn.linp.r <<- linp(E,f,G,h,Cost)
  X <- rn.linp.r$X
  k <- which(X[(1:Ns)]>0)
  out <- k
  if (verbose) {
    cat("needed species to be an organization: ",id[k],"\n")
    k <- which(X[(1:Ns)+Ns]>0)
    if (length(k)>0) cat("aditional overprodced species: ",id[k],"\n")
  }
  return(out)
}

rn.linp_sp <- function(crn,id=crn$sp.id,sp.add.i=integer(0),verbose=F) {
  Ns <- nrow(crn$mr); Nr <- ncol(crn$mr)
  E <- cbind(diag(Ns),crn$mp-crn$mr) 
  f <- rep(1,Ns)
  Cost <- rep(0,ncol(E)); Cost[setdiff(1:Ns,sp.add.i)] <- 1 
  rn.linp.r <- linp(E,f,Cost=Cost)
  rn.linp.r <<- rn.linp.r
  X <- rn.linp.r$X
  r <- union(which(X[1:Ns]<1),sp.add.i)
  if (verbose) {
    cat("overproduced species: ",id[r],"\n")
  }
  r
}


rn.linp_org_ov <- function(crn,id=crn$sp.id,opsp,verbose=T) {
  Ns <- nrow(crn$mr); Nr <- ncol(crn$mr)
  
  O <- matrix(0,nrow=length(id),ncol=length(opsp))
  k=1
  for(i in opsp){ 
    O[i,k] <- -1
    k <- k+1
  }
  E <- cbind(O,crn$mp-crn$mr)
  f <- rep(0,nrow(E))
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E))
  Cost <- rep(1.0,ncol(E));
  rn.linp.r <<- linp(E,f,G,h,Cost)
  
  X <- rn.linp.r$X[-(1:length(opsp))]
  if(rn.linp.r$IsError) X <- NULL
  
  return(X)
}


rn.sprod <- function(crn,i=1:nrow(crn$mr),sp.i=integer(0),verbose=F) {
  Ns <- length(i); Nr <- ncol(crn$mr[i,])
  rc <- c()
  for(j in 1:ncol(crn$mr)) if(all(crn$mp[-i,j]==0) && all(crn$mr[-i,j]==0)) rc <- c(rc,j) 
  sp <- rep(F,Ns)
  if (length(sp.i)>0) {
    Esp <- numeric(Ns); Esp[sp.i] <- 1;
    sp[sp.i] <- T
  }
  else Esp <- NULL
  E0 <- cbind(-diag(Ns),Esp,crn$mp[i,rc]-crn$mr[i,rc])
  f0 <- rep(0,Ns)
  Cost0 <- rep(0,ncol(E0))
  for (s in which(!sp)) {
    if (verbose) cat("considered species :",crn$sp.id[s],"...")
    E <- E0; E[s,s] <- 1
    f <- f0; f[s] <- 1
    Cost <- Cost0; Cost[s] <- 1
    r <- linp(E,f,Cost=Cost)
    if (r$solutionNorm==0) sp[s] <- T
    if (verbose) cat("cost:",r$solutionNorm,"\n")
    k <- which(r$X[1:Ns]>0)
    k <- k[k!=s]
    if (verbose && length(k)>0) cat("        overproduced:",crn$sp.id[k],"\n")
  }
  rn.linp.r <<- r
  if(any(sp))  return(i[which(sp)])
  else return(numeric(0))
}


rn.dcom <- function(crn,sp){
  opsp <- rn.sprod(crn,sp)
  
  c_m <- crn$mr[sp,]
  nc_m <- c_m
  for(i in 1:ncol(crn$mr)){ 
    c_m[,i] <- (crn$mp[sp,i]!=0 & crn$mr[sp,i]!=0) & (crn$mp[sp,i]==crn$mr[sp,i])
    nc_m[,i] <-(crn$mp[sp,i]!=crn$mr[sp,i])
  }
  csp <- sp[rowSums(nc_m)==0 & rowSums(c_m)!=0]
  
  fsp <- setdiff(sp,union(opsp,csp))
  rc <- c()
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  if(length(rc)==0) return()
  m <- crn$mp[fsp,rc,drop=F]-crn$mr[fsp,rc,drop=F]
  adj <- matrix(0,length(fsp),length(fsp))
  if(length(fsp)>0) for(i in 1:ncol(m)){
    k.a <- which(m[,i]!=0)
    adj[k.a,k.a] <- 1
    }
  adj <- adj+t(adj)
  adj[col(adj)==row(adj)] <-0 
  adj <- adj>0
  v <- fsp*0
  eci <- 0
  while(!is.na(ec <- match(0,v))){
    j <- 1
    while(T){
      ec<-c(ec,setdiff(which(adj[,ec[j]]),ec))
      j <- j+1
      if(j > length(ec)) break
    }
    eci <- eci+1
    v[ec] <- eci
  }
  return(list(opsp=opsp,csp=csp,fsp=fsp,ec=v))
}


rn.dcom.opsp <- function(crn,sp,opsp){
  
  c_m <- crn$mr[sp,]
  nc_m <- c_m
  for(i in 1:ncol(crn$mr)){ 
    c_m[,i] <- (crn$mp[sp,i]!=0 & crn$mr[sp,i]!=0) & (crn$mp[sp,i]==crn$mr[sp,i])
    nc_m[,i] <-(crn$mp[sp,i]!=crn$mr[sp,i])
  }
  csp <- sp[rowSums(nc_m)==0 & rowSums(c_m)!=0]
  
  fsp <- setdiff(sp,union(opsp,csp))
  rc <- c()
  
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  if(length(rc)==0) return()
  m <- crn$mp[fsp,rc,drop=F]-crn$mr[fsp,rc,drop=F]
  adj <- matrix(0,length(fsp),length(fsp))
  if(length(fsp)>0) for(i in 1:ncol(m)){
    k.a <- which(m[,i]!=0)
    adj[k.a,k.a] <- 1
  
  }
  adj <- adj+t(adj)
  adj[col(adj)==row(adj)] <-0 
  adj <- adj>0
  v <- fsp*0
  eci <- 0
  while(!is.na(ec <- match(0,v))){
    j <- 1
    while(T){
      ec<-c(ec,setdiff(which(adj[,ec[j]]),ec))
      j <- j+1
      if(j > length(ec)) break
    }
    eci <- eci+1
    v[ec] <- eci
  }
  return(list(opsp=opsp,csp=csp,fsp=fsp,ec=v))
}

# discrete dynamic concentration calculation function
rn.dina_step_stoi <- function(crn,sp,x,v,n){
  
  rc <- c()
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  P_M <- matrix(0,ncol=(n),nrow=length(x))
  C_M <- matrix(0,ncol=(n+1),nrow=length(x))
  V_M <- matrix(0,ncol=(n+1),nrow=length(v))
  R_L <- list()
  rownames(C_M) <- crn$sp.id[sp]
  rownames(P_M) <- crn$sp.id[sp]
  if(length(rc)==0){
    sapply(1:(n+1),function(i) C_M[,i] <- x)
    R_L <- lapply(1:n,function (x) return(1:ncol(crn$mr)))
    return(list(C_M=C_M,R_L=R_L))
  }
  
  S <- crn$mp[sp,rc,drop=F]-crn$mr[sp,rc,drop=F]
  Sr <- crn$mr[sp,rc,drop=F]
  R_L <<- list()
  
  vp <- function(x_i){
    
    asp <- which(Sr%*%v > x_i)
    if(length(asp)>0){
      
      r_p <- x_i/Sr%*%v
      p <- rep(1,length(v))
      ind <- c()
      for (i in 1:ncol(S)){
        nsp <- which(Sr[,i]!=0)  
        if(length(nsp)!=0) p[i] <- max(0,min(r_p[nsp]))
        
      }
      R_L <- c(R_L,list(which(p==0)))
      vo <- v
      vo[p!=1] <- v[p!=1]*p[p!=1]
      return(vo)
    }
    else{return(v)}
  }
  for (i in 1:(n+1)){
    
    if(i==1){ C_M[,i] <- x; V_M[,i] <- v}
    else{ 
      C_M[,i] <- round(sapply(S%*%vp(C_M[,i-1]) + C_M[,i-1],function (x) max(0,x)),2)
      V_M[,i] <- round(vp(C_M[,i-1]),2)
      P_M[,i-1] <- round(S%*%vp(C_M[,i-1]),2)
    }
  }
  return(list(C_M=C_M,V_M=V_M,P_M=P_M))
}

rn.dcom2sp <- function(rn,decom){
  opsp <- rn$sp.id[decom$opsp]
  csp <- rn$sp.id[decom$csp]
  fcs <- list()
  for(i in unique(decom$ec))
  {
    fcsp <- decom$fsp[decom$ec==i]
    fcs <- c(fcs,list(rn$sp.id[fcsp]))
  }
  return(list(opsp=opsp,fcs=fcs,csp=csp))
} 

