#' ---
#' title: "project source code: Community detection via bayesian approach"
#' author: Lei
#' description: the implementation of Viariational Bayesian estimation follows "Jake M. Hofman, Chris H. Wiggins, A Bayesian Approach to Network Modularity, arXiv:0709.3512v3, 2008."
#' ---

##generate random graph using SBM model 
#pi~the vector of probabilities for community assignment
library(Matrix)
sbmrnd=function(N, pintra, pinter, pi)
{
  K=length(pi)
  assignment=sample(1:K,size=N,replace=T,prob=pi)
  Moduleassign=Matrix(0, nrow=N, ncol=K,sparse=TRUE)
  Moduleassign[cbind(1:N,assignment)]=1
  #print(as.matrix(Moduleassign))
  ##i and j in the same community with zij==1
  a=which(Moduleassign!=0,arr.ind=T)
  a=a[,2]
  Z=matrix(rep(a,N),ncol=N)==matrix(rep(a,N),ncol=N,byrow=TRUE)
  #print(as.matrix(Z))
  m1=matrix(runif(N*N),nrow=N)<pintra*Z+pinter*(1-Z)
  m1[lower.tri(m1, diag = TRUE)]=0
  Ajancencymatrix= Matrix(m1|t(m1), sparse =TRUE)
  Massign=sparseMatrix(1:N,a,x=1)
  return(list(A=Ajancencymatrix, MA=Massign))
}

tr=function(A)
  return(sum(diag(A)))

#Variational Bayesian Approach
##input data, module number K
VBestimate=function(input,K,priors=NULL, maxiteration=15, tolerance=1e-4) {
  #K-number of communities
  #N-number of nodes
  #Q-mean filed parameters Q[i,j]
  #a, b - distribution hyperparameter for p_intra 
  #c, d - distribution hyperpaprameter for p_inter
  # n - hyperparameter for pi
  
  #initilization
  N=dim(input)[1]
  A=input
  if(!is.null(priors)) {
    a0=priors$a0
    b0=priors$b0
    c0=priors$c0
    d0=priors$d0
    n0=priors$n0
  }
  else {
    a0=2
    b0=1
    c0=1
    d0=2
    n0=matrix(1,1,K)
  }
  
  a=a0
  b=b0
  c=c0
  d=d0
  n=n0
  
  Q=matrix(runif(N*K),ncol = K,nrow = N)
  Q=Q/rowSums(Q)
  
  Fmeasure=numeric(maxiteration)
  #Iteration for variational Bayesian
  for(r in 1:maxiteration) {
    ##p_intra~ Beta(a,b)
    ##p_inter~Beta(c,d)
    
    ##Variational E-step
    ###Expecation of ln(Beta(a,b)) and ln(Beta(c,d))
    psia=psigamma(a)
    psib=psigamma(b)
    psic=psigamma(c)
    psid=psigamma(d)
    psiab=psigamma(a+b)
    psicd=psigamma(c+d)
    ###Expectation of ln(Dirichlet(n)) n-vector
    nhmu=-psigamma(sum(n))+psigamma(n)
    JG=psid-psicd-psib+psiab
    JL=psia-psib-psic+psid
  
    ###Evaluation expectation of Q distribution
    for(i in 1:N) {
      Q[i,]=0
      sQ=colSums(Q)
      lnQ=JL*(A[i,]%*%Q)-JG*sQ+nhmu
      qrowi=exp(lnQ-max(lnQ))
      Q[i,]=qrowi/sum(qrowi)
    }
    
    ##Variational M-step
    M=0.5*sum(A)
    C=N*(N-1)/2
    nc=colSums(Q)
    ta=0.5*tr(t(Q)%*%A%*%Q)
    temp=matrix(rep(nc,N),nrow = N, byrow=T)
    tb=0.5*tr(t(Q)%*%(temp-Q))-ta
    tc=M-ta
    td=C-M-tb
    ##update hyperparameter
    a=ta+a0
    b=tb+b0
    c=tc+c0
    d=td+d0
    n=nc+n0
    
    #compute approximation to the negative loglikelihood
    Qnz=Q[abs(Q)>1.1102e-16]
    Fmeasure[r]=lbeta(a,b)-lbeta(a0,b0)+lbeta(c,d)-lbeta(c0,d0)+
      sum(lgamma(n))-lgamma(sum(n))-
      (sum(lgamma(n0))-lgamma(sum(n0)))-
      sum(Qnz*log(Qnz))
    Fmeasure[r]=-Fmeasure[r]
    
    if(r>1) {
      dF=Fmeasure[r]-Fmeasure[r-1];
      if(dF>tolerance)
        warning('F increased by ', dF)
      if(abs(dF)<tolerance)
        break
    }
    
  }
  
  return(list(Fmeasure=Fmeasure[r], Q=Q, K=K, pintra_a=a, pintra_b=b, pinter_c=c, pinter_d=d, pin=n, time=Sys.time()))
}

#Running multiple times with different initial values to get the best
Vbmestimate=function(input,Kvec,priors=NULL, rounds=5, maxiteration=20) {
  bestresultlists=list()
  bestFlist=numeric(length(Kvec))
  for(k in 1:length(Kvec)) {
    Flist=numeric(rounds)
    resultlists=list()
    for(i in 1:rounds) {
      result=VBestimate(input,Kvec[k],priors,maxiteration)
      resultlists[[i]]=result
      Flist[i]=result$Fmeasure
      cat(Kvec[k],'round ',i, '\n')
    }
    mini=which(Flist==min(Flist))
    mini=mini[1]
    bestresultlists[[k]]=resultlists[[mini]]
    bestFlist[k]=bestresultlists[[k]]$Fmeasure
  }
  ##find best K
  bestkind=which(bestFlist==min(bestFlist))
  return(list(bestforek=bestresultlists, bestkindex=bestkind[1]))
}
  
#Dataset experiment

##real dataset
library(igraph)
library(Matrix)
links = read.table("facebook_combined.txt") 
el=as.matrix(links);
el[,1]=as.character(el[,1])
el[,2]=as.character(el[,2])
net=graph.edgelist(el,directed=FALSE)
A=get.adjacency(net,type="both")
#plot(net,vertex.label=NA)
image(A)

Kvec=seq(10,50,by=10)
start <- Sys.time ()
result=Vbmestimate(A, Kvec)
Sys.time () - start
finalr=result$bestforek[[result$bestkindex]]
finalr$K

Fkresult=numeric(length(Kvec))
for(i in 1:length(Kvec)){
  Fkresult[i]=(result$bestforek[[i]])$Fmeasure
}
plot(Kvec, Fkresult, xlab ='k', ylab = '-L(q)', type="b" )

##synthesized data
library(Matrix)
sbmdata=sbmrnd(1000,0.8,0.2,c(0.2,0.3,0.1,0.2,0.2))
as.matrix(sbmdata$A)
as.matrix(sbmdata$MA)
image(sbmdata$A)

Kvec=c(2:8)
start <- Sys.time ()
result=Vbmestimate(sbmdata$A, Kvec)
Sys.time () - start
finalr=result$bestforek[[result$bestkindex]]
finalr$K

Fkresult=numeric(length(Kvec))
for(i in 1:length(Kvec)){
  Fkresult[i]=(result$bestforek[[i]])$Fmeasure
}
plot(Kvec, Fkresult, xlab ='k', ylab = '-L(q)', type="b" )

#library(matlab)
#imagesc(as.matrix(sbmdata$MA))
#imagesc(finalr$Q, col=rev(gray.colors(5)), xlab ='group index', ylab = 'nodes')

##Gibbssampling method
library(CIDnetworks)
gibbsresultlist=list()
i=0
start <- Sys.time ()
for(k in Kvec) {
  i=i+1
  gibbsresultlist[[i]]=CID.Gibbs(as.matrix(0+sbmdata$A),components = list(SBM(k)))
}
Sys.time () - start
gibbslikelihood=numeric(length(Kvec))
for(i in 1:length(Kvec)) {
  gibbslikelihood[i]=gibbsresultlist[[i]]$CID.mean$log.likelihood
}
which(-gibbslikelihood==min(-gibbslikelihood))
