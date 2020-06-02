lb.gee.criter<-function(object,u,criter,x,y,p,id,corstrv) {
  # RJ
  if(criter==1 | criter==6 | criter==7) {
    # Robust VC
    Vr<-object$vbeta
    # Naive VC
    V<-object$vbeta.naiv
    gamma<-Vr %*% MASS::ginv(V)
    rj1<-sum(diag(gamma))/p
    rj2<-sum(diag(gamma)^2)/p
    # cat("\n",rj1,rj2,"\n")
    if(criter==1) res<-sqrt((1-rj1)^2+(1-rj2)^2)
    else if(criter==6) res<-abs(1-rj1)
    else res<-abs(1-rj2)
  }
  # QIC and CIC
  if(criter==2 | criter==3) {
    res<-lb.gee.QIC(model.R=object, x=x,y=y, mu=u,id=id)
  }
  # SC
  if(criter==4 | criter==5) {
    res<-lb.gee.ess(y,p,id,u,object$gamma,object$alpha,corstrv,criter)
  }
  return(res)
}
# library(MASS)
lb.gee.ess <- function(y,p=NULL,id,u,scale,
                    alpha,corstrv,criter) {
  lvl<-unique(id)
  num<-length(lvl)
  # print(num)
  res_t<-sapply(1:num, function(t) {
    u_t<-round(u[which(id==lvl[t])],digits=15)
    y_t<-y[which(id==lvl[t])]
    # x_t<-as.matrix(x[which(id==lvl[t]),])
    if(length(u_t)!=1L) {
      # print(length(u_t))
      corr_mat<-corStr(waves=length(u_t),corstrv=corstrv,alpha=alpha)
      v_t<-diag(u_t*(1-u_t))
      V_t<-sqrt(v_t)%*%corr_mat%*%sqrt(v_t)*scale
      ess_t<-t(y_t-u_t)%*%MASS::ginv(V_t)%*%(y_t-u_t)
      # ess_t<-t(y_t-u_t)%*%solve(V_t)%*%(y_t-u_t)
      if(criter==5) {
        gp_t<--0.5*(ess_t+log(det(V_t)))
      }
    } else {
      v<-u_t*(1-u_t)*scale
      ess_t<-(y_t-u_t)*(1/v)%*%t(y_t-u_t)
      if(criter==5) {
        gp_t<--0.5*(ess_t+log(abs(v)))
      }
    }
    if(criter==5) return(gp_t)
      else return(ess_t)
  })
  # print(res_t)
  # print(sum(res_t))
  # GP
  if(criter==5) res<-Reduce("+",res_t)
    # SC
    else res<-Reduce("+",res_t)/(length(y)-p-length(alpha))
  return(res)
}

corStr<-function (waves, corstrv,alpha)
{
  if (corstrv == 1) {
    # In
    cormat<-diag(waves)
  } else if (corstrv == 2) {
    # Ex
    cormat<-(diag(waves)+1)*alpha
    # print(cormat)
    diag(cormat)<-1
  } else if (corstrv == 3) {
    # AR-1
    cormat<-diag(waves)
    mat_temp<-sapply(1:(waves-1),function(m) sapply(which(upper.tri(cormat)[m,])-m, function(n) alpha^n))
    for(i in 1:(waves-1)) {
      cormat[i,which(upper.tri(cormat)[i,])]<-as.vector(mat_temp[[i]])
    }
    cormat[lower.tri(cormat)]<-t(cormat)[lower.tri(cormat)]
  } else {
    if(!length(alpha)>1) stop("\nThe length of alpha is wrong.\n")
    col_names<-sapply(1:(waves-1),function(m) crossutri(waves,m))
    znames <- lapply(col_names, function(m) paste("alpha", m, sep = "."))
    # Unstr
    cormat<-diag(waves)
    for(i in 1:(waves-1)) {
      cormat[i,which(upper.tri(cormat)[i,])]<-as.vector(alpha[znames[[i]]])
    }
    cormat[lower.tri(cormat)]<-t(cormat)[lower.tri(cormat)]
  }
  return(cormat)
}

crossutri <- function(waves,n) {
  # ans<-rep(0,waves-n)
  # print(ans)
  for (i in (n+1):waves) {
    if(i==n+1) ans<-paste(n, i, sep=":")
      else ans <- c(ans, paste(n, i, sep=":"))
  }
  ans
}

lb.gee.QIC <- function(model.R=NULL, x,y, mu=NULL,id)
{
  # beta<-model.R$beta
  # Robust VC
  Vr<-model.R$vbeta
  # print(beta)
  # print(Vr)
  # print(vbeta.indep)
  # if(is.list(model.R)) print("uuu") else print("vvv")
  # Quasilikelihood
  # names(mu.R)<-row.names(x)
  # print(format(mu[which(mu>0.9)],nsmall=10))
  scale=model.R$gamma
  quasi.R<-sum(ifelse((1-y)==0,log(mu),log(1-mu)))

  # # Trace Term (penalty for model complexity)
  # x.t<<-x
  # id.t<<-id
  # mu.t<<-mu
  # scale.t<<-scale
  omegaI<-lb.gee.omegaI(x,id,mu,scale)
  trace.R<-sum(diag(omegaI %*% Vr)) # CIC

  # QIC
  QIC<-2*(trace.R - quasi.R)
  # QICu<-2*(p - quasi.R)    # Approximation assuming model structured correctly
  # qic<-c(QIC, QICu, quasi.R, trace.R, p)
  qic<-c(QIC, quasi.R, trace.R)

  # print(qic)
  # names(qic) <- c('QIC','QICu', 'Quasi-Lik', 'Trace', 'px')
  # print(vbeta.indep)
  names(qic) <- c('QIC', 'Quasi-Lik', 'Trace')
  output<-qic
  return(output)
}

lb.gee.omegaI<-function(x,id,u,scale)
{
  lvl<-unique(id)
  num<-length(lvl)
  omega_t<-lapply(c(1:num), function(t){
    u_t<-u[which(id==lvl[t])]
    if(length(u_t)!=1L) {
      # inv_v<-diag(ifelse(u_t>=0.99999,1/u_t,1/(u_t*(1-u_t))))
      inv_v<-diag(1/(scale*u_t*(1-u_t)))
      D<-as.matrix(x[which(id==lvl[t]),])*as.vector(u_t)
      omega_temp<-t(D)%*%as.matrix(inv_v)%*%as.matrix(D)
    } else {
      inv_v<-1/(u_t*(1-u_t))
      D<-as.matrix(x[which(id==lvl[t]),])*as.vector(u_t)
      omega_temp<-(D*inv_v)%*%t(D)
    }
    # cat(t,"\n")
  })
  omega<-Reduce("+",omega_t)
  return(omega)
}


#
# lb.gee.naive<-function(X,clusz,u,scale,alpha,corstrv)
# {
#   lvl<-unique(id)
#   num<-length(lvl)
#   v_t<-lapply(clusz, function(t){
#     u_t<-u[which(id==lvl[t])]
#     if(length(u_t)!=1L) {
#       corr_mat<-corStr(waves=clusz[t],corstr=corstr,alpha=alpha)
#       v<-diag(u_t*(1-u_t))
#       V<-sqrt(v)%*%corr_mat%*%sqrt(v)*scale
#     } else {
#       v<-u_t*(1-u_t)
#     }
#     x<-as.matrix(X[which(id==lvl[t]),])*as.vector(u_t)
#     H<-t(x)%*%solve(V)%*%as.matrix(x)
#   })
#   v<-Reduce("+",v_t)
#   return(v)
# }
