lb.gee.reformed <- function(x,y, t = NULL, id=NULL,corstrv,p=NULL,
                            num_bvectors = NULL, indicator = NULL,
                            start=NULL, control,criter,...)
{
  sub_x<-subset.matrix(x,select=-1)
  rownames.sub_x<-row.names(sub_x)
  # t <- as.data.frame(bvectors)
  rownames.bv<-row.names(t)
  colnames.bv<-colnames(t)
  t.repa<-as.data.frame(matrix(NA,ncol=ncol(t),
                        nrow=nrow(t),byrow=T,
                        dimnames=dimnames(t)))
  for(r in 1:num_bvectors)
  {
    if(r == 1) {
      sub_x<-matrix(unlist(apply(sub_x,1,function(z) z-t[r,])),ncol=ncol(sub_x),byrow=T)
      dimnames(sub_x)<-list(rownames.sub_x,colnames.bv)
      t.repa[r,]<-temp<-t[r, ]
      t<-matrix(unlist(apply(t,1,function(z) z-temp)),ncol=ncol(t),byrow=T)
      dimnames(t)<-list(rownames.bv,colnames.bv)
    } else {
      q <- r-1
      if(t[r,q]!=0) {
        temp <- t[r, ]/t[r, q]
      } else {
        # Find nonzero member in current vector
        nonzero<-which(!t[r,1:ncol(t)]%in%0)
        # print(nonzero)
        nonzero.col<-nonzero[1]
        # Switch the position of current column with the nonzero column
        # print(colnames.bv)
        name.temp<-colnames.bv[q]
        colnames.bv[q]<-colnames.bv[nonzero.col]
        colnames.bv[nonzero.col]<-name.temp
        # print(colnames.bv)
        sub_x<-sub_x[,colnames.bv]
        t<-t[,colnames.bv]

        t.repa<-t.repa[,colnames.bv]
        temp <- t[r, ]/t[r, q]
      }
      t.repa[r,]<-temp
      sub_x<-matrix(unlist(apply(sub_x,1,function(z) z-temp*z[q])),ncol=ncol(sub_x),byrow=T)
      dimnames(sub_x)<-list(rownames.sub_x,colnames.bv)
      t<-matrix(unlist(apply(t,1,function(z) z-temp*z[q])),ncol=ncol(t),byrow=T)
      dimnames(t)<-list(rownames.bv,colnames.bv)
      # print(sub_x)
    }
  }
  temp_x <- cbind(indicator=indicator, sub_x)
  temp_x <- temp_x[-which(temp_x[,"indicator"]==1), ]
  y<-y[-which(indicator==1)]
  id<-id[-which(indicator==1)]
  sub_x<-subset(temp_x,select=-c(1:num_bvectors))

  # print(formula)
  # Fit the modle with reparameterized dataset.
  # Adjusting the starting values to fit the model without constant.
  # init<-lb.gee.init(x=sub_x,y=y,noconstant=T)
  # # print(init)
  # # if(is.list(init))print(111) else print(222)
  # if(init$SAS==0) {
  #   res.glm<-suppressWarnings(stats::glm.fit(x=sub_x,y=y,
  #                 family = binomial(link = log),
  #                 mustart=init$mustart))
  # } else {
  #   res.glm<-suppressWarnings(stats::glm.fit(x=sub_x,y=y,
  #                 family = binomial(link = log),
  #                 start=init$start))
  # }
  #

  # print(format(start,nsmall=15))
  # if(is.list(res.glm))print(333) else print(444)
  res.glm<-suppressWarnings(stats::glm.fit(x=sub_x,y=y,
                    family = binomial(link = log),
                    start=start[colnames(sub_x)]))
  start <- res.glm$coefficients
  corstr<-switch(corstrv,"independence", "exchangeable", "ar1",
                "unstructured", "userdefined", "fixed")
  results_temp <- suppressWarnings(geepack::geese.fit(x=sub_x,y=y,
                                  family = binomial(link = log),
                                  id=id,b=start,corstr=corstr,
                                  control=control,...))
  # if(is.list(results_temp))print(111) else print(222)
  # print(results_temp$alpha)
  # print(results_temp$gamma)
  mu<-exp(sub_x %*% results_temp$beta)
  names(mu)<-row.names(sub_x)
  # print(mu)
  # quasi<-sum(ifelse((1-y)==0,log(mu),log(1-mu)))
  # print(quasi)
  # start<-results_temp$beta
  # criterion<-NULL
  # print(corstrv)
  # print(criter)
  criterion<-lb.gee.criter(object=results_temp,u=mu,
                criter=criter,x=sub_x,y=y,p=p,id=id,corstrv=corstrv)
  # if(is.null(criterion)) print(555) else print(666)
  res.members<-c("beta", "alpha","vbeta",
                "valpha", "valpha.stab",
                "gamma","vgamma",
                "vbeta.naiv", "valpha.naiv",
                "vbeta.ajs", "valpha.ajs","vgamma.ajs",
                "vbeta.j1s", "valpha.j1s","vgamma.j1s",
                "vbeta.fij", "valpha.fij","vgamma.fij",
                "error", "model", "control")
  results_temp<-results_temp[match(res.members,names(results_temp))]
  dimnames(results_temp$vbeta)<-
    dimnames(results_temp$vbeta.naiv)<-
    dimnames(results_temp$vbeta.ajs)<-
    dimnames(results_temp$vbeta.j1s)<-
    dimnames(results_temp$vbeta.fij)<-
      list(colnames(sub_x),colnames(sub_x))
    # names(results_temp$beta)<-name[2:length(name)]

  res <- c(results_temp,list(criterion=criterion,t.repa=t.repa,mu=mu))
  # if(is.list(res)) print("1") else print("0")
  return(res)
}
