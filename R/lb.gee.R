# library(geepack)
lb.gee <- function(formula = formula(data),contrasts = NULL,subset,na.action,
                  data = parent.frame(),resporder=NULL,corstr,lfv=0.95,
                  id,control=lb.gee.control(),criter="cic",ext=TRUE,ext.lvl=0.999999,...)
{
  scall <- match.call()
  corstr<-substitute(corstr)
  if(!is.character(corstr)) corstr <- deparse(corstr)
  CORSTRS <- c("independence", "exchangeable", "ar1", "unstructured", "userdefined", "fixed")
  corstrv <- pmatch(tolower(corstr), CORSTRS, 0)
  criter<-substitute(criter)
  if(!is.character(criter)) criter <- deparse(criter)
  Criter<-c("RJ","QIC","CIC","SC","GP","RJ1","RJ2")
  criteri<-pmatch(toupper(criter),Criter,0)
  if(criteri==0) stop("\nOnly the following listed criterions are allowed:\nRJ, QIC, CIC, SC, GP.\n")

  if (corstrv == 0) stop("invalid corstr.")
  mnames <- c("", "formula", "data", "subset", "na.action", "id")
  cnames <- names(scall)
  cnames <- cnames[match(mnames,cnames,0L)]
  mcall <- scall[cnames]
  if(is.null(mcall$id)) mcall$id <- as.name("id")
  mcall[[1]] <- as.name("model.frame")
  # data <- model.frame(formula = formula, data = data, id=id)
  data<-eval(mcall, parent.frame())
  # dat.t<<-data
  na.action <- attr(data, "na.action")
  terms<-attr(data,"terms")
  # print(terms)
  y<-model.extract(data, "response")
  factor_index<-NULL
  if(is.factor(y)) {
    if(is.null(resporder))lvl<-levels(y) else lvl<-as.factor(resporder)
    if(length(lvl)!=2L) stop("Response variable has to be binary.")
    y_temp<-y
    y<-as.numeric(y)
    y[which(y_temp==lvl[1])]<-0
    y[which(y_temp==lvl[2])]<-1
    factor_index<-data.frame(Levels=lvl, Num=c(0,1))
  }
  id_old<-model.extract(data, id)
  if(is.null(id_old)) stop("id variable not found.")
  id_sort<-order(id_old)
  # print(id_sort)
  id<-cumsum(!duplicated(id_old[id_sort]))
  clusnew <- c(which(diff(as.numeric(id)) != 0), length(id))
  clusz <- c(clusnew[1], diff(clusnew))

  x<-model.matrix(terms,data=data,contrasts)
  contrasts<-attr(x,"contrasts")
  x<-x[id_sort,]
  y<-y[id_sort]
  # print(rformula)
  start<-lb.gee.init(x=x,y=y)
  # print(start)
  if(any(is.na(start$start)==TRUE))
    stop("Please check the model. Some covariates may be 100 percent explained by others.
        \rIf interaction term is included, please reform the data and generate them and
        \rrerun the function.")
  res.glm<-suppressWarnings(stats::glm.fit(x=x,y=y,
                family=binomial(link=log),
                start=start$start,control=glm.control()))
  # res.glm2<-suppressWarnings(stats::glm.fit(x=x,y=y,
  #               family=binomial(link=log),
  #               start=c(-1,rep(0,ncol(x)-1)),control=glm.control(maxit=100)))
  # print(res.glm1$coefficients)
  # print(res.glm2$coefficients)
  # start<-ifelse(res.glm1$deviance<res.glm2$deviance,
  #   res.glm1$coefficients,res.glm2$coefficients)
  # print(start)
  start<-res.glm$coefficients
  # str.t<<-res.glm$coefficients
  # x.t<<-x
  # y.t<<-y
  res.gee<-geepack::geese.fit(x,y,family=binomial(link=log),
                          b=start,corstr=corstr,
                          id=id,control=control,...)
  # if(is.list(res.gee)) print("aaa") else print("bbb")
  # print(res.gee$control$epsilon)
  # print(format(res.glm$deviance/2,nsmall=10))
  # print(format(res.qic$qic,nsmall=10))

  # print(res.gee$beta)
  mu<-as.vector(exp(x %*% res.gee$beta))
  # print(format(mu[which(mu>0.95)],nsmall=15))

  silly<-NULL
  # To active the triger to accept exact solution anyway.
  if(ext & any(mu>=ext.lvl)) silly<-TRUE
  # quasi<-sum(ifelse((1-y)==0,log(mu),log(1-mu)))
  start<-res.gee$beta
  criterion<-lb.gee.criter(object=res.gee,u=mu,
            criter=criteri,x=x,y=y,p=ncol(x),id=id,corstrv=corstrv)
  # print(criterion)
  x<-cbind(x,indicator=0, fitvalues=0)
  x[,"fitvalues"]<-mu
  x[,"indicator"][which(mu > lfv & y!=0)]<-1

  indicator<-x[,"indicator"]
  b_vectors <- as.data.frame(x)[which(x[,"indicator"]==1), ]
  # Remove the first column which is intercept
  b_vectors<-subset.data.frame(b_vectors,select=-1)
  x<-subset.matrix(x,select=-c(indicator,fitvalues))
  b_vectors<-b_vectors[order(b_vectors$fitvalues, decreasing = T),]
# print(b_vectors)
  colname<-colnames(b_vectors)[!is.na(match(colnames(b_vectors),colnames(x)))]
  b_vectors<-subset.data.frame(b_vectors,select=colname)

  # print(colnames(x))
  res<-c(res.gee,list(criterion=criterion,bv=0))
#
  # res.temp <-lb.gee.fit(x,y,indicator = indicator,
  #   bvectors = b_vectors, num_bvectors = nrow(b_vectors),
  #   start=start,id=id,corstrv=corstrv,criter=criteri,criterion=criterion,
  #   control=control,silly=silly)
# print(criterion)
  if(nrow(b_vectors) != 0) {
    res.temp <- tryCatch({
      suppressWarnings(lb.gee.fit(x,y,indicator = indicator,
        bvectors = b_vectors, num_bvectors = nrow(b_vectors),
        start=start,id=id,corstrv=corstrv,criter=criteri,criterion=criterion,
        control=control,silly=silly,...))
    },
      error=function(err) {
        return(ind<-0)
    })
    # if(is.list(res.temp))print("aaa") else print("bbb")
    if(is.list(res.temp))
      res <- c(res.temp, list(bv=1))
  }

  res <- c(res, list(factor=factor_index,formula = formula,criter=criteri,
          contrasts=contrasts,na.action=na.action,response=y))
  # print(res$bvectors)
  res$call<-scall
  res$clusz<-clusz
  # print(res$bvectors)
  if(!is.null(res$bvectors)) {
    res$bvectors<-data[row.names(res$bvectors),]
    # print(colnames(res$bvectors))
  }
  # print(res$qic$qic)
  class(res)<-"lb.gee"
  return(res)
}

lb.gee.fit<-function(x,y,id,bvectors = NULL,corstrv,indicator,control,criter,
                    criterion=NULL,num_bvectors = NULL,start=NULL,silly,...)
{
  p<-ncol(x)
  name_covariates<-colnames(x)
  if(is.null(num_bvectors))
    stop("No boundary vectors in the model.")
    # GP is looking for a positive maximum values.
    # So turn it to negative to be easy to encode.
  # print(criterion)
  # print(bvectors)
  # Block a silly results of qic when solution is on the boundary.
  if(is.null(silly)) silly<-FALSE
  if(criter==2 | criter==3) {
    if(abs(criterion[3])>500) silly<-TRUE
  }

  crit_old<-switch(criter,criterion,criterion[1],
                  criterion[3],criterion,-criterion,
                  criterion,criterion)
  # print((crit_old))
  # if(criter==2) crit_old<--criterion

  num_var<-ncol(x)-1

  names<-row.names(bvectors)
  bvectors_1e<-floor(bvectors*1e+5)
  ind_dup<-cumsum(!duplicated(bvectors_1e))
  bvectors<-cbind(bvectors,ind_dup=ind_dup)
  bvectors_temp <- bvectors
  # bvectors_temp <- bvectors[!duplicated(bvectors_1e[,1:(ncol(bvectors_1e)-1)]), ]
  # print(bvectors_1e)
  # print(duplicated(bvectors_1e[,1:(ncol(bvectors_1e)-1)]))
  max_num<-min(num_var,nrow(bvectors_temp))
  improved<-FALSE
  # print(max_num)

  for(j in 1:max_num)
  {
    cat("\nSearching for the admissible combinations of",j,"boundary vectors.\n")
    # print(bvectors_temp)
    # if(improved) print(bv_num)
    if(!improved) ind_matrix<-t(utils::combn(nrow(bvectors_temp),j))
      else ind_matrix<-t(utils::combn(nrow(bvectors_temp)-bv_num,j-bv_num))
      # print(ind_matrix)
    for(i in 1:nrow(ind_matrix))
    {
      # print(i)
      cat(".",sep="")
      if(i %% 10==0) cat("|",sep="")
      if(i %% 5==0 & i %% 10!=0) cat(" ",sep="")
      if(i %% 50==0) cat("\n",sep="")
      x_temp<-cbind(x, indicator=indicator)

      # print(colnames(x_temp))
      if(!improved) temp_bvectors<-bvectors_temp[ind_matrix[i,],]
        else temp_bvectors<-rbind(bvectors_temp[1:bv_num,],
          bvectors_temp[(bv_num+1):nrow(bvectors_temp),][ind_matrix[i,],])
      x_temp[,"indicator"][!(row.names(x_temp) %in%
          names[which(sapply(bvectors$ind_dup,
            function(t){t %in% temp_bvectors$ind_dup}))])]<-0
      indicator_temp<-x_temp[,"indicator"]
      x_temp<-subset(x_temp, select=-indicator)
      temp_bvectors<-subset(temp_bvectors,select=-ind_dup)
      # print(temp_bvectors)
      # print(colnames(x_temp))
      conv<-tryCatch({
        lb.gee.reformed(x = x_temp,y=y, t = temp_bvectors,
          num_bvectors = nrow(temp_bvectors),id=id,p=p,
          indicator = indicator_temp, start=start,
          control=control,corstrv=corstrv,criter=criter,...)
      },
        error=function(err){
          return(ind<-0)
      })
      if(is.list(conv)){
        # print(conv$criterion)
        if(criter==2 | criter==3) {
          if(abs(conv$criterion[3])>100) conv<-0
        }
      }

      # conv<-lb.gee.reformed(x = x_temp,y=y, t = temp_bvectors,
      #     num_bvectors = nrow(temp_bvectors),id=id,
      #     indicator = indicator_temp, start=start,
      #     control=control,corstrv=corstrv,criter=criter)

      # if(is.list(conv)) {
      #   if(conv$deviance<=deviance) break
      # }
      if(i==control$maxit) break

      # if(is.list(conv)) cat("\nBPS conv\n") else cat("\nBPS noconv\n")
      if(is.list(conv)) {

        # print(format(conv$qic$qic[3],nsmall=10))
          # 1 RJ, 2 QIC, 3 CIC, 4 SC, 5 GP.
        crit_new<-switch(criter,conv$criterion,conv$criterion[1],
                        conv$criterion[3],conv$criterion,-conv$criterion,
                        conv$criterion,conv$criterion)

        if(silly & !(improved)) {
          crit_old<-crit_new
          res<-conv
          improved<-TRUE
          bv_num<-j
          bvectors_new<-temp_bvectors

          if(nrow(temp_bvectors)!=nrow(bvectors)) {
            temp_fit<-conv$mu[row.names(bvectors_temp)[-which(row.names(bvectors_temp)%in%
                row.names(temp_bvectors))]]
            # Reorder the boundary vectors left in the bvectors_temp by their fitted probabilities.
            temp_fit_name<-names(temp_fit[order(temp_fit,decreasing = T)])
            temp_fit_name<-c(row.names(temp_bvectors),temp_fit_name)
            bvectors_temp<-bvectors_temp[temp_fit_name,]
            # print(temp_fit_name)
          }
          # if(corstrv==1) cat("\nMinimum QICu improved.\n")
          #   else cat(ifelse(criter==2,"\nMaximum ","\nMinimum "),
          #         switch(criter,"QIC","GP","SC")," improved.\n",sep = "")
          # # cat("\nCriterion improved.\n")
          #
          # if(i==nrow(ind_matrix)) cat(ifelse(criter==5,"\nMaximum ","\nMinimum "),
          #   switch(criter,"RJ","QIC","CIC","SC","GP"), " improved. ",
          #   switch(criter,"RJ","QIC","CIC","SC","GP"),"=",
          #   format(ifelse(criter==5,-crit_old,crit_old),nsmall=8),"\n",sep = "")
          break
        }
        # print(crit_old)
        # print(crit_new)
        if(crit_new<crit_old) {
          crit_old<-crit_new
          res<-conv
          improved<-TRUE
          bv_num<-j
          bvectors_new<-temp_bvectors
          if(nrow(temp_bvectors)!=nrow(bvectors)) {
            temp_fit<-conv$mu[row.names(bvectors_temp)[-which(row.names(bvectors_temp)%in%
                              row.names(temp_bvectors))]]
            # Reorder the boundary vectors left in the bvectors_temp by their fitted probabilities.
            temp_fit_name<-names(temp_fit[order(temp_fit,decreasing = T)])
            temp_fit_name<-c(row.names(temp_bvectors),temp_fit_name)
            bvectors_temp<-bvectors_temp[temp_fit_name,]
            # print(bvectors_temp)
          }
          # if(corstrv==1) cat("\nMinimum QICu improved.\n")
          #   else cat(ifelse(criter==2,"\nMaximum ","\nMinimum "),
          #         switch(criter,"QIC","GP","SC")," improved.\n",sep = "")
          # cat("\nCriterion improved.\n")
          break
        }
      }
    }
    if(j==max_num & !(improved))
      stop("No admissible pairs of boundary vectors could be found.")
  }
  if(improved) cat(ifelse(criter==5,"\nMaximum ","\nMinimum "),
    switch(criter,"RJ","QIC","CIC","SC","GP","RJ1","RJ2"), " improved. ",
    switch(criter,"RJ","QIC","CIC","SC","GP","RJ1","RJ2"),"=",
    format(ifelse(criter==5,-crit_old,crit_old),nsmall=8),"\n",sep = "")
  # print(res$beta)
  # res.t<<-res
  if(is.list(res)) {
    # print(res$beta)
    res.temp<-lb.gee.makeup(beta=res$beta,
                            vcov=res$vbeta,
                            t.repa=res$t.repa)
    # Reorganize the order
    res$beta<-res.temp$beta[name_covariates]
    res$vbeta<-res.temp$vcov[name_covariates,
                            name_covariates]
    vbeta.naiv<-lb.gee.makeup(vcov=res$vbeta.naiv,
                              t.repa=res$t.repa)
    res$vbeta.naiv<-vbeta.naiv[name_covariates,
                              name_covariates]
    # print(res$qic)
    # if(!is.null(vbeta.indep)) print(777) else print(888)
    if(isTRUE(control$jack)) {
      vbeta.ajs<-lb.gee.makeup(vcov=res$vbeta.ajs,
                              t.repa=res$t.repa)
      res$vbeta.ajs<-vbeta.ajs[name_covariates,
                              name_covariates]
    } else {
      num<-match("vbeta.ajs",names(res))
      res<-res[-c(num,num+1)]
    }
    if(isTRUE(control$j1s)) {
      vbeta.j1s<-lb.gee.makeup(vcov=res$vbeta.j1s,
                              t.repa=res$t.repa)
      res$vbeta.j1s<-vbeta.j1s[name_covariates,
                              name_covariates]
    } else {
      num<-match("vbeta.j1s",names(res))
      res<-res[-c(num,num+1)]
    }
    if(isTRUE(control$fij)) {
      vbeta.fij<-lb.gee.makeup(vcov=res$vbeta.fij,
                              t.repa=res$t.repa)
      res$vbeta.fij<-vbeta.fij[name_covariates,
                              name_covariates]
    } else {
      num<-match("vbeta.fij",names(res))
      res<-res[-c(num,num+1)]
    }
    # if(corstrv!=1 & (criter==1 | criter==4)) {
    #   vbeta.indep<-lb.gee.makeup(vcov=res$criterion$vbeta.indep,
    #                             t.repa=res$t.repa)
    #   vbeta.indep<-vbeta.indep[name_covariates,name_covariates]
    #   criterion<-lb.gee.qic(x=x,y=y,vbeta.indep=vbeta.indep,
    #                   beta=res$beta,Vr=res$vbeta)
    #   res$criterion<-criterion
    # }

    # if(is.list(qic)) print("ggg") else print("hhh")
    # res<-c(res, list(bvectors=bvectors_new))
    res$bvectors<-bvectors_new
    # print(res)
    return(res)
  } else {
    stop("No admissible combination of boundary vectors in model.")
  }
}

lb.gee.control <- function(epsilon = 1e-06, maxit = 100,
                          scale.fix = FALSE, jack = FALSE,
                          j1s = FALSE, fij = FALSE)
{
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of epsilon must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  list(trace = FALSE,jack=jack,j1s=j1s,fij=fij,
      maxit = as.integer(maxit), epsilon = epsilon)
}

lb.gee.makeup<-function(beta=NULL,vcov,t.repa)
{
  t<-cbind(0,t.repa)
  name.vcov<-c("(Intercept)",names(t.repa))
  if(!is.null(beta)) beta.complete<-rep(0,length(name.vcov))
  names(t)<-name.vcov
  if(!is.null(beta)) names(beta.complete)<-name.vcov

  vcov.complete<-matrix(0,nrow=length(name.vcov),
    ncol=length(name.vcov),
    dimnames=list(name.vcov,name.vcov))
  vcov.complete[row.names(vcov),colnames(vcov)]<-as.matrix(vcov)
  if(!is.null(beta)) beta.complete[names(beta)]<-beta
  # Makeup estimates and variance-covariance matrix.
  for (i in nrow(t.repa):1)
  {
    cov.vector<-apply(vcov.complete,1,function(x) {
      cov<--sum(x*t[i,])
      return(cov)
    })
    var.temp<-sum(diag(t[i,])%*%vcov.complete%*%t(t[i,]))
    vcov.complete[i,]<-vcov.complete[,i]<-cov.vector
    vcov.complete[i,i]<-var.temp
    if(!is.null(beta)) beta.complete[i]<--sum(t[i,]*beta.complete)
  }
  if(!(is.null(beta))) {
    return(list(beta=beta.complete,
              vcov=vcov.complete))
  } else return(vcov=vcov.complete)
}

summary.lb.gee <- function(object,CF.lvl=0.95, wald=FALSE,...)
{
  if(0.5>=CF.lvl || CF.lvl>=1)
    stop("CF.lvl should be a value between 0.5 and 1.")
  alpha<-1-CF.lvl

  mean.sum <- data.frame(estimate = object$beta,
                        san.se = sqrt(diag(object$vbeta)))
  mean.sum$stat <- mean.sum$estimate / mean.sum$san.se
  mean.sum$p <- 2*pnorm(-abs(mean.sum$stat))
  confint<-apply(mean.sum[,1:2], 1, function(t){
    lower<-t[1]+qnorm(alpha/2)*t[2]
    upper<-t[1]+qnorm(1-alpha/2)*t[2]
    dat<-cbind(lower,upper)
  })
  if(!wald) mean.sum<-cbind(mean.sum,t(confint))
  ## Borrow idea from the Stata to present confidence interval.
  confint_name<-c(paste("[",(1-alpha)*100,"% ","Conf.",sep=""),"Interval]")
  dn<-c("Estimate", "San.StdErr","z value","Pr(>|z|)")
  dn_mean.sum<-c(dn,confint_name)
  corr.sum <- data.frame(estimate = object$alpha,
                        san.se = sqrt(diag(object$valpha)))
  corr.sum$stat <- corr.sum$estimate / corr.sum$san.se
  corr.sum$p <- 2*pnorm(-abs(corr.sum$stat))

  scale.sum <- data.frame(estimate = object$gamma,
    san.se = sqrt(diag(object$vgamma)))
  scale.sum$stat <- scale.sum$estimate / scale.sum$san.se
  scale.sum$p <- 2*pnorm(-abs(scale.sum$stat))
  if(wald) {

    dn_mean.sum<-
      dn<-c("Estimate", "San.StdErr","Wald","P-value")
    mean.sum$stat <-mean.sum$stat^2
    mean.sum$p <- 1 - pchisq(mean.sum$stat, df=1)
    corr.sum$stat<-corr.sum$stat^2
    corr.sum$p <- 1 - pchisq(corr.sum$stat, df=1)
    scale.sum$stat<-scale.sum$stat^2
    scale.sum$p <- 1 - pchisq(scale.sum$stat, df=1)
  }

  if(object$control$jack) {
    cbind(mean.sum,ajs.se = sqrt(diag(object$vbeta.ajs)))
    cbind(corr.sum,ajs.se = sqrt(diag(object$valpha.ajs)))
    cbind(scale.sum,ajs.se = sqrt(diag(object$vgamma.ajs)))
  }
  if(object$control$j1s) {
    cbind(mean.sum,j1s.se = sqrt(diag(object$vbeta.j1s)))
    cbind(corr.sum,j1s.se = sqrt(diag(object$valpha.j1s)))
    cbind(scale.sum,j1s.se = sqrt(diag(object$vgamma.j1s)))
  }
  if(object$control$fij) {
    cbind(mean.sum,fij.se = sqrt(diag(object$vbeta.fij)))
    cbind(corr.sum,fij.se = sqrt(diag(object$valpha.fij)))
    cbind(scale.sum,fij.se = sqrt(diag(object$vgamma.fij)))
  }
  if(wald) {
    dimnames(mean.sum) <- list(names(object$beta),dn)
  } else dimnames(mean.sum) <- list(names(object$beta),dn_mean.sum)
  if(nrow(corr.sum) > 0)
    dimnames(corr.sum) <- list(names(object$alpha),dn)
  if(nrow(scale.sum) > 0)
    dimnames(scale.sum) <- list("scale",dn)

  ans <- list(mean=mean.sum, correlation=corr.sum, model=object$model,
              factor=object$factor,call=object$call, clusz=object$clusz,
              control=object$control,error=object$err,scale=scale.sum,
              bvector=object$bvector,criter=object$criter,criterion=object$criterion)
  class(ans) <- "summary.lb.gee"
  ans
}

print.summary.lb.gee <- function(x, digits = max(5L, getOption("digits")-1L),
                                quote = FALSE, prefix = "",... )
{
  cat("\nCall:\n")
  dput(x$call)
  cat("\nMean Model:\n")
  cat(" Mean Link:                ", x$model$mean.link, "\n")
  cat(" Variance to Mean Relation:", x$model$variance, "\n")
  cat("\nBoundary Vectors:\n")
  print(x$bvector, digits=digits-2L)
  if(!is.null(x$factor)) {
    cat("\nIndex of response variable:\n")
    factor_names<-as.character(x$factor[,1])
    factor_names<-format(factor_names,
      width=max(nchar(factor_names)),
      justify="right")
    cat(factor_names[1],"$ ",x$factor[,2][1],"\n",
      factor_names[2],"$ ",x$factor[,2][2],"\n",sep="")
  }
  cat("\nCoefficients:\n")
  print(x$mean, digits = digits)

  cat("\nScale Model:\n")
  cat(" Scale Link:    ", x$model$sca.link, "\n")
  cat("\nEstimated Scale Parameters:\n")
  print(x$scale, digits = digits)

  cat("\nCorrelation Model:\n")
  cat(" Correlation Structure:    ", x$model$corstr, "\n")
  if (pmatch(x$model$corstr, c("independence", "fixed"), 0) == 0) {
    cat(" Correlation Link:         ", x$model$cor.link, "\n")
    cat("\n Estimated Correlation Parameters:\n")
    print(x$corr, digits = digits)
  }
  criter<-switch(x$criter,"RJ","QIC","CIC","SC","GP","RJ1","RJ2")
  cat("\n",criter,": ",switch(criter, x$criterion,
      QIC=x$criterion[1],CIC=x$criterion[3]),"\n",sep="")
  ##cat("\nNumber of observations : ", x$nobs, "\n")
  ##cat("\nMaximum cluster size   : ", x$max.id, "\n")

  cat("\nReturned Error Value:    ")
  cat(x$error, "\n")
  cat("Number of clusters:  ", length(x$clusz), "  Maximum cluster size:", max(x$clusz), "\n\n")
  invisible(x)
}

print.lb.gee <- function(x, digits = max(3L, getOption("digits")-1L),
                        quote = FALSE, prefix = "", ...)
{
  cat("\nCall:\n")
  dput(x$call)
  cat("\nMean Model:\n")
  cat(" Mean Link:                ", x$model$mean.link, "\n")
  cat(" Variance to Mean Relation:", x$model$variance, "\n")
  cat("\nBoundary Vectors:\n")
  print(x$bvector, digits=digits-2L)
  if(!is.null(x$factor)) {
    cat("\nIndex of response variable:\n")
    factor_names<-as.character(x$factor[,1])
    factor_names<-format(factor_names,
      width=max(nchar(factor_names)),
      justify="right")
    cat(factor_names[1],"$ ",x$factor[,2][1],"\n",
      factor_names[2],"$ ",x$factor[,2][2],"\n",sep="")
  }
  cat("\n Coefficients:\n")
  print(unclass(x$beta), digits = digits)

  cat("\nScale Model:\n")
  cat(" Scale Link:               ", x$model$sca.link, "\n")
  cat("\nCorrelation Model:\n")
  cat(" Correlation Structure:    ", x$model$corstr, "\n")
  if (pmatch(x$model$corstr, c("independence", "fixed"), 0) == 0) {
    cat(" Correlation Link:         ", x$model$cor.link, "\n")
    cat("\n Estimated Correlation Parameters:\n")
    print(unclass(x$alpha), digits = digits)
  }
  criter<-switch(x$criter,"RJ","QIC","CIC","SC","GP","RJ1","RJ2")
  cat("\n",criter,": ",switch(criter, x$criterion,
      QIC=x$criterion[1],CIC=x$criterion[3]),"\n",sep="")

  ##cat("\nNumber of observations : ", x$nobs, "\n")
  ##cat("\nMaximum cluster size   : ", x$max.id, "\n")

  cat("\nReturned Error Value:  ")
  cat(x$error, "\n")
  cat("Number of clusters:  ", length(x$clusz), "  Maximum cluster size:", max(x$clusz), "\n\n")
  invisible(x)
}
