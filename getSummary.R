############ getSummary ##############################################
###### Arguments
### model
# The stanfit fit using data returned by getStanDataInfo and
# the stan code in Kimball, Shantz, Eager, and Roy (2016)
#
### datainfo
# The list returned by getStanDataInfo
#
###### Value: a list with the following elements
### 'formula'
# The formula passed to getStanDataInfo
#
### 'response'
# A character vector with the name of the response variable
# and its levels
#
### 'fixed'
# A list with the fixed regression estimates which includes the
# posterior probability that each estimate is greater than zero,
# plus the fixed effect covariance and correlation matrices
#
### 'random'
# A list with a sub-list for subjects and a sub-list for items.
# Each list includes the estimates for individual subjs/items
# (ranef), the betas summed with these estimates (coef),
# the standard deviation in each effect (sd), and the effect
# covariance and correlation matrices (cov and cor)
#
### 'pars'
# A list of matrices, each being the summary for a given set
# parameters from stan's summary() function, but renamed to
# match the data provided to getStanDataInfo.
#
### 'parnames'
# A list that matches the names in the stanfit to the names
# in the summary, along with dimnames which can be assigned
# to extracted parameters

getSummary <- function(model,datainfo){
    d <- datainfo$data
    i <- datainfo$info

    # individual parameter summaries
    pn <- list(
        beta = list(stanname="beta",
            dimnames=list(NULL,i$P)),
        gamma_subj = list(stanname="gamma_subj",
            dimnames=list(NULL,i$S,i$QS)),
        sigma_subj = list(stanname="sigma_subj",
            dimnames=list(NULL,i$QS)),
        omega_subj = list(stanname="omega_subj",
            dimnames=list(NULL,i$QS,i$QS)),
        gamma_item = list(stanname="gamma_item",
            dimnames=list(NULL,i$I,i$QI)),
        sigma_item = list(stanname="sigma_item",
            dimnames=list(NULL,i$QI)),
        omega_item = list(stanname="omega_item",
            dimnames=list(NULL,i$QI,i$QI)),
        predicted = list(stanname="y_hat",
            dimnames=list(NULL,paste(1:d$N))))
    names(pn)[2:7] <- c(paste("gamma",i$subj,sep="_"),
        paste("sigma",i$subj,sep="_"),paste("omega",i$subj,sep="_"),
        paste("gamma",i$item,sep="_"),paste("sigma",i$item,sep="_"),
        paste("omega",i$item,sep="_"))

    ps <- list()
    for(p in names(pn)){
        sn <- pn[[p]]$stanname
        dn <- pn[[p]]$dimnames
        rn <- expand.grid(dn[length(dn):2])
        if(ncol(rn)>1) rn <- rn[,ncol(rn):1]
        rn <- apply(rn,1,function(x) paste(x,collapse=","))

        ps[[p]] <- summary(model,probs=c(.025,.975),pars=sn)$summary
        pn[[p]]$rownames <- data.frame(stan = rownames(ps[[p]]),
            named = paste(p,"[",rn,"]",sep=""), stringsAsFactors=F)
        rownames(ps[[p]]) <- pn[[p]]$rownames$named
    }

    # remove upper tri and diag from correlation matrices
    for(o in c(4,7)){
        keep <- as.vector(lower.tri(diag(sqrt(nrow(ps[[o]])))))
        ps[[o]] <- ps[[o]][keep,]
        if(!is.matrix(ps[[o]])){  # if there is only one slope
            ps[[o]] <- matrix(ps[[o]],1,7)
            rownames(ps[[o]]) <- pn[[o]]$rownames$named[3]
        }
    }

    # fixed
    b <- extract(model,pars="beta")$beta
    ppos <- apply(b,2,function(x) mean(x>0))
    bcov <- var(b)
    bcor <- cor(b)
    dimnames(bcov) <- dimnames(bcor) <- list(i$P,i$P)
    b <- cbind(i$P,data.frame(ps$beta[,c(1,3:5)]))
    rownames(b) <- NULL
    colnames(b) <- c("Name","Estimate","SD","2.5%","97.5%")
    b[,"P(B>0)"] <- ppos
    f <- list(coef = b, cov = bcov, cor = bcor)

    # random
    subj <- list()
    subj$ranef <- matrix(ps[[2]][,1],d$S,d$QS,byrow=T,
        dimnames=list(i$S,i$QS))
    subj$coef <- matrix(ps$beta[,1],d$S,d$P,byrow=T,
        dimnames=list(i$S,i$P))
    subj$coef[,i$QS] <- subj$coef[,i$QS] + subj$ranef
    subj$sd <- ps[[3]][,1]
    names(subj$sd) <- i$QS
    subj$cor <- matrix(1,d$QS,d$QS,dimnames=list(i$QS,i$QS))
    subj$cor[lower.tri(subj$cor)] <- ps[[4]][,1]
    subj$cor[upper.tri(subj$cor)] <- ps[[4]][,1]
    subj$cov <- diag(subj$sd) %*% subj$cor %*% t(diag(subj$sd))
    dimnames(subj$cov) <- dimnames(subj$cor)

    item <- list()
    item$ranef <- matrix(ps[[5]][,1],d$I,d$QI,byrow=T,
        dimnames=list(i$I,i$QI))
    item$coef <- matrix(ps$beta[,1],d$I,d$P,byrow=T,
        dimnames=list(i$I,i$P))
    item$coef[,i$QI] <- item$coef[,i$QI] + item$ranef
    item$sd <- ps[[6]][,1]
    names(item$sd) <- i$QI
    item$cor <- matrix(1,d$QI,d$QI,dimnames=list(i$QI,i$QI))
    item$cor[lower.tri(item$cor)] <- ps[[7]][,1]
    item$cor[upper.tri(item$cor)] <- ps[[7]][,1]
    item$cov <- diag(item$sd) %*% item$cor %*% t(diag(item$sd))
    dimnames(item$cov) <- dimnames(item$cor)

    r <- list(subj = subj, item = item)
    names(r) <- c(i$subj,i$item)

    # some diagnostics
    ndiv <- sum(do.call(rbind, args = get_sampler_params(
        model, inc_warmup = FALSE))[,"divergent__"])
    allpars <- do.call(rbind, args = ps)
    rhat <- sum(allpars[,"Rhat"]>1.1)
    neff <- sum(allpars[,"n_eff"]<100)
    d$U <- 2*d$P+1 + (d$S+1)*d$QS+choose(d$QS,2) + 
		(d$I+1)*d$QI+choose(d$QI,2)

    mess <- character()
    if(ndiv>0) mess <- c(mess,paste("WARNING:",ndiv,
        "divergent transitions post-warmup"))
    if(rhat>0) mess <- c(mess,paste("WARNING:",rhat,
        "parameters have Rhat values greater than 1.1"))
    if(neff>0) mess <- c(mess,paste("WARNING:",neff,
        "parameters have effective sample sizes under 100"))
    if(d$U>d$N) mess <- c(mess,paste("WARNING: model has more",
        "unconstrained parameters than observations"))

    return(list(formula = i$formula, response = i$y, fixed = f,
        random = r, pars = ps, parnames = pn,
        dims = d[c("N","S","I","P","QS","QI","U")], messages = mess))
}

######################################################################
