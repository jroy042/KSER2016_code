############ getStanDataInfo #########################################
###### Arguments
### 'formula'
# A formula as would be passed to glmer()
# e.g. y ~ x1 + x2 + (1 + x1 + x2 | subject) + (1 + x1 + x2 | item)
# All numeric variables included in the formula should be scaled
# to mean 0 and standard deviation 1, and logicals
# (i.e. 0/1 predictors or T/F predictors) should be converted to
# factors prior to running the function.  The function automatically
# sets sum contrasts for unordered factors and polynomial contrasts
# for ordered factors.
#
### 'data'
# The data.frame which contains the data to apply the formula to.
#
### 'subj'
# Character string indicating which of the two groups in the formula
# above corresponds to subjects, e.g. "subject", "speaker", etc.
# This is required so the names of random effects parameters in the
# output match correctly.
#
### 'item'
# Character string indicating which of the two groups in the formula
# above corresponds to items, e.g. "item", "word", etc.
# This is required so the names of random effects parameters in the
# output match correctly.
#
###### Value: a list with two elements
### 'data'
# A list which can be passed as the 'data' argument in stan() using
# the stan code in Kimball, Shantz, Eager, and Roy (2016).
#
### 'info'
# A list which contains information about the names of the parameters
# that will be returned by stan().  It also contains a character
# vector 'keep' with the names of the transformed parameters which
# can be passed to the 'pars' argument in stan() to prevent the raw
# parameters from being returned as part of the model (saves space)

getStanDataInfo <- function(formula,data,subj,item){
    require(lme4)
    require(rstan)
    require(stringr)
    options(contrasts=c("contr.sum","contr.poly"))

    # get the response and fixed effects model matrix
    fixd <- nobars(formula)
    fr <- model.frame(fixd,data,drop.unused.levels=T)
    y <- factor(model.response(fr))
    ynames <- c(colnames(fr)[1],levels(y))
    y <- as.numeric(y) - 1
    x <- model.matrix(fixd,fr)
    N <- nrow(x)
    P <- length(Pnames <- colnames(x))

    # get the random effects matrices
    rand <- mkReTrms(findbars(formula),data)
    S <- length(Snames <- levels(rand$flist[,subj]))
    QS <- length(QSnames <- rand$cnms[[subj]])
    I <- length(Inames <- levels(rand$flist[,item]))
    QI <- length(QInames <- rand$cnms[[item]])

    # in case interaction order is different between P/Q
    Psplit <- lapply(str_split(Pnames,":"),sort)
    fixs <- which(!(QSnames %in% Pnames))
    fixi <- which(!(QInames %in% Pnames))
    for(i in fixs){
        j <- sort(str_split(QSnames[i],":")[[1]])
        p <- which(vapply(Psplit,function(x){
            ifelse(length(x)==length(j),
            all(x==j),FALSE)},logical(1)))
        if(length(p)!=1) stop("random effect not in fixed effects")
        QSnames[i] <- Pnames[p]
    }
    for(i in fixi){
        j <- sort(str_split(QInames[i],":")[[1]])
        p <- which(vapply(Psplit,function(x){
            ifelse(length(x)==length(j),
            all(x==j),FALSE)},logical(1)))
        if(length(p)!=1) stop("random effect not in fixed effects")
        QInames[i] <- Pnames[p]
    }

    OSnames <- expand.grid(QSnames,QSnames)[,2:1]
    OInames <- expand.grid(QInames,QInames)[,2:1]
    GSnames <- expand.grid(QSnames,Snames)[,2:1]
    GInames <- expand.grid(QInames,Inames)[,2:1]

    # combine fixed and random matrices and convert to CSR
    x <- extract_sparse_parts(cbind(x,
        t(rand$Ztlist[[which(names(rand$cnms)==subj)]]),
        t(rand$Ztlist[[which(names(rand$cnms)==item)]])))
    nz <- length(x$w)

    # return data and info
    return(list(
        data = list(N = N, S = S, I = I, P = P, QS = QS, QI = QI,
            y =y, nz = nz, x_w = x$w, x_v = x$v, x_u = x$u),
        info = list(formula = formula, subj = subj, item = item,
            y = ynames, P = Pnames, S = Snames, QS = QSnames,
            I = Inames, QI = QInames, OS = OSnames, OI = OInames,
            GS = GSnames, GI = GInames, keep = c("beta","y_hat",
            "gamma_subj","sigma_subj","omega_subj","gamma_item",
            "sigma_item","omega_item"))))
}

######################################################################
