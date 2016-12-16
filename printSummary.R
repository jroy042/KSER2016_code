############ printSummary ############################################
###### Argument: list returned getSummary
###### Value: prints the major parts of the summary

printSummary <- function(summ){
    w <- getOption("width")
    options(width = 200)

    cat("Mixed effects logistic regression fit with Stan\n\nFormula:")
    show(summ$formula)
    cat("\n")

    cat(paste("Obervations:",summ$dims$N),
        paste("Unconstrained Parameters:",summ$dims$U),sep="\n")

    y <- summ$response
    if(any(y[-1]!=c("0","1"))){
        y <- paste(y[1]," (",paste(paste(y[-1],"->",c("0","1")),
            collapse="; "),")",sep="")
    } else y <- y[1]
    cat(paste("Response:",y),"\n\n")

    cat("Random Effects:\n\n")
    for(r in 1:2){
        G <- ifelse(r==1,summ$dims$S,summ$dims$I)
        cat(paste(names(summ$random)[r]," (",G,")",sep=""),"\n")
        rsd <- format(round(summ$random[[r]]$sd,3),nsmall=3)
        rsc <- round(data.frame(summ$random[[r]]$cor),2)
        firstneg <- any(rsc[,1]<0)
        rsc <- format(rsc,nsmall=2)
        rsc[upper.tri(rsc,diag=T)] <- ""
        rsc <- cbind(names(rsd),rsd,rsc)
        rsc <- rsc[,-ncol(rsc)]
        colnames(rsc)[1:2] <- c("Name","SD")
        colnames(rsc)[3] <- ifelse(firstneg," Corr","Corr")
        if(ncol(rsc)>3) colnames(rsc)[4:ncol(rsc)] <- ""
        print(rsc,right=F,row.names=F)
        cat("\n")
    }

    cat("Fixed Effects:\n")
    b <- summ$fixed$coef
    b[,2:5] <- format(round(b[,2:5],3),nsmall=3)
    p <- b[,6] <- round(b[,6],4)
    b[,6] <- substr(format(b[,6],nsmall=4),2,6)
    for(j in 2:5){
        if(any(round(summ$fixed$coef[,j],3)<0)){
            colnames(b)[j] <- paste(" ",colnames(b)[j],sep="")
        }
    }
    b[p==0,6] <- "<.0001"
    b[p==1,6] <- ">.9999"
    b[p!=0&p!=1,6] <- paste(" ",b[p!=0&p!=1,6],sep="")
    colnames(b)[6] <- " P(B>0)"
    print(b,row.names=F,right=F)
	
	allpars <- do.call(rbind,args=summ$pars)
	cat("","Effect Sample Sizes:",sep="\n")
	show(summary(allpars[,"n_eff"]))
	cat("","Rhat:",sep="\n")
	show(summary(allpars[,"Rhat"]))

    if(length(summ$messages)>0) cat("",summ$messages,sep="\n")

    options(width = w)
}

######################################################################
