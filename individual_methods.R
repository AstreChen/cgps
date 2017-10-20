library(limma)
library(Biobase)
library(EnrichmentBrowser)
library(GSVA)
library(dplyr)

#### DATA WE NEED ###
#* GENE sets
#* GENE NETWORK *
#* svm *
#* train data and test data * 

get.gs.grn.from.kegg <- function(spe){
    ### either use file or use spe to get KEGG PATHWAY ###
    message(paste("loading KEGG PATHWAY of",spe))
    pwys <- EnrichmentBrowser::download.kegg.pathways(spe)
    grn <- EnrichmentBrowser::compile.grn.from.kegg(pwys)
    gsl <- EnrichmentBrowser::get.kegg.genesets(pwys)
    return (list(gsl=gsl, grn=grn))
}
read.gs.from.file <- function(gmtfile){
  if(file.exists(gmtfile)) stop(paste0(gmtfile," file not exist."))
  tb <- read.delim(gmtfile)
  gsl <- lapply(seq_len(nrow(tb)), FUN=function(x){return (as.character(tb[x,]))})
  names(gsl) = rownames(tb)
  return (gsl)
}
  
# write genesets to file in GMT format
write.gmt <- function(gs, gmt.file)
{
  ## collapse geneset members to one tab separated string 
  gs.strings <- sapply(gs, function(x) paste(x,collapse="\t"))
  
  ## paste an not annotated second column (artifact of gmt format)
  ann <- paste(names(gs), rep(NA,length(gs)), sep="\t")
  
  ## paste all together
  all <- paste(ann, gs.strings, sep="\t")
  
  ## collapse all together to a single newline separated string
  all.str <- paste(all, collapse="\n")
  all.str <- paste(all, "\n", sep="")
  
  ## write the gs in gmt format
  cat(all.str, file=gmt.file, sep="")
}


read.exp <- function(exp){
    expdf = read.delim(exp, stringsAsFactors=F, header=T)
    expdata = as.matrix(expdf[,2:ncol(expdf)])
    rownames(expdata) <- expdf[,1]
    return(expdata)
}

read.phe <- function(phe.f){
    phe <- read.table(phe.f)
    phe <- phe$V1
    samples <- colnames(expdata)
    pdata <- data.frame(SAMPLE= samples, GROUP=phe, stringsAsFactors=F)
    rownames(pdata) <- samples
    pdata <- Biobase::AnnotatedDataFrame(pdata)
    return (pdata)
}

####### NORMALIZE  AND  MAKE EXSET #######
norm.exp <- function(expdata,pdata,dtype){
    ### need to normalize ###
    if (dtype=='ma'){
        expdata <- limma::normalizeBetweenArrays(expdata, method="quantile")
    } else if(dtype =='rseq'){
        cpm.count <- cpm(expdata,log=T, prior.count=3)
        expdata <- cpm.count
    }else stop("\'dtype\' must be either \'rseq\' or \'ma\' ")
    exset <- Biobase::ExpressionSet(assayData=exp.matrix, phenoData=phdata)
    return (exset)
}


norm.de.exp <- function(exset, datatype){
    if (datatype=='ma'){
        exprs(exset) <- limma::normalizeBetweenArrays(expdata, method="quantile")
        exset <- de.ana(exset, de.method="limma")  ### only saved in fData, but not change EXPRESSION DATA
    }
    else if (datatype=='rseq') {
        is.wholenumber <- function(x, tol=.Machine$double.eps^0.5) abs(x-round(x)) < tol
        auto.detect.normalize <- function(expr)
            ifelse(all(is.wholenumber(expr), na.rm=TRUE), "unnormed", "normed")
        if (auto.detect.normalize(exprs(exset))=="unnormed"){
            cpm.count <- cpm(exprs(exset),log=T, prior.count=3)
            exprs(exset) <- cpm.count
        }
        else{
            exprs(exset) <- log2(exprs(exset) + 3.0)
        }
        exprs(exset) <- limma::normalizeBetweenArrays(expdata, method="quantile")
        exset <- de.ana(exset,de.method="limma")
    }
    else {
        stop('none acceptable data type.')
    }
    return (exset)
}


get.exset <- function(expdata,pdata){
    ### no need to normalize ###
    phdata <- Biobase::AnnotatedDataFrame(phdata)
    exset <- Biobase::ExpressionSet(assayData=exp.matrix, phenoData=phdata)
    return (exset)
}


###### RUN EVERY METHOD ####
run.plage <- function(exset,gs){
    pdata <- Biobase::pData(exset)
    design <- model.matrix(~pdata$GROUP)
    plage.res <- GSVA::gsva(exset, gs, method="plage", rnaseq=F)
    fit <- limma::lmFit(plage.res$es.obs, design)
    fit <- limma::eBayes(fit)
    plage.tbl <- limma::toptable(fit, coef=2, number=length(gs))
    plage.tbl <- dplyr::add_rownames(plage.tbl,"GENE.SET")
    #write.table(plage.tbl, file=paste0(path,"_","plage",".tsv"), sep="\t",row.names=F)
    return(plage.tbl)
}
run.gage <- function(exset, gsets){
    ### phenotype word ###
    expdata = Biobase::exprs(exset)
    phdata = Biobase::pData(exset)
    ctrl = grep('0',phdata$GROUP)
    case = grep('1',phdata$GROUP)

    if (length(ctrl)==length(case)){
        gage.res <- gage::gage(exprs=expdata, gsets=gsets, ref=ctrl, samp=case, same_dir=F)
    }else{
        gage.res <- gage::gage(exprs=expdata, gsets=gsets, ref=ctrl, samp=case, same_dir=F, compare='unpaired')
    }

    #### GET UP AND DOWN RESULTS ###
    grdf <- data.frame(gage.res$greater[, 3:4])
    colnames(grdf) <- c("pval.gr","qval.gr")
    grdf$GENE.SET <- rownames(grdf)
    ledf <- data.frame(gage.res$less[, 3:4])
    colnames(ledf) <- c("pval.le","qval.le")
    ledf$GENE.SET <- rownames(ledf)
    ### MERGE ALLDF ###
    alldf <- merge(grdf, ledf)
    alldf$pval.min <- apply(alldf[,c(2,4)],1,min)
    alldf$qval.min <- apply(alldf[,c(3,5)],1,min)
    alldf$direction <- "up"
    alldf$direction[with(alldf, pval.min==pval.le)] <- "down"
    ### ORDER BY P-VALUE ###
    alldf <- alldf[order(alldf$qval.min),]
    #### FILTER THE NULL GENE SETS ###
    alldf <- alldf[complete.cases(alldf),]
    #write.table(alldf, file=paste0("gage",".tsv"), sep="\t",row.names=F)
    return (alldf)
}
run.gsea <- function(exset,gene.sets){
    res.tbl <- EnrichmentBrowser::sbea(method='gsea', eset=exset, gs=gene.sets, perm=1000)$res.tbl
    return(res.tbl)
}
run.safe <- function(exset,gene.sets){
    res.tbl <- EnrichmentBrowser::sbea(method='safe', eset=exset, gs=gene.sets, perm=1000)$res.tbl
    return(res.tbl)
}
run.globaltest <- function(exset,gene.sets){
    res.tbl <- EnrichmentBrowser::sbea(method='globaltest', eset=exset, gs=gene.sets, perm=1000)$res.tbl
    return(res.tbl)
}
run.gsa <- function(exset,gene.sets){
    res.tbl <- EnrichmentBrowser::sbea(method='gsa', eset=exset, gs=gene.sets, perm=1000)$res.tbl
    return(res.tbl)    
}

run.padog <- function(exset,gene.sets){
    res.tbl <- EnrichmentBrowser::sbea(method='padog', eset=exset, gs=gene.sets, perm=1000)$res.tbl
    return(res.tbl)   
}

run.cepa <- function(exset, gene.sets,grn){
    res.tbl <- EnrichmentBrowser::nbea(method='cepa', eset=exset, gs=gene.sets,grn=grn, perm=1000)$res.tbl
    return(res.tbl) 
}
run.ganpa <- function(exset, gene.sets, grn){
    res.tbl <- EnrichmentBrowser::nbea(method='ganpa', eset=exset, gs=gene.sets,grn=grn, perm=1000)$res.tbl
    return(res.tbl)    
}

##### GET EVERY METHODS RESULTS ####


#####  USE SVM TO PREDICT ####



##### give the result ######
#R score
#sort by R score 
#give the result table 



