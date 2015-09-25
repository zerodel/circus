SetQuantSFTableHead <- function(default.table.head="Transcript	Length	TPM	RPKM	KPKM	EstimatedNumKmers	EstimatedNumReads"){
        assign("HEAD.COMMON", strsplit(default.table.head, split = "\t")[[1]], 
               pos = .GlobalEnv)
}

BASIC.RESULT.CLASS.TYPE <- "Circ.Simulation.Quant"

#HEAD.COMMON <- strsplit("Transcript	Length	TPM	RPKM	KPKM	EstimatedNumKmers	EstimatedNumReads", split = "\t")[[1]] 

######
ReadQuantDotSF <- function(dotsf){
    
    if (!exists("HEAD.COMMON")) {
        SetQuantSFTableHead()
    }
    
    if (!file.exists(dotsf)) {
        stop(paste("Quant file : ", dotsf," NOT exist !") )
    }
    
    tmp.table.data <- read.table(file = dotsf, comment.char = "#" , stringsAsFactors = F,
                                 colClasses = c("character", rep("numeric", times = length(HEAD.COMMON) - 1)))
    names(tmp.table.data) <- HEAD.COMMON # global_variable
    
    return(tmp.table.data)
}

GenerateQuantFileNames <- function(repeat.times, par.dir, target.quant.file="quant.sf"){
    repeat.patter.middle <- paste("/","0", 1:repeat.times, "/",sep = "")
    paste(par.dir, repeat.patter.middle, target.quant.file, sep = "")
}

Order.by  <- function(id, by) {
    if (any(by == id)) {
        return(which(by == id)[1])
    }
    return(0)
}

LoadMainOriginal <- function(sf.original, id="Transcript", data.of.interest="EstimatedNumReads"){

    
    origin.data <- ReadQuantDotSF(sf.original)
    
    result.data <- list()
    result.data$Transcript <- with(origin.data, get(id))
    result.data$NumberReads <- with(origin.data, get(data.of.interest))
    
    library(hash)
    HASH.MAIN.ORIGIN <<- hash()
    values(HASH.MAIN.ORIGIN, keys = result.data$Transcript) <- result.data$NumberReads

    return(data.frame(result.data))
}

ExtractQuants <- function(dot.sf.file.list, by="Transcript", on.which="EstimatedNumReads"){
    # hey,  result is a list contain a matrix and a dataframe
    result <- list()
    class(result) <- c(BASIC.RESULT.CLASS.TYPE)
    result$by.id <- by
    result$on.id <- on.which
    
    #check whether file 
    for (i in dot.sf.file.list) {
        if (!file.exists(i)) {
            stop(paste("No Such File : \"", i,"\"\n"))
        }
    }
    
    get.single.quant.sf <- function(sf, on){
        tmp.single.quant.sf <- ReadQuantDotSF(sf)
        
        if (!on %in% names(tmp.single.quant.sf)) {
            stop("given \"on\" not exists ")
        } # error rise when you give wrong "by" string. 
        
        return(tmp.single.quant.sf[on])
    }
    
    # get the data . 
    by.matrix <- sapply(dot.sf.file.list, get.single.quant.sf, on = by)
    on.matrix <- sapply(dot.sf.file.list, get.single.quant.sf, on = on.which)
    
    # check by.matrix is identical on each column 
    if (all(sapply(by.matrix, function(x, data){return(all(x == data))}, data = by.matrix[[1]])))
    {
        result$iso.list <- by.matrix[[1]]
        rm(by.matrix)        
    } else {
        stop("No the same index variables in data")
    }
    
    if (length(unique(sapply(on.matrix, length))) == 1) {
        result$data <- data.matrix(data.frame(on.matrix)) 
        colnames(result$data) <- paste("X", 1:length(on.matrix), sep = "")                            
    } else {
        stop("Not the same length of each quant file")
    }
    
    
    return(result)
} # end of whole function

AddOriginalInformation <- function(simultaion.result, assignment.polyest){
    #error happen when result is not complete
    if (!class(simultaion.result) == BASIC.RESULT.CLASS.TYPE) {
        stop("you should given a result produced by method :ExtractQuants")
    }
    isoforms <- simultaion.result$iso.list
    index.in.assignment <- sapply(isoforms, Order.by, by = assignment.polyest$Transcript)
    corresponding.quant <- ifelse(index.in.assignment > 0 , assignment.polyest$NumberReads[index.in.assignment],0)
    simultaion.result$original.reads <- corresponding.quant
    
    return(simultaion.result)
}

GetCorrespoindingGene <- function(isoform.shared, map.file){
    if (!file.exists(map.file)) {
        stop("No Such File")
    }
    table.iso.gene.num <- read.table(map.file, header = T,
                                     stringsAsFactors = F)
    
    index.of.iso <- sapply(isoform.shared, Order.by, by = table.iso.gene.num$Transcript)
    
    gene.corresponding <- table.iso.gene.num$gene[index.of.iso]
    
    return(gene.corresponding)
} # end of method : GetGeneIsoMap


Trim.EB.two.condition <- function(results){
    if (any(NULL == sapply(results, function(result){
        return(with(result, get("iso.list")))}))) 
    {
        stop("Some of the result not valid")
    }
    isoform.shared <- results[[1]]$iso.list
    for (result.single in results) {
        isoform.shared <- intersect(isoform.shared, result.single$iso.list)
    }    
    
    extracted <- lapply(results, 
                        function(result,ids){
                            index.in.data <- sapply(ids,Order.by, by = result$iso.list)
                            return(result$data[index.in.data,])}, 
                        ids = isoform.shared)
    
    value.return <- list()
    value.return$quant.mat <- Reduce(cbind, extracted)
    row.names(value.return$quant.mat) <- isoform.shared
    value.return$iso.shared <- isoform.shared
    value.return$divide.factor <- as.factor(rep(1:length(extracted), each = sapply(extracted, function(x){dim(x)[2]})))
    return(value.return)
} # end of function "Trim.EB.two.condition"

PrepareForEBseq <- function(results, map.file){
    trimed <- Trim.EB.two.condition(results)
    trimed$gene <- GetCorrespoindingGene(trimed$iso.shared, map.file)
    return(trimed)
}


DoEBSeq <- function(obj.eb, FDR.cer=0.05){
    library(EBSeq)
    iso.size.factor.list <- MedianNorm(obj.eb$quant.mat)
    trans.gene <- GetNg(obj.eb$iso.shared, obj.eb$gene)
    iso.ng.trun <- trans.gene$IsoformNgTrun
    
    iso.output <- EBTest(Data = obj.eb$quant.mat,
                         NgVector = iso.ng.trun,
                         Conditions = obj.eb$divide.factor,
                         sizeFactors = iso.size.factor.list,
                         maxround = 5)
    
    iso.output.res <- GetDEResults(iso.output, FDR = FDR.cer)
    
    iso.fc <- PostFC(iso.output)
    
    
    output = list()
    output$result <- iso.output
    output$res <- iso.output.res
    output$fc <- iso.fc
    return(output)
} # end of DoEBSeq


GetOriginalNum <- function(id, data = HASH.MAIN.ORIGIN){

    if (exists(as.character(substitute(data))) && is.hash(data)) {
        if (has.key(key = id, hash = data)){
            return(data[[id]])
        }
        else {
            return(0)
        }
    } else {
        stop("you should load the original data first")
    }
}



SumUpQuant <- function(isoform.ids, ...) {
    if (length(isoform.ids) > 0) {
        return(sum(c(sapply(isoform.ids, GetOriginalNum, ...)), na.rm = T))
    } else {
        # here we assume if there is no overlap isoform , take it as zero
        return(0)
    }
}

SplitLongString <- function(long.string, sep=","){
    return(strsplit(long.string, split = sep)[[1]])
}    


GetSpliceSiteOverlapMapping <- function(overlap.description.file, string.operation.function = SplitLongString) {
    overlap.circ.linear <- read.table(file = overlap.description.file, 
                                      stringsAsFactors = F, header = F,
                                      sep = "\t")
    
    overlap.mapping <- list()
    class(overlap.mapping) <- c("overlap.list.5and3")
    overlap.mapping$circRNA <- overlap.circ.linear$V1
    # here, we can get a list of list , I can not decide whether it is a good idea [15/09/19/:17:21]
    overlap.mapping$linear.5 <- sapply(overlap.circ.linear$V2, string.operation.function, sep = ",")
    overlap.mapping$linear.3 <- sapply(overlap.circ.linear$V3, string.operation.function, sep = ",")
    
    overlap.mapping$all <- mapply(intersect, overlap.mapping$linear.5, overlap.mapping$linear.3)
    
    return(overlap.mapping)
}

IsoAndHisFriends <- function(...){
    UseMethod("IsoAndHisFriends")
}

IsoAndHisFriends.overlap.list.5and3 <- function(overlap.mapping){
    some.ha <- hash()
    for (i in seq_along(overlap.mapping$circRNA)) {
        trans.5 = overlap.mapping$linear.5[i][[1]]
        trans.3 = overlap.mapping$linear.3[i][[1]]
        for (iso in trans.3) {
            if (!iso == overlap.mapping$circRNA[i])
            {some.ha[[iso]] <-  overlap.mapping$circRNA[i]}
        }
        for (iso in trans.5) {
            if (!iso == overlap.mapping$circRNA[i]) {
                some.ha[[iso]] <- overlap.mapping$circRNA[i]}
        }
    }
    
    return(some.ha)
}

GetNumSpliceSite <- function(overlap.mapping, ...){
    num.circ <- sapply(overlap.mapping$circRNA, GetOriginalNum, ...)
    name.circ <- overlap.mapping$circRNA
    num.overlap.all <- sapply(overlap.mapping$all, SumUpQuant, ...)
    num.overlap.all[num.overlap.all == 0] = 1
    
    output <- as.data.frame(list(name.circ, num.circ, num.overlap.all))
    names(output) <- c("name.circ", "num.circ", "num.linear.splice.site")
    return(output)
}


CalculateRateOnSpliceSite <- function(overlap.mapping, CalculateRate="/", iso.de.ee, HASH.OVERLAP, ...)
{

    num.splice <- GetNumSpliceSite(overlap.mapping, ...)
    num.circ <- num.splice$num.circ
    num.overlap.all <- num.splice$num.linear.splice.site
    
    GetQuantOfOverlap <- function(id, data){
        if (has.key(hash = HASH.OVERLAP, key = id)) {
            circ.id = HASH.OVERLAP[[id]]
            index.of.interest <- Order.by(id = circ.id, by = overlap.mapping$circRNA)
            return(data[index.of.interest])
        }
        return(0)
    }
    
    de.ee.circ.correspoding <- sapply(iso.de.ee, GetQuantOfOverlap, data = num.circ)
    de.ee.linear.correspoding <- sapply(iso.de.ee, GetQuantOfOverlap, data = num.overlap.all)
    
    return(mapply(FUN = CalculateRate, de.ee.circ.correspoding, de.ee.linear.correspoding))
}

SetPostDataHash <- function(...) {
    require(hash)
    if (exists("HASH.MAIN.POST")){
        clear(HASH.MAIN.POST)
        rm(list = c("HASH.MAIN.POST"))
    }
    
    HASH.MAIN.POST <<- hash(...)
}

