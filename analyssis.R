# set the working directory
#setwd("~/workspace/quants/")


# identifiers like "NUMBER.BER" are global variables. 

# get the corresponding file in subdirs . 



# get the transcipt length
original.quant <- ReadQuantDotSF("~/workspace/ref/quant.sf")


trans.DEseq <- names(iso.output.res$Status) # get all the isoform id . 

trans.DEseq.index <- sapply(trans.DEseq, function(tran, data){which(data == tran) }, data = original.quant$Transcript)

trans.DEseq.len <- original.quant$Length[trans.DEseq.index]
trans.DEseq.expression <- original.quant$EstimatedNumReads[trans.DEseq.index]

# do a logistic regression
# make a vector of T/F

is.DE.expressed <- rep(F, length(trans.DEseq.expression))
is.DE.expressed[iso.output.res$Status == "DE"] = T

index.de.and.ee <- iso.output.res$Status == "DE" | iso.output.res$Status == "EE"

pick.from.original <- data.frame(is.DE.expressed[index.de.and.ee], trans.DEseq.len[index.de.and.ee], trans.DEseq.expression[index.de.and.ee])

fitDE.binomalmodel <- glm(is.DE.expressed ~ trans.DEseq.len + trans.DEseq.expression, 
                          family = binomial(),data = pick.from.original)

# this shows the length and reads number assigned to the simulator software has no significant impact on the expression on the transcripts. 

# let me add some variables ...
# for example, the overlap information ? 

overlaps <- read.table(file = "~/workspace/ref/splicenum.txt",
                         header = F, sep = "\t",
                         stringsAsFactors = F)
names(overlaps) <- c("circ", "end5", "end3")

overlap.linear.circ <- read.table(file = "~/workspace/ref/linearOverlapCirc.txt", header = T, sep = "\t", stringsAsFactors = F)

trans.de.overlap.index = sapply(trans.DEseq, AllInAll, data = overlap.linear.circ$linear)

trans.de.hasoverlap = trans.de.overlap.index > 0

trans.de.share = sapply(trans.de.overlap.index, function(overlap.id, data, where){
    if (overlap.id > 0) {
        return(any(where == data[overlap.id]))
            }else{
        return(F)
    }
    } , data = overlap.linear.circ$circ, where = cRnall$Transcript)

# since no other circRNA isoform overlaps with those linear isoform ...

picked <- cbind(pick.from.original, trans.de.share[index.de.and.ee])

fit.overlap <- glm(formula = is.DE.expressed ~ trans.DEseq.len + trans.DEseq.expression + trans.de.share,
                   family = binomial(),
                   data = picked)
summary(fit.overlap)

# this result shows that length and share is important , so we make a new model.
fit.overlap.len <- glm(formula = is.DE.expressed ~ trans.DEseq.len + trans.de.share,
                   family = binomial(),
                   data = picked)
summary(fit.overlap.len)


fit.overlap.only.share <- glm(formula = is.DE.expressed ~ trans.de.share,
                          family = binomial(),
                          data = picked)
summary(fit.overlap.only.share)

# these two model has no significant difference
anova(fit.overlap,fit.overlap.len, test = "Chisq")

anova(fit.overlap,fit.overlap.only.share, test = "Chisq")


save(list = ls(), file = "~/Cloud/Work/reports/0905/data.rda")

sure.de <- is.DE.expressed[index.de.and.ee]
sure.share <- trans.de.share[index.de.and.ee]

de.share.table <- matrix(c(
length(which(sure.de & sure.share)),
length(which(sure.de & !sure.share)),
length(which(!sure.de & sure.share)),
length(which(!sure.de & !sure.share))),
nrow = 2, ncol = 2, byrow = T,
dimnames = list(c("DE","EE"),c("share","pure-linear"))
)

chisq.test(de.share.table, correct = T)



tmpC <- data.frame(cRnall$X1[first.in.second][index.de.and.ee], cRnall$NumRead[first.in.second][index.de.and.ee])
names(tmpC) <- c("out", "input")

fit.c <- lm(out ~ input, data = tmpC)
summary(fit.c)

tmpL <- data.frame(lRnall$X1[index.de.and.ee], lRnall$NumRead[index.de.and.ee])
names(tmpL) <- c("out", "input")
fit.l <- lm(out ~ input, data = tmpL)
summary(fit.l)




# get linear id , trans.DEseq
#
#lrnall.corres.circ <- overlap.linear.circ$circ[trans.de.overlap.index]
# library(hash)
# allha <- hash(key = original.quant$Transcript,
#      values = original.quant$EstimatedNumReads )
# clear(allha)
# rm(allha)


linear.fc <- names(iso.foldchange$PostFC)
linear.de.or.ee <- trans.DEseq[index.de.and.ee]


# ! important
linear.worth.analysis <- intersect(linear.fc, linear.de.or.ee)

GetCorrespondingCirc <- function(iso, map.list.linear, map.list.circ, picked.circ,data){
    index.in.quant.sf <- 0
    if (any(map.list.linear == iso))
    {
        index.in.quant.sf <- which(map.list.linear == iso)
        circ.id = map.list.circ[index.in.quant.sf]
        if (any(picked.circ == circ.id))
        {
            return(data[which(picked.circ == circ.id)])
        }
    } 
    return(NA)

}

linear.corres.circ.expression <- sapply(linear.worth.analysis,
                                        GetCorrespondingCirc,
                                        map.list.linear = overlap.linear.circ$linear,
                                        map.list.circ = overlap.linear.circ$circ,
                                        picked.circ = cRnall$Transcript,
                                        data = cRnall$NumRead
)

index.circ.decent.expressed <- !is.na(linear.corres.circ.expression) & linear.corres.circ.expression > 100

linear.decent <- linear.worth.analysis[index.circ.decent.expressed]


# get the expression of linear RNA
index.linear.decent.expression <- sapply(linear.decent, AllInAll, data = cRnall$Transcript)
linear.decent.expressed <- cRnall$NumRead[index.linear.decent.expression]


index.fc.linear.decent <- sapply(linear.decent, AllInAll, data = linear.fc)

fc.linear.decent <- iso.foldchange$PostFC[index.fc.linear.decent]
circ.decent.expressed <- linear.corres.circ.expression[index.circ.decent.expressed]

data4lm <- data.frame(fc.linear.decent, circ.decent.expressed, linear.decent.expressed/circ.decent.expressed)
names(data4lm) <- c("fc", "circ", "ratiol2c")

fit.fc.circ <- lm(fc ~ circ, data = data4lm)
summary(fit.fc.circ)



# ok , [15/09/10/:13:50] we check the overlap circ and linear 
#first , list of isoforms sperated by comma should be retrived . and get the expression information . 

GetLinearNumberOnSplieSpot <- function(comma.divide.str, data.isoform, data.num){
    isoform.list <- strsplit(comma.divide.str, split = ",")[[1]]
    index.isoform <- sapply(isoform.list, AllInAll, data = data.isoform)
    return(data.num[index.isoform])
}



