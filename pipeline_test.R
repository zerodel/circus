rm(list = ls())
source("methods.R")

result.path <- "/Users/zerodel/workspace/quants/qki/1000"

overlap.description.file <- "~/workspace/ref/splicenum.txt"
quant.file.head.new <- "Name\tLength\tTPM\tNumReads"
SetQuantSFTableHead(quant.file.head.new)

whole.quant.file <- "~/workspace/quants/qki/out/oc/quant.sf"
original.sum <- LoadMainOriginal(whole.quant.file, id = "Name", data.of.interest = "NumReads")
overlap.mapping <- GetSpliceSiteOverlapMapping(overlap.description.file)
HASH.OVERLAP <- IsoAndHisFriends(overlap.mapping)

# read the quant result . 


l.qki <- ExtractQuants(GenerateQuantFileNames(1, paste(result.path, "/out/l",sep = "")),by="Name", on="NumReads")
c.qki <- ExtractQuants(GenerateQuantFileNames(1, paste(result.path, "/out/c",sep = "")), by="Name", on="NumReads")



# calculate rate
# 
# CalRate <- function(x,y){
#     if (y == 0){
#         
#     }
# } 

circ.iso <- Filter(function(id){substr(id, 1,3) == "chr"},
                              as.character(c.qki$iso.list))
length(unique(circ.iso))

rate.qki <- CalculateRateOnSpliceSite(overlap.mapping, function(x,y){x/(x+y)}, circ.iso, HASH.OVERLAP)

num.circ <- sapply(circ.iso, GetOriginalNum)

SetPostDataHash(keys = c.qki$iso.list , values = c.qki$data)


rate.qki.post <- CalculateRateOnSpliceSite(overlap.mapping, function(x,y){x/(x+y)}, circ.iso, HASH.OVERLAP, data = HASH.MAIN.POST)
# na , maybe from 0/0

circ.rate.na <- is.na(rate.qki.post)

hist(num.circ[circ.rate.na])
summary(num.circ[circ.rate.na])

circ.rate.postive <- (! circ.rate.na) & ( num.circ > 0)

plot(log(num.circ[circ.rate.postive]), rate.qki.post[circ.rate.postive])
