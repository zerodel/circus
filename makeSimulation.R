rm(list = ls())
source("methods.R")

quant.file.head.new <- "Name\tLength\tTPM\tNumReads"
SetQuantSFTableHead(quant.file.head.new)


whole.quant.file <- "~/workspace/quants/qki/out/oc/quant.sf"

overlap.description.file <- "~/workspace/ref/splicenum.txt"
original.sum <- LoadMainOriginal(whole.quant.file, id = "Name", data.of.interest = "NumReads")
overlap.mapping <- GetSpliceSiteOverlapMapping(overlap.description.file)
HASH.OVERLAP <- IsoAndHisFriends(overlap.mapping)


# lc.num.splice.site.global <- GetNumSpliceSite(overlap.mapping)
# 
# ratio.raw <- mapply(function(x, y){return(x/(x + y))}, lc.num.splice.site.global$num.circ, lc.num.splice.site.global$num.linear.splice.site)
# 
# summary(ratio.raw)

# table(ratio.raw == 0 )
# 
# ratio.raw[ratio.raw ==0] = 0.001
# ratio.adv <- log(ratio.raw/(1-ratio.raw))
# hist(ratio.adv, breaks = 50)
# 
# MimicNorm <- function(x){
#     return(rnorm(n = length(x), mean = mean(x, na.rm = T, trim = 0.45), sd = sd(x, na.rm = T)))
# }
# 
# qqplot(ratio.adv, MimicNorm(ratio.adv))
# 

# if we want sampling as the orignal distribution . try rank 



index.to.pickout <- sample(1: length(overlap.mapping$circRNA), size = 1000)

pool.iso.all <- Reduce(f = union, x = overlap.mapping$all[index.to.pickout], init = overlap.mapping$circRNA[index.to.pickout])
length(unique(pool.iso.all))

iso.num <- sapply(pool.iso.all, GetOriginalNum)

length(iso.num) == length(pool.iso.all)
summary(iso.num)

simu.assignment <- data.frame(pool.iso.all, iso.num)
names(simu.assignment) <- c("Transcript", "NumReads")

write.table(x=simu.assignment, file = "~/workspace/quants/assignment/all1000/iso.txt"
        , sep = "\t", row.names = F, col.names = T)

