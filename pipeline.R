source("methods.R")

par.dir <- "~/workspace/quants/"

l.sf <- GenerateQuantFileNames(5, paste(par.dir,"l",sep = ""))
c.sf <- GenerateQuantFileNames(5, paste(par.dir,"c",sep = ""))

l.result <- ExtractQuants(l.sf)
c.result <- ExtractQuants(c.sf)

Polyester.res <- list(l.result,c.result)
obj.eb <- PrepareForEBseq(Polyester.res, "~/workspace/quants/mixed.txt")
out.eb <- DoEBSeq(obj.eb)

#CalcRatioByCircID <- function(whole.quant.file, overlap.description.file){
whole.quant.file <- "~/workspace/ref/quant.sf"
overlap.description.file <- "~/workspace/ref/splicenum.txt"

original.sum <- LoadMainOriginal(whole.quant.file)

overlap.mapping <- GetSpliceSiteOverlapMapping(overlap.description.file)

HASH.OVERLAP <- IsoAndHisFriends(overlap.mapping)

lc.num.splice.site.global <- GetNumSpliceSite(overlap.mapping)

ratio.raw <- mapply(function(x, y){return(x/(x + y))}, lc.num.splice.site.global$num.circ, lc.num.splice.site.global$num.linear.splice.site)

summary(ratio.raw)
table(ratio.raw > 0.5)


isoform.each.circ <- mapply(intersect, overlap.mapping$linear.5, overlap.mapping$linear.3)

isoform.pool <- overlap.mapping$circRNA[ratio.raw > 0.5]

# for (iso.circ.index in which(ratio.raw > 0.5)) {
#     isoform.pool <- c(isoform.pool, isoform.each.circ[iso.circ.index])
# }

pool.overlapping.isoforms <- Reduce(f = union, x = isoform.each.circ[ratio.raw > 0.5], init = isoform.pool)

isoforms.control <- setdiff(original.sum$Transcript, pool.overlapping.isoforms)

set.seed(1010)
pool.control.isoforms <- sample(isoforms.control, size = length(pool.overlapping.isoforms))

pool.final.isoforms <- unique(union(pool.control.isoforms, pool.overlapping.isoforms))

num.overlap.isoforms <- sapply(pool.overlapping.isoforms, GetOriginalNum)
data.overlap <- data.frame(pool.overlapping.isoforms, num.overlap.isoforms)
names(data.overlap) <- c("Transcript", "NumReads")
write.table(data.overlap, file = "~/workspace/quants/0920/iso.txt",
                   sep = "\t", col.names = T, row.names = F)

num.final.isoforms <- sapply(pool.final.isoforms, GetOriginalNum)

output <- data.frame(pool.final.isoforms, num.final.isoforms)
names(output) <- c("Transcript", "NumReads")

write.table(output, file = "~/workspace/quants/0921/iso.txt",
      sep = "\t", col.names = T, row.names = F)

aa <- read.table("~/workspace/quants/0920/mixed.txt", header = T, stringsAsFactors = F)

######--------------after EBseq --------------------------------
iso.de.ee <- names(out.eb$fc$RealFC)


rate.circ.de.ee <- CalculateRateOnSpliceSite(overlap.mapping, function(x,y){x/(x + y)},iso.de.ee, HASH.OVERLAP)




# this is the evaluation and plot step 


plot(y = out.eb$fc$RealFC, x = rate.circ.de.ee, ylim = c(0,3), xlim = c(0,5))

plot(y = out.eb$fc$RealFC[rate.circ.de.ee > 0 ], x = rate.circ.de.ee[rate.circ.de.ee > 0],
     xlim = c(0,5), ylim = c(0,5))
lines(lowess(y = out.eb$fc$RealFC[rate.circ.de.ee > 0 ], x = rate.circ.de.ee[rate.circ.de.ee > 0]))
     
tmp.fit.data <- data.frame(out.eb$fc$RealFC[rate.circ.de.ee > 0] - 1,rate.circ.de.ee[rate.circ.de.ee > 0])
names(tmp.fit.data) <- c("fc", "rate")
aaa <- lm(fc ~ rate, data = tmp.fit.data)
summary(aaa)

system.time(CommonIso2(overlap.mapping$linear.5[1:200], overlap.mapping$linear.3[1:200]))

system.time(CommonIso(overlap.mapping$linear.5[1:200], overlap.mapping$linear.3[1:200]))
