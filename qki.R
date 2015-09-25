overlap.description.file <- "~/workspace/ref/splicenum.txt"

quant.file.head.new <- "Name\tLength\tTPM\tNumReads"
SetQuantSFTableHead(quant.file.head.new)

whole.quant.file <- "~/workspace/quants/qki/out/oc/quant.sf"
original.sum <- LoadMainOriginal(whole.quant.file, id = "Name", data.of.interest = "NumReads")
overlap.mapping <- GetSpliceSiteOverlapMapping(overlap.description.file)
HASH.OVERLAP <- IsoAndHisFriends(overlap.mapping)



l.qki <- ExtractQuants(c("~/workspace/quants/qki/out/ol/quant.sf"),by="Name", on="NumReads")
c.qki <- ExtractQuants(c("~/workspace/quants/qki/out/oc/quant.sf"),by="Name", on="NumReads")

qki.re <- list(l.qki, c.qki)

obj.eb.qki <- PrepareForEBseq(full.re, "~/workspace/quants/0920/bigRate.mapsf")
out.eb.qki <- DoEBSeq(obj.eb.full)



linear.isoforms.qki <- Filter(function(id){!substr(id, 1,3) == "chr"},
                              as.character(original.sum$Transcript))



rate.qki <- CalculateRateOnSpliceSite(overlap.mapping, function(x,y){x/(x+y)}, linear.isoforms.qki, HASH.OVERLAP)

num.qki <- sapply(linear.isoforms.qki, GetOriginalNum)
#num.qki[num.qki < 1] = 1


save(list=c("rate.qki", "num.qki"), file="qki.rda")

plot(rate.qki, log(num.qki), pch = 20, col = "gray")
abline(h=5,lwd=2,lty=3)


lines(lowess(rate.qki, log(num.qki),f = 1), lwd = 3, col = "blue")

tmp.lm <- data.frame(log(num.qki), rate.qki)
names(tmp.lm) <- c("num","rate")
summary(lm( num ~ rate,data = tmp.lm))
