source("methods.R")

whole.quant.file <- "~/workspace/ref/quant.sf"
overlap.description.file <- "~/workspace/ref/splicenum.txt"
original.sum <- LoadMainOriginal(whole.quant.file)
overlap.mapping <- GetSpliceSiteOverlapMapping(overlap.description.file)
HASH.OVERLAP <- IsoAndHisFriends(overlap.mapping)

par.dir <- "~/workspace/quants/0920/repeated5/"

quant.file.head.new <- "Name\tLength\tTPM\tNumReads"
SetQuantSFTableHead(quant.file.head.new)

quant.file.c <- GenerateQuantFileNames(5, paste(par.dir, "out5/oc", sep = ""))
quant.file.l <- GenerateQuantFileNames(5, paste(par.dir, "out5/ol", sep = ""))

quant.file.full.c <- GenerateQuantFileNames(5, paste(par.dir, "fullout5/oc", sep = ""))
quant.file.full.l <- GenerateQuantFileNames(5, paste(par.dir, "fullout5/ol", sep = ""))


l.full.re <- ExtractQuants(quant.file.full.l, by="Name", on="NumReads")
c.full.re <- ExtractQuants(quant.file.full.c, by="Name", on="NumReads")

full.re <- list(l.full.re, c.full.re)
obj.eb.full <- PrepareForEBseq(full.re, "~/workspace/quants/0920/bigRate.mapsf")
out.eb.full <- DoEBSeq(obj.eb.full)

######=======
l.re <- ExtractQuants(quant.file.l, by="Name", on="NumReads")
c.re <- ExtractQuants(quant.file.c, by="Name", on="NumReads")

full.re <- list(l.re, c.re)
obj.eb <- PrepareForEBseq(full.re, "~/workspace/quants/0920/bigRate.mapsf")
out.eb <- DoEBSeq(obj.eb)




####-------after ebseq

iso.de.ee.full <- names(out.eb.full$fc$RealFC)

rate.circ.de.ee.full <- CalculateRateOnSpliceSite(overlap.mapping, function(x,y){x/(x+y)},iso.de.ee.full, HASH.OVERLAP)

plot(rate.circ.de.ee.full, out.eb.full$fc$PostFC, ylim = c(0, 5))


iso.de.ee <- names(out.eb$fc$PostFC)
rate.circ.de.ee <- CalculateRateOnSpliceSite(overlap.mapping, function(x,y){x/(x+y)}, iso.de.ee, HASH.OVERLAP)

plot(rate.circ.de.ee, log(out.eb$fc$PostFC))
table(out.eb$fc$PostFC < 1)

# should remove those points with zero linear isoforms
num.overlap <- GetNumSpliceSite(overlap.mapping)
num.circ <- num.overlap$num.circ
num.linear.ss <- num.overlap$num.linear.splice.site


original.iso <- sapply(iso.de.ee, GetOriginalNum)
l.estimated.iso <- sapply(iso.de.ee, function(id, by, data){data[Order.by(id,by)]}, by = l.re$iso.list, data = l.re$data[,1])
l.estimated.iso[l.estimated.iso >0 & l.estimated.iso <1] = 1
c.estimated.iso <- sapply(iso.de.ee, function(id, by, data){data[Order.by(id,by)]}, by = c.re$iso.list, data = c.re$data[,1])
c.estimated.iso[c.estimated.iso>0 & c.estimated.iso<1] = 1

cores.circ.iso <- sapply(iso.de.ee, function(linear, mapping.hash){GetOriginalNum(mapping.hash[[linear]])}, mapping.hash = HASH.OVERLAP)


#sapply(list(num.circ, num.linear.ss, original.iso, l.estimated.iso, c.estimated.iso), length)

expression.hasCir <- data.frame(original.iso,l.estimated.iso,c.estimated.iso,rate.circ.de.ee, cores.circ.iso)

with(expression.hasCir, 
     plot(log10(original.iso), log10(c.estimated.iso))
) # end of with

with(expression.hasCir, 
    plot(log10(original.iso), log10(l.estimated.iso) )     
     )# end with


require(MASS)

fit.all <- lm(c.estimated.iso ~ . - l.estimated.iso, data = expression.hasCir)

stepAIC(fit.all, direction = "backward")

fit.final <- lm(c.estimated.iso ~ original.iso + rate.circ.de.ee, data = expression.hasCir)

p1 <- ggplot(expression.hasCir, aes(rate.circ.de.ee, log10(c.estimated.iso)))
p1 + geom_point() + geom_smooth(method = "loess")


+ geom_smooth(data = expression.hasCir[rate.circ.de.ee > 0.5,])
+ geom_smooth(data= expression.hasCir[rate.circ.de.ee < 0.5,])

qplot(data = expression.hasCir, original.iso,c.estimated.iso, geom = "point", log = "xy")
      
qplot()

AIC(fit.all, fit.final)

save(list = ls(),file = "~/Cloud/Work/reports/0921/main.rda")



with(expression.hasCir, 
{     plot(rate.circ.de.ee, log10(c.estimated.iso)
          ,pch=20, col = densCols(rate.circ.de.ee,c.estimated.iso
                                  #,colramp = colorRampPalette(c("gray","black"))
                                  ),
               )
    abline(v=.5, lwd=4, lty=2, col="orange")
} 
    ) # end with

with(expression.hasCir[rate.circ.de.ee>0.5 & c.estimated.iso >1,], 
     {
         lines(lowess(rate.circ.de.ee ,log10(c.estimated.iso), f=1),
               col="green", lwd = 3)
     }
     )

with(expression.hasCir[rate.circ.de.ee < 0.5 & c.estimated.iso >1,], 
     {
         lines(lowess(rate.circ.de.ee ,log10(c.estimated.iso), f=1),
               col="blue", lwd = 3)
     }
)

step(lm(c.estimated.iso ~ . - l.estimated.iso, data = expression.hasCir))


plot(rate.circ.de.ee, log(out.eb$fc$PostFC),
     pch=20, 
     col = densCols(rate.circ.de.ee, log(out.eb$fc$PostFC),
#                    colramp = colorRampPalette(c("gray","black"))
    ),
xlab = "circulation rate", ylab="log(fold change)", main =""
) # end of plot
abline(h=0, lwd=3, lty=1, col="orange")
