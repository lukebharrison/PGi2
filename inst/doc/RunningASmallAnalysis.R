### R code from vignette source 'RunningASmallAnalysis.Rnw'

###################################################
### code chunk number 1: RunningASmallAnalysis.Rnw:16-17
###################################################
options(width = 80, prompt = "> ")


###################################################
### code chunk number 2: RunningASmallAnalysis.Rnw:74-75
###################################################
writeLines(readLines("Velhagen1997.nex"))


###################################################
### code chunk number 3: RunningASmallAnalysis.Rnw:83-85
###################################################
library(pgi2)
velhagen <- pgi.read.nexus("Velhagen1997.nex")


###################################################
### code chunk number 4: RunningASmallAnalysis.Rnw:91-92
###################################################
summary(velhagen)


###################################################
### code chunk number 5: RunningASmallAnalysis.Rnw:95-96
###################################################
plot(velhagen,show.tip.seq=TRUE)


###################################################
### code chunk number 6: RunningASmallAnalysis.Rnw:109-110
###################################################
velhagen.con.trees<-pgi(velhagen,nruns=c(4,2),inf.params=list(heuristic="pgi",cycles=100,replicates=100,ret.anc.seq=100,edit.cost.func="parsimov",simultaneity=TRUE),con.params=list(con.type="semi-exhaustive",semi.ex.con.max.n=5000,edit.cost.func="parsimov"),verbosity=1)


###################################################
### code chunk number 7: RunningASmallAnalysis.Rnw:116-117
###################################################
summary(velhagen.con.trees)


###################################################
### code chunk number 8: RunningASmallAnalysis.Rnw:123-124
###################################################
plot(velhagen.con.trees[[3]],show.tip.seq=TRUE,show.anc.seq=TRUE)


###################################################
### code chunk number 9: RunningASmallAnalysis.Rnw:131-132
###################################################
plot(velhagen.con.trees,show.tip.seq=FALSE,show.anc.seq=FALSE)


###################################################
### code chunk number 10: RunningASmallAnalysis.Rnw:141-142
###################################################
velhagen.supercon<-pgi.supercon(velhagen.con.trees,tol=0,verbosity=1)


###################################################
### code chunk number 11: RunningASmallAnalysis.Rnw:144-145
###################################################
plot(velhagen.supercon,show.tip.seq=TRUE,show.anc.seq=TRUE)


###################################################
### code chunk number 12: RunningASmallAnalysis.Rnw:151-152
###################################################
plot(velhagen.supercon,seq.het.thres=0,print.support=TRUE,show.anc.seq=TRUE,show.tip.seq=TRUE)


