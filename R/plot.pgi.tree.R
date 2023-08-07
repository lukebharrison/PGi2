plot.pgi.tree <-
function(pgi.tree,seq.het.thres=0.5,print.support=FALSE,show.anc.seq=FALSE,show.tip.seq=FALSE,show.bls=FALSE,...) {
	phylo_t<-pgi.to.phylo(pgi.tree)
	if(show.tip.seq == TRUE) {
		for(i in 1:length(phylo_t$tip.labe)) phylo_t$tip.label[i]<-paste(phylo_t$tip.label[i]," (",phylo_t$tip.seqs[i],")",sep="")
	}
	if(show.anc.seq) nodepos<-1 else nodepos<-2
	p<-plot.phylo(phylo_t,show.node.label=F,type="phylogram",node.pos=nodepos,use.edge.length=show.bls,...)
	if(show.anc.seq == T) nodelabels(phylo_t$node.label,cex=0.6,frame="rect",adj=c(0,0.5),font=2,bg="white")
	#if(pgi.tree$type == "superconsensus") print.support<-TRUE
 	if(seq.het.thres < -1 || seq.het.thres > 1) stop("Error: Please specify a consensus strictness for sequence heterochronies between 0 and 1, (-1 disables printing of heterochronies, 0 = display all, 0.5 = majority rule, 1 = strict)")
	if(seq.het.thres == 0) {
		cat("Sequence heterochrony threshold is 0, automatically printing support values\n")
		print.support<-TRUE
	}
	if(seq.het.thres >= 0 && any(pgi.tree$type==c("consensus","superconsensus"))) {
		# therefore print heterochronies
		# code block borrow from ape::edgelabels
		lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		sel <- 1:dim(lastPP$edge)[1]
		subedge <- lastPP$edge
		if (lastPP$type == "phylogram") {
			if (lastPP$direction %in% c("rightwards", "leftwards")) {
				#XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[,2]])/2
				XX<-lastPP$xx[subedge[,1]]
				YY <- lastPP$yy[subedge[, 2]]
			} else {
				XX <- lastPP$xx[subedge[, 2]]
				YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[,2]])/2
			}
		}
		# end code block borrowed from ape::edgelabels

		if(print.support) {
			edge.labs.a<-paste("A:",apply(phylo_t$edge.a,1,
			function(x) { 
				ret<-c()
				for(i in 1:length(x)) {
					if(x[i] > seq.het.thres) {
						ret<-c(ret,paste(i,"(",format(x[i],digits=2),")",sep=""))
					} else {
						ret<-c(ret,character(0))
					}
				}
			return(paste(ret,collapse=","))
			}))
			edge.labs.d<-paste("D:",apply(phylo_t$edge.d,1,
			function(x) { 
				ret<-c()
				for(i in 1:length(x)) {
					if(x[i] > seq.het.thres) {
						ret<-c(ret,paste(i,"(",format(x[i],digits=2),")",sep=""))
					} else {
						ret<-c(ret,character(0))
					}
				}
				return(paste(ret,collapse=","))
			}))
		} else {
			edge.labs.a<-paste("A:",apply(phylo_t$edge.a,1,function(x) { return(paste(which(x > seq.het.thres),collapse=",")) }))
			edge.labs.d<-paste("D:",apply(phylo_t$edge.d,1,function(x) { return(paste(which(x > seq.het.thres),collapse=",")) }))
		}
		text(XX+0.1,YY,edge.labs.a,adj=c(0,-0.4),cex=0.53)
		text(XX+0.1,YY-0.05,edge.labs.d,adj=c(0,1),cex=0.53)
	}
}

