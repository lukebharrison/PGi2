### Functions to convert a PGi Tree formatted object into a nexus data object
### Author: Luke B Harrison
### Currently simply ouputs a nexus.data object (as in APE) with event pairs
pgi.tree.to.nexus.data<-function(pgi.tree, event.pair=TRUE) {
	if(!event.pair) {
		stop("Error: Non event-pair sequence output not implemented (yet)\n")
	}
	out.list<-list()
	seq.mat<-extract.seqs.from.tree(pgi.tree)
	for(i in 1:dim(seq.mat)[1]) {
		t.rds<-seq.mat[i,]
		t.evpmat<-create.event.pair.matrix(t.rds)
 		t.evpseq<-c()
	  for(z in 1:dim(t.evpmat)[1]) {
		  for(j in z:dim(t.evpmat)[1]) {
  			t.evpseq<-c(t.evpseq,t.evpmat[j,z])
			}
		}
		out.list<-append(out.list,list(t.evpseq))
  }
	names(out.list)<-rownames(seq.mat)
	return(out.list)
}

