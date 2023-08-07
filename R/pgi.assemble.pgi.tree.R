pgi.assemble.pgi.tree <-
function(old.pgi.tree) {
	pgi.tree<-list(type="empty",seq.len=.pgi.old.get.seq.len(old.pgi.tree),nnodes=.pgi.old.nnodes(old.pgi.tree),ntaxa.=pgi.old.ntaxa(old.pgi.tree),tree=old.pgi.tree)
	class(pgi.tree)<-"pgi.tree"
	return(pgi.tree)
}

