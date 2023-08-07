pgi.to.phylo <-
function(pgi.tree) {
	if(class(pgi.tree) != "pgi.tree") stop("This tree isn't in pgi.tree format (use pgi.assemble.tree() or pgi.read.nexus() to convert/load)")
	require(ape)
	seq_len<-pgi.tree$seq.len
	ntaxa<-pgi.tree$ntaxa
	phylo<-rtree(ntaxa)
	edges<-c()
	tips<-c()
	nodes<-0
	edge.lengths<-c()
	taxa_ctr<-1
	node_ctr<-ntaxa
	node_lbl<-c()
	edge_lbl<-c()
	edge.a<-c()
	edge.d<-c()
	tip.seqs<-c()

 # check to see what kind of tree: filled, empty or censensus
	#this is an empty tree
	if(any(c("empty","filled") %in% pgi.tree$type)) {
		phylo_st<-list(edges,tips,tip.seqs,nodes,taxa_ctr,node_ctr)
		phylo_st<-.pgi.to.phylo.fill(pgi.tree$tree,phylo_st)
		phylo$edge<-phylo_st[[1]]
		phylo$Nnode<-phylo_st[[4]]
		phylo$tip.label<-phylo_st[[2]]
		phylo$tip.seqs<-phylo_st[[3]]
		phylo$edge.length<-NULL
	} else if(any(c("consensus","superconsensus") %in% pgi.tree$type)) {
		# this is a filled tree
		phylo_st<-list(edges,tips,tip.seqs,nodes,taxa_ctr,node_ctr,edge.lengths,node_lbl,edge_lbl,edge.a,edge.d)
		phylo_st<-.pgi.con.to.phylo.fill(pgi.tree$consensus,phylo_st,seq_len)
		phylo$edge<-phylo_st[[1]]
		phylo$Nnode<-phylo_st[[4]]
		phylo$tip.label<-phylo_st[[2]]
		phylo$edge.length<-phylo_st[[7]]
		phylo$node.label<-phylo_st[[8]]
		phylo$edge.label<-phylo_st[[9]]
		phylo$tip.seqs<-phylo_st[[3]]
		phylo_st[[10]]<-matrix(phylo_st[[10]],dim(phylo$edge)[1],seq_len,byrow=TRUE)
		phylo$edge.a<-phylo_st[[10]]
		phylo_st[[11]]<-matrix(phylo_st[[11]],dim(phylo$edge)[1],seq_len,byrow=TRUE)
		phylo$edge.d<-phylo_st[[11]]
	}
	return(phylo)
}

