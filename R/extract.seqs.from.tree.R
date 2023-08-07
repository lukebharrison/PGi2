extract.seqs.from.tree <-
function(tree,seq_len=NULL,seqs=c()) {
 if(is.null(seq_len)) {
	## First iteration	
	seq_len<-tree$seq.len
	tree<-tree$tree
 } 
 t.seqs<-c()
 for(i in 2:length(tree)) {
  if(length(tree[[i]][[1]]) > 1) {
   #tip
	 t.rownames<-rownames(t.seqs)
   t.seqs<-rbind(t.seqs,tree[[i]][[1]][1:seq_len])
	 rownames(t.seqs)<-c(t.rownames,rownames(tree[[i]][[1]]))
  } else {
   #node
   ret.seqs<-extract.seqs.from.tree(tree[[i]],seq_len=seq_len)
	 t.rownames<-rownames(t.seqs)
   t.seqs<-rbind(t.seqs,ret.seqs)
	 rownames(t.seqs)<-c(t.rownames,rownames(ret.seqs))
  }
 }
 return(t.seqs)
}

