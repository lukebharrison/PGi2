create.edit.cost.matrix <-
function(pgi.tree,inf.params=NULL,verbosity=0) {
	require(e1071)
	if(is.null(inf.params)) {
		if(is.null(pgi.tree$inf.params)) stop("No inference parameters specified in the tree or explicitely. Stopping")
	} else {
		pgi.tree$inf.params<-inf.params
	}
	all_dev_seq_mat<-unique(rbind(create.all.dev.seq.mat(pgi.tree$seq.len,pgi.tree$inf.params$simultaneity,verbosity=verbosity),extract.seqs.from.tree(pgi.tree$tree,pgi.tree$seq.len,c())))
	dims<-dim(all_dev_seq_mat)[1]
	if(verbosity > 0) cat("Creating Edit-Cost (Step) Matrix, dimentions:",dims,"x",dims,".... \n")
	#editcostmat<-matrix(-1,dims,dims)
	editcostmat<-rep(-1,times=dims*dims)
	dim(editcostmat)<-c(dims,dims)
	if(verbosity == 1) pb<-txtProgressBar(min=0,max=dims,style=3)
	for(i in 1:dims) {
		for(j in 1:dims) {
			if(i == j) editcostmat[i,j]<-0
			if(j < i) {
				editcostmat[i,j]<-edit.cost(all_dev_seq_mat[i,],all_dev_seq_mat[j,],0,2,pgi.tree$inf.params$edit.cost.func,pgi.tree$seq.len)
				editcostmat[j,i]<-editcostmat[i,j]
			}   
		}
		if(verbosity == 1) setTxtProgressBar(pb,i)
	}
	cat("Created Edit-Cost (Step) Matrix dimentions:",dims,"x",dims,"\n")
	if(verbosity > 1) cat("Done.\n")
	if(verbosity == 1) close(pb)
	return(list(all_dev_seq_mat,editcostmat))
}

