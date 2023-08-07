pgi.read.nexus <-
function(file=NULL,interactive=FALSE,strip.single.quotes=TRUE) {
	if(is.null(file)) {
		if(interactive) {
			cat("\n***Load a dataset from file***\n\n")
			cat("Type the file name of the NEXUS file in the working directory (",getwd(),") or the fully qualified file name and press <Enter> or type only <Enter> to open a GUI prompt \n",sep="")
			file<-readline(prompt="File name > ")
			if(file == "") {
				require(tcltk)
				file<-tk_choose.files(caption="Select NEXUS file with a PGi tree")
			}
		}
	}
	library(ape)
	tree<-read.nexus(file)
	data<-read.nexus.data(file)
	## Strip single quotes, if any
	if(strip.single.quotes) names(data)<-gsub("'","",names(data),fixed=TRUE)
	for(i in 1:length(data)) {
		data[[i]]<-as.integer(.trans.chars(data[[i]]))
		data[[i]]<-compress.rds(data[[i]])
	}
	ntaxa<-length(tree$tip.label)
	seq_len<-length(data[[1]])
	#edges<-tree[1][[1]]
	edges<-tree$edge
	fin_tree<-.pgi.build.tree(tree,edges,data,ntaxa,seq_len)
	pgi.tree<-list(type="empty",seq.len=seq_len,nnodes=tree$Nnode,ntaxa=ntaxa,tree=fin_tree)
	class(pgi.tree)<-"pgi.tree"
	return(pgi.tree)
}

