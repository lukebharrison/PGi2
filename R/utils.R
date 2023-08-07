## Define classes used in PGi
setClass("pgi.tree",slots=c(type="character",seq.len="numeric",nnodes="numeric",ntaxa="numeric",tree="list"))
setClass("multi.pgi.tree",slots="list")

## Now define methods for summary and plot for pgi.tree and multi.pgi.tree objects

setMethod("summary",signature("multi.pgi.tree"),summary.multi.pgi.tree)
setMethod("summary",signature("pgi.tree"),summary.pgi.tree)
setMethod("plot",signature(x="pgi.tree"),
	  function(x,y,...) {
		  plot.pgi.tree(pgi.tree=x,...)
	  }
)
setMethod("plot",signature(x="multi.pgi.tree"),
	  function(x,y,...) {
		  plot.multi.pgi.tree(multi.pgi.tree=x,...)
	  }
)
