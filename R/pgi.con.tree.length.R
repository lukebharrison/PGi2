pgi.con.tree.length <-
function(con.tree) {
	t.l<-length(con.tree)
	this.node<-0
	if(length(con.tree[[1]]) == 1) {
		return(0)
	} else {
		for(i in 1:(length(con.tree)-1)) {
			this.node<-this.node+sum(con.tree[[1]][[2*i+2]][[1]])+sum(con.tree[[1]][[2*i+2]][[2]])+pgi.con.tree.length(con.tree[[i+1]])
		}
	}
	return(this.node)
}

