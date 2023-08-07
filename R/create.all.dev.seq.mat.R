create.all.dev.seq.mat <-
function(seq.len,simul,verbosity=0) {
	if(seq.len > 8) stop("Sequence length too large for R for now. Quiting")
 require(e1071)
 if(verbosity > 1) cat("Creating Matrix of all Possible Sequences... ")
 if(simul) {
  all_dev_seq_mat<-.combi(seq(1,seq.len))
 } else {
  all_dev_seq_mat<-permutations(seq.len)
 }
 if(verbosity > 1) cat("Done.\n")
 return(all_dev_seq_mat)
}

