edit.cost <-
function(anc_dev_seq,der_dev_seq,output,def_penalty=2,method,seq_len) {
	#asses penalties for loss/gain
	idxs<-.clean_idx(anc_dev_seq,der_dev_seq,seq_len)
	missing_elements<-.identify_penalty(anc_dev_seq,der_dev_seq,seq_len)
	penalty<-length(missing_elements)*def_penalty
	#parse+clean sequences
	anc_dev_seq_t<-.clean_seq(anc_dev_seq,der_dev_seq,seq_len,1)
	der_dev_seq_t<-.clean_seq(anc_dev_seq,der_dev_seq,seq_len,2)
	# call the appropriate method
	if(method == "parsimov") {
		sel_evn<-.jeffery.parsimov(anc_dev_seq_t,der_dev_seq_t,idxs)
		if(output == 0) { return(length(sel_evn)+penalty) }
		if(output == 1) { return(length(sel_evn)) }
		if(output == 2) { return(sel_evn) }
		if(output == 3) { return(c(sel_evn,missing_elements)) }
	}
	if(method == "exhaustive-parsimov") {
		require(e1071)
		sel_evn<-.exhaustive.parsimov(anc_dev_seq_t,der_dev_seq_t)
		if(output == 0) { return(length(sel_evn)+penalty) }
		if(output == 1) { return(length(sel_evn)) }
		if(output == 2) { return(sel_evn) }
		if(output == 3) { return(c(sel_evn,missing_elements)) }
	}
	if(method == "event-pair") {
		score<-.event.pair.dist(anc_dev_seq_t,der_dev_seq_t)
		if(output == 0) { return(score+penalty) }
		if(output == 1) { return(score) }
	}
}

