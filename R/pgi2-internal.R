.c.pgi.tree<-function(...) {
	args<-list(...)
	class(args)<-"multi.pgi.tree"
	return(args)
}

.c.multi.pgi.tree<-function(...) {
	args<-list(...)
	new.multi<-list()
	for(i in 1:length(args)) {
		for(j in 1:length(args[[i]])) {
			new.multi<-c(new.multi,list(args[[i]][[j]]))
		}
	}
	class(new.multi)<-"multi.pgi.tree"
	return(new.multi)
}

.calc_trc_v <-
function(event_pair, order_v) { 
 trc_v <- rep(0,times=length(order_v))
 for(j in 1:length(order_v)) {
  if(order_v[j] > 1) {
   trc_v[j] <- sum(event_pair[(order_v[j]-1),])
  }
 }
 for(i in 1:length(order_v)) {
  if(order_v[i] < length(order_v)) {
   trc_v[i] <- trc_v[i] - sum(event_pair[,order_v[i]])
  }
 }
 trc_v_final<-rbind(trc_v,order_v)
 return(trc_v_final)
}

.calculate.consensus <-
function(filled_tree,seq_len,verbosity=0) {
	if(verbosity > 2) print("Calculating percentages on the completed consensus tree")
	num_child<-(length(filled_tree[[1]][[1]][1,]) - seq_len)-1
	num_anc<-filled_tree[[1]][[2]][[1]][seq_len+1]
	for(i in 2:(num_child+1)) {
		filled_tree[[1]][[2*i]][[1]]<-filled_tree[[1]][[2*i]][[1]]/filled_tree[[1]][[2*i-1]]
		filled_tree[[1]][[2*i]][[2]]<-filled_tree[[1]][[2*i]][[2]]/filled_tree[[1]][[2*i-1]]
		if(filled_tree[[i]][[1]][[1]][1,][seq_len+1] != -1) {
			filled_tree[[i]]<-.calculate.consensus(filled_tree[[i]],seq_len,verbosity=verbosity)
		}
	}
	for(i in 1:seq_len) {
		if(filled_tree[[1]][[2]][[2]][i] != 0) {
			if(filled_tree[[1]][[2]][[2]][i]/num_anc > 0.5) {
				filled_tree[[1]][[2]][[1]][i]<- -2
			}
		}
		if(filled_tree[[1]][[2]][[1]][i] != -2) {
			filled_tree[[1]][[2]][[1]][i]<-filled_tree[[1]][[2]][[1]][i]/(num_anc-filled_tree[[1]][[2]][[2]][i])
		}
	}
	filled_tree[[1]][[2]][[1]]<-round(filled_tree[[1]][[2]][[1]])
	filled_tree[[1]][[2]][[1]][1:seq_len]<-compress.rds(filled_tree[[1]][[2]][[1]][1:seq_len])
	return(filled_tree)
}

.clean.edge.mat <-
function(edges,anc,numc) {
 new.edges<-c()
 count<-0
 for(i in 1:dim(edges)[1]) {
  if(edges[i,1] == anc) {
   count<-count+1
  }
  if(count != numc) {
   new.edges<-rbind(new.edges,edges[i,])
  }
 }
 dim(new.edges)<-c(length(new.edges)/2,2)
 return(new.edges)
}

.clean_idx <-
function(anc_dev_seq, der_dev_seq, seq_len) {
 true_idx<-c()
 for(i in 1:seq_len) {
  if(anc_dev_seq[i] == -2 || der_dev_seq[i] == -2) {
   true_idx<-c(true_idx,i)
  }
 }
 return(true_idx)
}

.clean_seq <-
function(anc_dev_seq,der_dev_seq,seq_len,idx) {
 dims<-c(1,seq_len)
 i<-1
 while(i <= dims[2]) {
  if(anc_dev_seq[i] == -2) {
   if(der_dev_seq[i] == -2) {
    anc_dev_seq<-anc_dev_seq[-i]
    der_dev_seq<-der_dev_seq[-i]
    dims[2]<-dims[2]-1
   } else {
    anc_dev_seq<-anc_dev_seq[-i]
    der_dev_seq<-.cut_rds(der_dev_seq,i,dims[2])
    dims[2]<-dims[2]-1
   }
   next
  }
  if(der_dev_seq[i] == -2 && anc_dev_seq[i] != -2) {
   der_dev_seq<-der_dev_seq[-i]
   anc_dev_seq<-.cut_rds(anc_dev_seq,i,dims[2])
   dims[2]<-dims[2]-1
   next
  }
  i<-i+1
 }
 dim(anc_dev_seq)<-dims
 dim(der_dev_seq)<-dims
 if(idx == 1) {
  return(anc_dev_seq)
 }
 if(idx == 2) {
  return(der_dev_seq)
 }
}

.combi <-
function(x) {
    n<-length(x)
    nn<-n^n
    cm<-matrix(-1,nn,n)
    for (i in 1:n) cm[,i]<-rep(x,each=n^(i-1))
    return (cm)
}

.create.event.pair.dif.mat <-
function(v1,v2) {
  event_pair1 <- array(0, dim=c(length(v1)-1,length(v1)-1))
  for (i in 2:length(v1)-1) { 
   for (j in (i+1):length(v1)) {
    if(v1[j] > v1[i]) event_pair1[j-1,i] <- 2 
    if(v1[j] < v1[i]) event_pair1[j-1,i] <- 0
    if(v1[j] == v1[i]) event_pair1[j-1,i] <- 1
   }
   j <- 0
  }
  event_pair2 <- array(0, dim=c(length(v2)-1,length(v2)-1))
  for (i in 2:length(v2)-1) {
   for (j in (i+1):length(v2)) {
    if(v2[j] > v2[i]) event_pair2[j-1,i] <- 2 
    if(v2[j] < v2[i]) event_pair2[j-1,i] <- 0
    if(v2[j] == v2[i]) event_pair2[j-1,i] <- 1
   }
  j <- 0
  }
 event_pair <- event_pair2 - event_pair1
 return(event_pair)
}

.cut_rds <-
function(rds_dev_seq,idx,seq_len) {
 i<-1
 temp_seq<-rds.to.seq(rds_dev_seq)
 ln<-length(temp_seq)
 while(i <= ln) {
  j<-1
  ln1<-length(temp_seq[[i]])
  while(j <= ln1) {
   if(temp_seq[[i]][j] == idx) {
    if(length(temp_seq[[i]]) < 2) {
     temp_seq[[i]]<-NULL
     ln<-ln-1
    } else {
     temp_seq[[i]]<-temp_seq[[i]][-j]
     ln1<-ln1-1
    }
   }
   j<-j+1
  }
  i<-i+1
 }
 rds_dev_seq<-seq.to.rds(temp_seq,seq_len)
 rds_dev_seq<-rds_dev_seq[-idx]
 return(rds_dev_seq)
}

.downsweep <-
function(partial.pgi.tree,der_dev_seqs,progress=NULL,verbosity=0) {
 # Case where internal node; recurse to a node that is calculable
	if(der_dev_seqs[[1]][1] == -1) {
		for(i in 2:length(der_dev_seqs)) {
			if(length(der_dev_seqs[[i]][[1]]) == 1) {
 				if(der_dev_seqs[[i]][[1]] == -1 || der_dev_seqs[[i]][[1]] == 0) {
					if(verbosity > 2) print("Recursing down the tree to a calculable node")
					t.ret<-.downsweep(partial.pgi.tree,der_dev_seqs[[i]],progress=progress,verbosity=verbosity)
					der_dev_seqs[[i]]<-t.ret[[2]]
					partial.pgi.tree<-t.ret[[1]]
				}
			}
		}
	}
	if(verbosity > 2) print("Found a calculable node, begining inference on this node")
	anc_dev_seq<-seq(1:partial.pgi.tree$seq.len)
	if(any(partial.pgi.tree$inf.params$heuristic == c("pgi","oldgenetic"))) {
		if(partial.pgi.tree$inf.params$edit.cost.func != "event-pair") {
			for(i in 1:partial.pgi.tree$seq.len) {
				prob<-0
				for(j in 2:length(der_dev_seqs)) {
					if(der_dev_seqs[[j]][[1]][i] == -2) {
						prob<-prob+1
					}
				}
				if(prob > 0) {
					#print("we have unscorable data here")
					prob<-prob/(length(der_dev_seqs)-1)
					if(runif(1) < prob) {
						# add a -2 @ this position
						dev_seq<-.cut_rds(anc_dev_seq,i,partial.pgi.tree$seq.len)
					} else {
						anc_dev_seq<-.gen_ran2(anc_dev_seq,partial.pgi.tree$inf.params$simultaneity)
					}
				}
			}
		}
		# a bit of code to seed a better anc_dev_seq for unscorable data.
		# basically if all dll derived seqs have teh same unscorable then prob so does anc
		# otherwise probability is the proportion that do
		seq_tosend <- c()
		for(i in 2:length(der_dev_seqs)) {
			seq_tosend<-c(seq_tosend,der_dev_seqs[[i]][1])
		}
		anc_dev_seq_mat<-.refine.anc.seqs(anc_dev_seq,seq_tosend,partial.pgi.tree$inf.params,0.15,verbosity=verbosity)

		# if the sampling heurisitc is applied: reduce the size of the ancestral seq matrix to sample_matrix
		if(partial.pgi.tree$inf.params$ret.anc.seq > 0 && (dim(anc_dev_seq_mat)[1] > 1)) {
			if(partial.pgi.tree$inf.params$ret.anc.seq > (dim(anc_dev_seq_mat)[1])) {
				smatrix_temp <- (dim(anc_dev_seq_mat)[1])
			} else {
				smatrix_temp <- partial.pgi.tree$inf.params$ret.anc.seq
			}
			anc_dev_seq_mat<-anc_dev_seq_mat[sample(seq(1:(dim(anc_dev_seq_mat)[1])),smatrix_temp,replace = FALSE),]
			anc_dev_seq_mat<-anc_dev_seq_mat[order(anc_dev_seq_mat[,partial.pgi.tree$seq.len+1],decreasing=TRUE),]
		}
		anc_dev_seq<-anc_dev_seq_mat[which.min(anc_dev_seq_mat[,(partial.pgi.tree$seq.len+1)]),]
		der_dev_seqs[[1]]<-anc_dev_seq_mat
		if(verbosity == 1) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
		if(partial.pgi.tree$inf.params$resume) {
			partial.pgi.tree$partial.tree<-.pgi.fill.partial.tree(partial.pgi.tree$partial.tree,der_dev_seqs)
			save(partial.pgi.tree,file=partial.pgi.tree$inf.params$resume.temp.file)
		}
		return(list(partial.pgi.tree,der_dev_seqs))
	} 

	# Exhaustive (non-heursitic) method
	if(partial.pgi.tree$inf.params$heuristic == "exhaustive") {
		seq_tosend <- c()
		for(i in 2:length(der_dev_seqs)) {
			seq_tosend<-c(seq_tosend,der_dev_seqs[[i]][1])
		}
		anc_dev_seq_mat<-.exhaustive.refine.anc.seqs(seq_tosend,partial.pgi.tree$inf.params,partial.pgi.tree$seq.len,verbosity=verbosity)
		der_dev_seqs[[1]]<-anc_dev_seq_mat
		if(verbosity == 1) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
		if(partial.pgi.tree$inf.params$resume) {
			partial.pgi.tree$partial.tree<-.pgi.fill.partial.tree(partial.pgi.tree$partial.tree,der_dev_seqs)			
			save(partial.pgi.tree,file=partial.pgi.tree$inf.params$resume.temp.file)
		}
		return(list(partial.pgi.tree,der_dev_seqs))
	}

	# SW04 search-based optimization
	if(partial.pgi.tree$inf.params$heuristic == "sw04") {
		# recode replicates as the edit cost matrix, cycles as the permu matrix	
		seq_tosend <- c()
		for(i in 2:length(der_dev_seqs)) {
			seq_tosend<-c(seq_tosend,der_dev_seqs[[i]][1])
		}
		#print("Hypothesized Ancestral Seq:")
		anc_dev_seq_mat<-.exhaustive.refine.anc.seqs(seq_tosend,partial.pgi.tree$inf.params,partial.pgi.tree$seq.len,verbosity=verbosity)
		der_dev_seqs[[1]]<-anc_dev_seq_mat
		if(verbosity == 1) setTxtProgressBar(progress, getTxtProgressBar(progress)+1)
		if(partial.pgi.tree$inf.params$resume) {
			partial.pgi.tree$partial.tree<-.pgi.fill.partial.tree(partial.pgi.tree$partial.tree,der_dev_seqs)
			save(partial.pgi.tree,file=partial.pgi.tree$inf.params$resume.temp.file)
		}
		return(list(partial.pgi.tree,der_dev_seqs))
	}
}

.event.pair.dist <-
function(anc_dev_seq,der_dev_seq) {
 difs<-.create.event.pair.dif.mat(anc_dev_seq,der_dev_seq) 
 score<-sum(abs(difs))
 return(score)
}

.exhaustive.parsimov <-
function(anc_dev_seq, der_dev_seq) {
 all_dev_seq_mat<-permutations(length(anc_dev_seq))
 score_matrix<-rep(1000,times=factorial(length(anc_dev_seq)))
 event_pair<-.create.event.pair.dif.mat(anc_dev_seq,der_dev_seq)
 trc_v<-.calc_trc_v(event_pair,seq(1,length(anc_dev_seq)))
 trc_v<-trc_v[1,]
 for(i in 1:length(score_matrix)) {
  if(i %% 500 == 0) print(i/length(score_matrix)*100)
  score_matrix[i] <- length(parsimov.soln(.calc_trc_v(event_pair,all_dev_seq_mat[i,]),event_pair))
 }
 return(.parsimov.soln(.calc_trc_v(event_pair,all_dev_seq_mat[which.min(score_matrix),]),event_pair))
}


.exhaustive.refine.anc.seqs <-
function(der_dev_seqs,inf.params,seq_len,verbosity=0) {
	num_scores<-dim(inf.params$all.dev.seq.mat)[1]
	score_matrix<-rep(1000,times=num_scores)
	for(i in 1:length(der_dev_seqs)) {
		score_matrix<-cbind(score_matrix,rep(1000,times=num_scores))
	}
	if(verbosity == 2) pb.node<-txtProgressBar(min=1,max=num_scores,style=3)
 for(j in 1:num_scores) {
   temp_score<-0
   temp_best<-0
   for(z in 1:(length(der_dev_seqs))) {
    temp_score2<-0
    if(der_dev_seqs[[z]][1,seq_len+1] == -1) {
     q<-which(duplicated(rbind(der_dev_seqs[[z]][1,][1:seq_len],inf.params$all.dev.seq.mat)))-1
     offset <- 1
    } else {
     offset <- 0
     q<-1
    }
    #assess unscorable data penalties
    #2 = default penalty (can be changed)
    # assessing US Penalties
    if(inf.params$heuristic == "exhaustive") {
     temp_best<-(edit.cost(inf.params$all.dev.seq.mat[j,],inf.params$all.dev.seq.mat[q,],0,2,inf.params$edit.cost.func,seq_len))+der_dev_seqs[[z]][1,seq_len+1] + offset
    } else {
     # SW04
     temp_best<-(inf.params$edit.cost.matrix[j,q])+der_dev_seqs[[z]][1,seq_len+1] + offset
    }
    score_matrix[j,z+1]<-1
    # For number of derived sequences at the child node
    for(k in 1:dim(der_dev_seqs[[z]])[1]) {
     if(k == 1) {
      next
     }
     #For each seq @ each child node
     offset <- 0
     if(inf.params$heuristic == "exhaustive") {
      #Asesss the penalties for unscorable data + add to score
      temp_score2<-(edit.cost(inf.params$all.dev.seq.mat[j,],inf.params$all.dev.seq.mat[k,],0,2,inf.params$edit.cost.func,seq_len))+der_dev_seqs[[z]][k,seq_len+1] + offset
     } else {
      #SW04
      temp_score2<-(inf.params$edit.cost.matrix[j,k])+der_dev_seqs[[z]][k,seq_len+1] + offset
     }
     if(temp_score2 < temp_best) {
      temp_best<-temp_score2
      score_matrix[j,z+1]<-k
     }
    } # End mulptiple derived seqs

    # Add to the total from the other nodes
    temp_score<-temp_score+temp_best

   } # End the # of nodes
   score_matrix[j,1] <- temp_score
   if(verbosity > 2) cat("Sequence #",j,"of",num_scores,"- Current best score:",min(score_matrix[,1]),"\n")
		if(verbosity == 2) setTxtProgressBar(pb.node,j)
  }
		if(verbosity == 2) close(pb.node)
 min_idx<-which.min(score_matrix[,1])
 min<-score_matrix[min_idx,1]
 aug_dev_seq_mat<-cbind(inf.params$all.dev.seq.mat[seq(1:num_scores),],score_matrix)
 dim(aug_dev_seq_mat)<-c(num_scores,seq_len+1+length(der_dev_seqs))
 return(aug_dev_seq_mat)
}

.gen_ran <-
function(rds_dev_seq,replicates,simul) {
	dev_seq_mat<-.gen_ran2(rds_dev_seq,simul)
	for(i in 1:(replicates-1)) {
		temp_seq <- .gen_ran2(rds_dev_seq,simul)
		j<-0
		while(dim(unique(rbind(temp_seq,rds_dev_seq)))[1] < 2 && j < replicates/2) { 
			j<-j+1
			if((!simul) && (j > factorial(length(rds_dev_seq)/4))) break
			temp_seq<-.gen_ran2(rds_dev_seq,simul)
		} 
		dev_seq_mat<-rbind(dev_seq_mat,temp_seq)
	}
	return(dev_seq_mat)
}

.gen_ran2 <-
function(rds_dev_seq,simul) {
 seq_len<-length(rds_dev_seq)
 us<-0
 dev_seq<-rds.to.seq(rds_dev_seq)
 if(simul) {
  collision_chance <- 0.5
 } else {
  collision_chance <- 0
 }
 # Random index to shift
 ran_idx <- 0
 while(ran_idx < 1) {
  ran_idx<-round(runif(1)*length(rds_dev_seq))
 }

 if(rds_dev_seq[ran_idx] == -2) {
  us<-1
  if(runif(1) > 0.10) {
   # Cop out early - no need to move it with 90% chance
   return(rds_dev_seq)
  }
 } else {
  # Case were we have a 10% chance we create an unscorable data point where one did not exist 
  if((runif(1) < 0.01)) {# CHANGE TO DISABLE runif(1) < 0.02
   ln<-dev_seq[[rds_dev_seq[ran_idx]]]
   for(i in 1:length(ln)) {
    if(ln[i] == ran_idx) {
     if(length(ln) > 1) {
      ln<-ln[-i]
     }
     else {
      ln <- NULL
     }
     dev_seq[[rds_dev_seq[ran_idx]]]<-ln
     return(seq.to.rds(dev_seq,seq_len))
    }
   }
  }
 }
 ran_dist <- runif(1) * max((max(rds_dev_seq)-rds_dev_seq[ran_idx]),abs(1-rds_dev_seq[ran_idx]))
 # randomly choose a direction of shift
 if(runif(1) > 0.5) ran_dist<-(-1 * ran_dist)
 ran_dist <- round(ran_dist)
 # Here - if taget == -2, then we have a 10% chance of losing it
 # Which means that we
 # here we can also copute teh chance that a non -2 target will be converted to
 # unscorable at a .10 cahnce.
 if((ran_dist + rds_dev_seq[ran_idx]) > max(rds_dev_seq)) {
  ran_dist <- max(rds_dev_seq) - rds_dev_seq[ran_idx]
 }
 if((ran_dist + rds_dev_seq[ran_idx]) < 1) {
  ran_dist <- 1 - rds_dev_seq[ran_idx]
 }
 target<-rds_dev_seq[ran_idx]+ran_dist 
 if(us==1) {
  ran_dist<-round(runif(1)*max(rds_dev_seq))
  while(ran_dist<1) {
   ran_dist<-round(runif(1)*max(rds_dev_seq))
  }
  target<-ran_dist
 }
 if(target == rds_dev_seq[ran_idx]) return(rds_dev_seq)
 # Remove the original seq
 if(us!=1) {
  if(length(dev_seq[[rds_dev_seq[ran_idx]]]) > 1) {
   dev_seq[[rds_dev_seq[ran_idx]]][dev_seq[[rds_dev_seq[ran_idx]]] == ran_idx] <- -1
  } else {
   dev_seq[[rds_dev_seq[ran_idx]]] <- -1
  }
 }
 if((runif(1) < collision_chance) && simul) {
  # chain the selcted value to its new position
  dev_seq[[target]] <- c(dev_seq[[target]],ran_idx)
 } else {
  # Save the displaced values
  displaced <- dev_seq[[target]]
  # put in the new element
  dev_seq[[target]] <- ran_idx
  # 50% chance we insert before or after
  dev_seq <- c(dev_seq,c(-1))
  if(runif(1) < 0.5) {
   for(i in (target+1):length(dev_seq)) {
    displaced2<-dev_seq[[i]]
    dev_seq[[i]]<-displaced
    displaced<-displaced2
   }
  } else {
   for(i in target:length(dev_seq)) {
    displaced2<-dev_seq[[i]]
    dev_seq[[i]]<-displaced
    displaced<-displaced2
   }
  }
 }
 rds_dev_seq<-seq.to.rds(dev_seq,seq_len)
 return(rds_dev_seq)
}

.identify_penalty <-
function(anc_dev_seq,der_dev_seq,seq_len) {
 elements<-c()
 for(i in 1:seq_len) {
  if(anc_dev_seq[i] == -2) {
   if(der_dev_seq[i] != -2) {
    elements<-c(elements,list(i+0.1))
   }
  } else {
   if(der_dev_seq[i] == -2) {
    elements<-c(elements,(i+0.1)*-1)
   }
  }
 }
 return(elements)
}

.jeffery.parsimov <-
function(anc_dev_seq, der_dev_seq,idxs) {
	if(missing(idxs)) idxs<-.clean_idx(anc_dev_seq,der_dev_seq,length(anc_dev_seq))
 sel_evn<-c()
 event_pair<-.create.event.pair.dif.mat(anc_dev_seq, der_dev_seq)
 trc_v<-.calc_trc_v(event_pair,(seq(1,length(anc_dev_seq))))
 while(sum(abs(trc_v[1,] != 0))) {
  i<-which.max(abs(trc_v[1,]))
  if(trc_v[1,i] > 0) {
   sel_evn<-c(sel_evn,-trc_v[2,][i])
  } else {
   sel_evn<-c(sel_evn,trc_v[2,][i])
  }
  if(trc_v[2,i] < length(trc_v[2,])) {
   event_pair[,trc_v[2,i]] <- 0
  }
  if(trc_v[2,i] >= 2) {
   event_pair[(trc_v[2,i]-1),] <- 0
  }
  trc_v<-.calc_trc_v(event_pair, trc_v[2,])
 }
 if(length(sel_evn) < 1) {
  return(NULL)
 }
 for(i in 1:length(idxs)) {
  temp_sel_evn<-sel_evn[abs(sel_evn) >= idxs[i]] 
  if(length(temp_sel_evn) > 0) {
   for(j in 1:length(temp_sel_evn)) {
    if(temp_sel_evn[j] > 0) {
     temp_sel_evn[j] <- temp_sel_evn[j] + 1
    } else {
     temp_sel_evn[j] <- temp_sel_evn[j] - 1
    }
   }
  }
  sel_evn[abs(sel_evn) >= idxs[i]] <- temp_sel_evn
 }
 return(sel_evn)
}


.parsimov.soln <-
function(trc_v_aug, event_pair) {
 sel_evn<-c()
 # Need to change for() to look at all possible input orders.... need heuristic
 for(i in 1:length(trc_v_aug[1,])) {
  if(sum(abs(trc_v_aug[1,])) == 0) {
   break
  }
  #now we parse each event in order and check its changes
  if(trc_v_aug[1,i] != 0) {
   #ignore events involved in no shifts at all - usually terminal
   #Now we assume event i 'occurs' ie all changes w.r.to a are eliminated
   if(trc_v_aug[1,i] > 0) {
    sel_evn<-c(sel_evn,-trc_v_aug[2,i])
   } else {
   sel_evn<-c(sel_evn,trc_v_aug[2,i])
   }
   if(trc_v_aug[2,i] < length(trc_v_aug[2,])) {
    event_pair[,trc_v_aug[2,i]] <- 0
   }
   if(trc_v_aug[2,i] >= 2) {
   event_pair[(trc_v_aug[2,i]-1),] <- 0
   }
   trc_v_aug<-.calc_trc_v(event_pair, trc_v_aug[2,])
  }
 }
 return(sel_evn)
}

.pgi.build.tree <-
function(ape_tree,edges,data,ntaxa,seq_len) {
	root_edge<-edges[edges[,1] == min(edges[,1])]
	dim(root_edge)<-c(length(root_edge)/2,2)
	num_child<-dim(root_edge)[1]
	if(length(root_edge[root_edge[,2] <= ntaxa]) == length(root_edge)) {
		tree<-list(0)
		for(i in 1:num_child) { 
			#tree<-c(tree,list(list(.prep_seq(data[[ape_tree[[2]][[root_edge[i,2]]]]],ape_tree[[2]][[root_edge[i,2]]]))))
			tip.idx<-which(names(data) %in% ape_tree$tip.label[root_edge[i,2]])
			tree<-c(tree,list(list(.prep_seq(data[[tip.idx]],ape_tree$tip.label[root_edge[i,2]]))))
		}
	} else {
		tree<-list(-1)
		for(i in 1:num_child) {
			if(root_edge[i,2] <= ntaxa) {
				#tree<-c(tree, list( list(.prep_seq( data[[ ape_tree[[2]][[root_edge[i,2]]]  ]], ape_tree[[2]][[root_edge[i,2]]]))))
				tip.idx<-which(names(data) %in% ape_tree$tip.label[root_edge[i,2]])
				tree<-c(tree,list(list(.prep_seq(data[[tip.idx]],ape_tree$tip.label[root_edge[i,2]]))))
			} else {
				#node
				edges_temp<-edges[edges[,1] != root_edge[i,1]]
				dim(edges_temp)<-c(length(edges_temp)/2,2)
				tree<-c(tree,list(.pgi.build.tree(ape_tree,edges_temp,data,ntaxa,seq_len)))
				edges<-.clean.edge.mat(edges,root_edge[i,1],i) 
			}
		}
	}
	return(tree)
}

.pgi.check.nnodes <-
function (tree) {
	if(length(tree[[1]]) > 1) {
		# not an empty node either a tip or a filled node
		if(dim(tree[[1]])[1] > 1) {
			# filled node
			nnodes<-1
		} else {
			nnodes<-0
		}
	} else {
		nnodes<-0
	}
	for(i in 2:length(tree)) {
		if(length(tree[[i]][[1]]) > 1) {
			if(dim(tree[[i]][[1]])[1] == 1) {
				next
			} else {
				nnodes<-nnodes+.pgi.check.nnodes(tree[[i]])	
			}
		} else {
			if(tree[[i]][[1]] == 0) {
				#nothing todo
			} else if(tree[[i]][[1]] == -1) {
				nnodes<-nnodes+.pgi.check.nnodes(tree[[i]]) 
			}
		}
	}
	return(nnodes)
}

.pgi.con.to.phylo.fill <-
function(tree,st,seq_len) {
  st[[6]]<-st[[6]]+1
  st[[4]]<-st[[4]]+1
  st[[8]]<-c(st[[8]],print.rds(tree[[1]][[2]][[1]][1:seq_len]))
  node<-st[[6]]
  for(i in 2:length(tree)) {
   if(tree[[i]][[1]][[1]][1,seq_len+1] != -1) {
    #node
    edge1<-c(node,st[[6]]+1)
    accel<-c()
    decel<-c()
		for(z in 1:seq_len) {
   	  if(!is.na(tree[[1]][[2*i]][[1]][z]) && tree[[1]][[2*i]][[1]][z] > 0.5) {
    	  accel<-c(accel,z)
			}
     if(!is.na(tree[[1]][[2*i]][[2]][z]) && tree[[1]][[2*i]][[2]][z] > 0.5) {
      decel<-c(decel,z)
     }
    }
		st[[10]]<-c(st[[10]],rbind(tree[[1]][[2*i]][[1]]))
		st[[11]]<-c(st[[11]],rbind(tree[[1]][[2*i]][[2]]))
    elen<-length(accel)+length(decel)
    if(elen == 0) {
     elen <- 1
    }
    ch_str<-c()
    if(length(accel) > 0) {
			ch_str<-paste(ch_str,"a:",toString(accel))
		} 
    if(length(decel) > 0) {
     ch_str<-paste(ch_str,"d:",toString(decel))
		}
    if(length(ch_str) == 0) {
     ch_str<-" "
    }
    st[[9]]<-c(st[[9]],ch_str)
    st[[7]]<-c(st[[7]],elen)
    st[[1]]<-rbind(st[[1]],edge1,deparse.level=0)
    st<-.pgi.con.to.phylo.fill(tree[[i]],st,seq_len)
   } else {
    # extant seq
    #cat("compting internal#",node,"to extant taxa #:",st[[4]],"\n")
    edge1<-c(node,st[[5]])
		st[[10]]<-c(st[[10]],rbind(tree[[1]][[2*i]][[1]]))
		st[[11]]<-c(st[[11]],rbind(tree[[1]][[2*i]][[2]]))
    acc<-tree[[1]][[2*i]][[1]][tree[[1]][[2*i]][[1]] > 0.5]
    del<-tree[[1]][[2*i]][[1]][tree[[1]][[2*i]][[2]] > 0.5]
    elen<-length(acc[!is.na(acc)])+length(del[!is.na(del)])
    accel<-c()
    decel<-c()
    for(z in 1:seq_len) {
     if(!is.na(tree[[1]][[2*i]][[1]][z]) && tree[[1]][[2*i]][[1]][z] > 0.5) {
      accel<-c(accel,z)
     }
     if(!is.na(tree[[1]][[2*i]][[2]][z]) && tree[[1]][[2*i]][[2]][z] > 0.5) {
      decel<-c(decel,z)
     }
    }
     ch_str<-c()
    if(length(accel) > 0) {
			ch_str<-paste(ch_str,"a:",toString(accel))
		}
    if(length(decel) > 0) {
     ch_str<-paste(ch_str,"d:",toString(decel))
		}    
		if(length(ch_str) == 0) {
     ch_str<-" " 
    }
    st[[9]]<-c(st[[9]],ch_str)
    st[[7]]<-c(st[[7]],elen)
    #add edge len
    st[[5]]<-st[[5]]+1
		st[[2]]<-c(st[[2]],rownames(tree[[i]][[1]][[1]]))
    st[[3]]<-c(st[[3]],print.rds(tree[[i]][[1]][[1]][1,]))
    st[[1]]<-rbind(st[[1]],edge1,deparse.level=0)
   }
 }
 return(st)
}

.pgi.exhaustive.consensus <-
function(pgi.tree,verbosity=0) {
	seq_len<-pgi.tree$seq.len
	con_tree<-.prepare.con.tree.struct(pgi.tree$filled.tree,pgi.tree$seq.len,verbosity=verbosity)
	if(verbosity > 2) print("Prepared the empty pseudo consensus tree data structure")
	min_score<-min(con_tree[[1]][[1]][,seq_len+1])
	num.min.sequences<-sum(as.integer(con_tree[[1]][[1]][,seq_len+1] == min_score))
	if(verbosity > 0 && pgi.tree$con.params$semi.ex.con.max.n > 0) cat("\nSemi Exhaustive consensus: maximum number of total per branch transversals =",pgi.tree$con.params$semi.ex.con.max.n,"\n")
	if(pgi.tree$con.params$semi.ex.con.max.n > 0) per.run.sample.max<-round(pgi.tree$con.params$semi.ex.con.max.n/num.min.sequences) else per.run.sample.max<-0
	if(verbosity > 0) {
		pb<- txtProgressBar(min = 0, max = num.min.sequences, style = 3)
	}
	if(verbosity > 1 && pgi.tree$con.params$semi.ex.con.max.n > 0) cat("\nSemi exhaustive consensus: # of per branch transversals allowed per equally parsimonious root sequence:",per.run.sample.max,"\n")
	for(i in 1:length(con_tree[[1]][[1]][,seq_len+1])) {
		if(con_tree[[1]][[1]][,seq_len+1][i] == min_score) {
			con_tree<-.upsweep(con_tree,i,pgi.tree$con.params$edit.cost.func,pgi.tree$seq.len,sample.max=per.run.sample.max,verbosity=verbosity)
			if(verbosity > 0) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
		}
	}
	if(verbosity > 2)	print(con_tree)
	#return(con_tree)
	pgi.tree$consensus<-.calculate.consensus(con_tree,seq_len,verbosity=verbosity)
	pgi.tree$type="consensus"
	if(verbosity > 0) close(pb)
	return(pgi.tree)
}

# function to match a partial pgi tree to the complete empty tree
.pgi.fill.partial.tree<-function(empty.tree,partial.tree) {
	if(.pgi.tree.hash(empty.tree) == .pgi.tree.hash(partial.tree)) {
		empty.tree<-partial.tree
	}
	for(i in 2:length(empty.tree)) {
		if(length(empty.tree[[i]][[1]]) == 1) empty.tree[[i]]<-.pgi.fill.partial.tree(empty.tree[[i]],partial.tree)
	}
	return(empty.tree)
}

.pgi.getn.completed.nodes <-
function(pgi.tree) {
	if(pgi.tree$type == "empty") {
		return(0)
	} else if(any(c("filled","consensus") %in% pgi.tree$type == "filled")) {
	 return(pgi.tree$nnodes)
	} else if(pgi.tree$type == "partial") {
		return(.pgi.check.nnodes(pgi.tree$partial.tree))
	} else {
		return(NA)
	}
}


.pgi.old.get.seq.len <-
function(tree,verbosity=0) {
	if(is.list(tree[[1]])) {
		# con tree, note does not work on a filled, non-consensus tree
		if(verbosity > 2) print("Getting PGi tree sequence length, consensus tree")
		return(length(tree[[1]][[4]][[1]]))
	} else {
		t.l<-length(tree[[1]])
		# non-con tree
		if(t.l == 1) {
			if(verbosity > 2) print("Getting PGi tree sequence length, un-filled tree")
			l.l<-length(tree[[2]][[1]])
			if(l.l != 1) return(l.l-2)
			r.l<-length(tree[[3]][[1]])
			if(r.l != 1) return(r.l-2)
			return(pgi.old.get.seq.len(tree[[2]]))
		} else {
			if(verbosity > 2) print("Getting PGi tree sequence length, filled, non-consensus tree")
			if(length(tree) == 1) return(length(tree[[1]])-2) else .pgi.old.get.seq.len(tree[[2]])
		}
	}
}


.pgi.old.nnodes <-
function (pgi_tree) {
	nnodes<-1
	if(length(pgi_tree[[1]]) >= 6) {
		for(i in 2:length(pgi_tree)) {			
			if(pgi_tree[[i]][[1]][[1]][1,.pgi.old.get.seq.len(pgi_tree)+1] == -1) {
				next 
			} else {
				nnodes<-nnodes+.pgi.old.nnodes(pgi_tree[[i]])
			}
		}
	} else {
		for(i in 2:length(pgi_tree)) {
			if(length(pgi_tree[[i]][[1]]) > 1) next
			if(pgi_tree[[i]][[1]] == 0) {
				nnodes<-nnodes+1
			} else if(pgi_tree[[i]][[1]] == -1) {
				nnodes<-nnodes+.pgi.old.nnodes(pgi_tree[[i]]) 
			}
		}
	}
	return(nnodes)
}

.pgi.old.ntaxa <-
function (pgi_tree) {
	ntaxa<-0
	if(length(pgi_tree[[1]]) >= 6) {
		for(i in 2:length(pgi_tree)) {
			if(dim(pgi_tree[[i]][[1]][[1]])[1] == 1) {
				ntaxa<-ntaxa+1
			} else {
				ntaxa<-ntaxa+.pgi.old.ntaxa(pgi_tree[[i]])
			}
		}
	} else {
		for(i in 2:length(pgi_tree)) {
			if(length(pgi_tree[[i]][[1]]) > 1) {
				ntaxa<-ntaxa+1
			} else if(pgi_tree[[i]][[1]] == 0) {
				ntaxa<-ntaxa+2
			} else if(pgi_tree[[i]][[1]] == -1) {
				ntaxa<-ntaxa+.pgi.old.ntaxa(pgi_tree[[i]]) 
			}
		}
	}
	return(ntaxa)
}

.pgi.simple.con <-
function (pgi.tree, verbosity=0) {
	con_tree<-.prepare.con.tree.struct(pgi.tree$filled.tree,pgi.tree$seq.len)
	min_score<-min(con_tree[[1]][[1]][,pgi.tree$seq.len+1])
	num.min.scores<-sum(as.integer(con_tree[[1]][[1]][,pgi.tree$seq.len+1] == min_score))
	if(verbosity > 0) pb<-txtProgressBar(min=0,max=num.min.scores,style=3)
	for(i in 1:length(con_tree[[1]][[1]][,pgi.tree$seq.len+1])) {
		if(con_tree[[1]][[1]][,pgi.tree$seq.len+1][i] == min_score) {
			con_tree<-.upsweep.simple(con_tree,i,pgi.tree$con.params$edit.cost.func,pgi.tree$seq.len)
		}
		if(verbosity > 0) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
	}
	if(verbosity > 0) close(pb)
	if(verbosity > 2) print(con_tree)
	pgi.tree$consensus<-.calculate.consensus(con_tree,pgi.tree$seq.len)
	pgi.tree$type<-"consensus"
	return(pgi.tree)
}


.pgi.supercon.fill <-
function(supercontree,con_trees,num_trees,seq_len,verbosity=0) {
	num_child<-(length(con_trees[[1]][[1]][[1]][1,]) - seq_len)-1
	if(verbosity > 2)	cat("recovered # children:",num_child,"\n")
	supercontree[[1]][[2]][[1]]<-rep(0,times=seq_len+1)
	supercontree[[1]][[2]][[2]]<-rep(0,times=seq_len+1)
	for(i in 1:num_child) {
		supercontree[[1]][[2*i+2]][[1]]<-rep(0,times=seq_len)
		supercontree[[1]][[2*i+2]][[2]]<-rep(0,times=seq_len)
		dim(supercontree[[1]][[2*i+2]][[1]])<-c(1,seq_len)
		dim(supercontree[[1]][[2*i+2]][[2]])<-c(1,seq_len)
	}
	if(verbosity > 2) cat("Completed initialization of the superconsensus data matrix\n")	
 for(i in 1:num_trees){
  for(j in 1:num_child) { 
   #sequence
   supercontree[[1]][[2]][[1]]<-supercontree[[1]][[2]][[1]]+con_trees[[i]][[1]][[2]][[1]]
   supercontree[[1]][[2]][[2]]<-supercontree[[1]][[2]][[2]]+con_trees[[i]][[1]][[2]][[2]]
   #accels+deccels
   supercontree[[1]][[2*j+2]][[1]]<-supercontree[[1]][[2*j+2]][[1]]+con_trees[[i]][[1]][[2*j+2]][[1]]
   supercontree[[1]][[2*j+2]][[2]]<-supercontree[[1]][[2*j+2]][[2]]+con_trees[[i]][[1]][[2*j+2]][[2]]
  }
 }
 supercontree[[1]][[2]][[1]]<-supercontree[[1]][[2]][[1]]/num_trees
 supercontree[[1]][[2]][[2]]<-supercontree[[1]][[2]][[2]]/num_trees
 for(j in 1:seq_len) {
  if(supercontree[[1]][[2]][[2]][j] > 0.5) {
   supercontree[[1]][[2]][[1]][j]<- -2
  }
 }
 supercontree[[1]][[2]][[1]]<-round(supercontree[[1]][[2]][[1]])
 supercontree[[1]][[2]][[1]][1:seq_len]<-compress.rds(supercontree[[1]][[2]][[1]][1:seq_len])
 for(j in 1:num_child) {
  #accels+deccels
  supercontree[[1]][[2*j+2]][[1]]<-supercontree[[1]][[2*j+2]][[1]]/num_trees
  supercontree[[1]][[2*j+2]][[2]]<-supercontree[[1]][[2*j+2]][[2]]/num_trees
  if(supercontree[[j+1]][[1]][[1]][1,][seq_len+1] != -1) {
   con_trees_tosend<-c()
   for(i in 1:num_trees) {
    con_trees_tosend<-c(con_trees_tosend,list(con_trees[[i]][[j+1]]))
   }
   supercontree[[j+1]]<-.pgi.supercon.fill(supercontree[[j+1]],con_trees_tosend,num_trees,seq_len,verbosity)
  }
 }
 return(supercontree)
}

.pgi.to.phylo.fill <-
function(tree,st) {
	if(tree[[1]] == 0) {
		# node with 2 tips
		st[[6]]<-st[[6]]+1
  	st[[4]]<-st[[4]]+1
  	node<-st[[6]]
  	for(i in 2:length(tree)) {
  		edge1<-c(node,st[[5]])
			st[[5]]<-st[[5]]+1
			#new taxa, need name
			st[[2]]<-c(st[[2]],rownames(tree[[i]][[1]]))
	    st[[3]]<-c(st[[3]],print.rds(tree[[i]][[1]][1,]))
			st[[1]]<-rbind(st[[1]],edge1,deparse.level=0)
		}
	} else if(tree[[1]] == -1) { 
		st[[6]]<-st[[6]]+1
		st[[4]]<-st[[4]]+1
		node<-st[[6]]
		edges<-c()
		for(i in 2:length(tree)) {
			if(length(tree[[i]][[1]]) > 1) {
				edge1<-c(node,st[[5]])
				st[[5]]<-st[[5]]+1
				st[[2]]<-c(st[[2]],rownames(tree[[i]][[1]]))
		    st[[3]]<-c(st[[3]],print.rds(tree[[i]][[1]][1,]))
				st[[1]]<-rbind(st[[1]],edge1,deparse.level=0)
			} else {
				#node as child - increment node ctr
				edge1<-c(node,st[[6]]+1)
				st[[1]]<-rbind(st[[1]],edge1,deparse.level=0)
				st<-.pgi.to.phylo.fill(tree[[i]],st)
			}
		}
	}
	return(st)
}

.pgi.tree.hash<-function(tree) {
	t.sum<-0
	for(i in 2:length(tree)) {
		if(length(tree[[i]][[1]]) > 1) {
			if(dim(tree[[i]][[1]])[1] == 1) {
				#tip
				t.sum<-t.sum+(sum(tree[[i]][[1]])*(tree[[i]][[1]][1]+1))
			} else {
				#filled node
				t.sum<-t.sum+.pgi.tree.hash(tree[[i]])	
			}
		} else {
			t.sum<-t.sum+.pgi.tree.hash(tree[[i]]) 
		}
	}
	return(t.sum)
}

.prepare.con.tree.struct <-
function(filled_tree,seq_len,verbosity=0) {
	if(verbosity > 2) cat("Preparing the pseudoconsensus tree data structure, seq_len = ",seq_len,"\n",sep="")
	empty_matrix<-rep(0,times=(seq_len))
	dim(empty_matrix)<-c(1,seq_len)
	num_times<-0
	con_seq<-rep(0,times=seq_len+1)
	filled_tree[[1]]<-list(filled_tree[[1]])
	count<-(length(filled_tree[[1]][[1]][1,])-seq_len)-1
	filled_tree[[1]]<-c(filled_tree[[1]],list(list(con_seq,con_seq)))
	for(i in 1:count) {
		filled_tree[[1]]<-c(filled_tree[[1]],list(num_times,list(empty_matrix,empty_matrix)))
	}
	num_child<-(length(filled_tree[[1]][[1]][1,]) - seq_len)-1
	for(i in 2:(num_child+1)) {
		if(filled_tree[[i]][[1]][1,][seq_len+1] != -1) {
			filled_tree[[i]]<-.prepare.con.tree.struct(filled_tree[[i]],seq_len)
		} else {
			filled_tree[[i]]<-list(filled_tree[[i]])
		}
	}
	return(filled_tree)
}

.prep_seq <-
function(seq,name) {
 seq<-c(seq,-1,-1)
 dim(seq)<-c(1,length(seq))
 rownames(seq)<-name
 return(seq)
}

.refine.anc.seqs <-
function(anc_dev_seq, der_dev_seqs, inf.params ,cut_per,verbosity=0) {
 seq_len<-length(anc_dev_seq)
 anc_dev_seq_mat<-c()
 # counter for lack of progess (if no new seqs in cycles / 10 # of cycles, skip)
 if(inf.params$heuristic == "pgi") {
  no_prog_count<-0
  if(inf.params$cycles < 40) {
   no_prog_cut<-4
  } else {
   no_prog_cut<-round(inf.params$cycles/10)
  }
 }
	if(verbosity == 2) {
		pb.node<-txtProgressBar(min=1,max=inf.params$cycles,style=3)
	}
 for(i in 1:inf.params$cycles) {
  if(inf.params$heuristic == "pgi" && no_prog_count >= no_prog_cut) {
   if(verbosity > 2) print("No Progress - Stoping inference for this node")
   break
  }
  dev_seq_mat<-rbind(.gen_ran(anc_dev_seq,((inf.params$replicates)-1),inf.params$simultaneity),anc_dev_seq)
  dev_seq_mat<-unique(dev_seq_mat)
  if(length(anc_dev_seq_mat) > 0) {
   idx<-which(duplicated(rbind(anc_dev_seq_mat[,seq(1:seq_len)],dev_seq_mat))[-seq(1:dim(anc_dev_seq_mat)[1])])
   dev_seq_mat<-rbind(dev_seq_mat[-idx,],anc_dev_seq)
  }
  new_reps<-(dim(dev_seq_mat)[1])
  score_matrix<-rep(1000,times=new_reps)
  for(u in 1:length(der_dev_seqs)) {
   score_matrix<-cbind(score_matrix,rep(1000,times=new_reps))
  }
  for(j in 1:new_reps) {
   temp_score<-0
   temp_best<-0
   for(z in 1:(length(der_dev_seqs))) {
    # US penalty
    # For each child of the node
    temp_score2<-0
    if(der_dev_seqs[[z]][1,seq_len+1] == -1) {
     offset <- 1
    } else {
     offset <- 0
    }
    temp_best<-(edit.cost(dev_seq_mat[j,],der_dev_seqs[[z]][1,1:seq_len],0,2,inf.params$edit.cost.func,seq_len))+der_dev_seqs[[z]][1,length(anc_dev_seq)+1] + offset
    score_matrix[j,z+1]<-1
    # For number of derived sequences at the child node
    for(k in 1:dim(der_dev_seqs[[z]])[1]) {
     if(k == 1) {
      next
     }
     if(der_dev_seqs[[z]][1,length(anc_dev_seq)+1] == -1) { 
      offset <- 1
     } else {
      offset <- 0
     }
     temp_score2<-(edit.cost(dev_seq_mat[j,],der_dev_seqs[[z]][k,1:seq_len],0,2,inf.params$edit.cost.func,seq_len))+der_dev_seqs[[z]][k,length(anc_dev_seq)+1] + offset 
     if(temp_score2 < temp_best) {
      temp_best<-temp_score2
      score_matrix[j,z+1]<-k
     }  
    } # End mulptiple derived seqs

    # Add to the total from the other nodes
    temp_score<-temp_score+temp_best

   } # End the # of nodes
   score_matrix[j,1] <- temp_score
  }
  min_idx<-which.min(score_matrix[,1])
  min<-score_matrix[min_idx,1]
  aug_dev_seq_mat<-cbind(dev_seq_mat,score_matrix)
  temp_anc_dev_seq_mat<-aug_dev_seq_mat[aug_dev_seq_mat[,(length(anc_dev_seq)+1)] <= min]
  dim(temp_anc_dev_seq_mat)<-c((length(temp_anc_dev_seq_mat)/(length(anc_dev_seq)+length(der_dev_seqs)+1)),length(anc_dev_seq)+length(der_dev_seqs)+1)
  temp_anc_dev_seq_mat<-unique(temp_anc_dev_seq_mat)
  if(inf.params$heuristic == "pgi") {
   if(length(anc_dev_seq_mat) > 0) {
    before<-dim(anc_dev_seq_mat)[1]
   } else {
    before<-0
   }
  }
  anc_dev_seq_mat<-rbind(anc_dev_seq_mat,temp_anc_dev_seq_mat)
  anc_dev_seq_mat<-unique(anc_dev_seq_mat)
  if(inf.params$heuristic == "pgi") {
   if(dim(anc_dev_seq_mat)[1] > before) {
    if(verbosity > 2) cat("Cycle #:",i,"of",inf.params$cycles,"Progress: ",dim(anc_dev_seq_mat)[1]-before," new hypothetical ancestors for score:",min,"\n")
    no_prog_count<-0
   } else {
    no_prog_count<-no_prog_count+1
    if(verbosity > 2) cat("Cycle #:",i,"of",inf.params$cycles,"No progress\n")
   }
  } else if(inf.params$heuristic == "oldgenetic") {
    if(verbosity > 2) cat("Cycle #:",i,"of",inf.params$cycles,"Current best score:",min,"\n")
  }
	if(verbosity == 2) setTxtProgressBar(pb.node,i)
 
  # for the next round, choose a random seq from the set of the best ones
  best_seqs<-anc_dev_seq_mat[anc_dev_seq_mat[,seq_len+1] <= min]
  dim(best_seqs)<-c((length(best_seqs)/(seq_len+1+length(der_dev_seqs))),(seq_len+1+length(der_dev_seqs)))
  num_seqs<-dim(best_seqs)[1]
  min_idx<-round(runif(1,min=(dim(anc_dev_seq_mat)[1]-num_seqs)+1,max=dim(anc_dev_seq_mat)[1]))
  anc_dev_seq<-anc_dev_seq_mat[min_idx,][1:seq_len]
 }
	if(verbosity == 2) close(pb.node)
 width<-length(anc_dev_seq_mat[1,])
 f_min<-min(anc_dev_seq_mat[length(anc_dev_seq)+1])
 cutoff<-cut_per*f_min
 if(cutoff < 5) {
  cutoff<-5
 }
 n_anc_dev_seq_mat<-anc_dev_seq_mat[anc_dev_seq_mat[,length(anc_dev_seq)+1] < min+cutoff]
 height<-length(n_anc_dev_seq_mat)/width
 dim(n_anc_dev_seq_mat)<-c(height,width)
 return(n_anc_dev_seq_mat)
}

.trans.chars <-
function(chars) {
 char_mat<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y")
 chars[chars == "z"]<-"-2"
 for(i in 1:25) {
  chars[chars == char_mat[i]] <- i+9
 }
 return(chars)
}

.upsweep <-
function(filled_tree,idx,edit.cost.func,seq_len, sample.max=0,verbosity=0) {
	num_child<-(length(filled_tree[[1]][[1]][idx,]) - seq_len)-1
	ptrs<-filled_tree[[1]][[1]][idx,][(seq_len+2):((seq_len+2)+num_child-1)]
	anc_dev_seq<-filled_tree[[1]][[1]][idx,][1:seq_len]
	dim(anc_dev_seq)<-c(1,seq_len)
	best_score<-0
	child_dims<-c()
	if(edit.cost.func == "event-pair") t.pars<-"parsimov" else t.pars<-edit.cost.func
	if(verbosity > 2) cat("\nThis node has ",num_child," children\n",sep="")
	for(i in 2:(num_child+1)) {
		child_dims<-c(child_dims,dim(filled_tree[[i]][[1]][[1]])[1])
		best_score<-best_score+edit.cost(anc_dev_seq,filled_tree[[i]][[1]][[1]][ptrs[i-1],(1:seq_len)],0,2,edit.cost.func,seq_len)
	}
	if(verbosity > 2) cat("the best score for this node is ",best_score,"\n",sep="")
	if(verbosity > 2) cat("Generating matrix of pointers for this node... ")
	child.dim.list<-list()
	for(i in 1:num_child) child.dim.list<-c(child.dim.list,list(seq(1,child_dims[i])))
	ptr.matrix<-as.matrix(expand.grid(child.dim.list))
	if(verbosity > 2) cat("\n\nDimensions of ptr matrix: ",dim(ptr.matrix)," child dims = ",child_dims,"\n\n",sep="")
	if(verbosity > 2) cat("Done.\n")
	for(i in 1:dim(ptr.matrix)[1]) {
		t.score<-0
		for(j in 1:num_child) {
			t.score<-t.score+edit.cost(anc_dev_seq,filled_tree[[j+1]][[1]][[1]][ptr.matrix[i,j],(1:seq_len)],0,2,edit.cost.func,seq_len)
		}
		if(t.score == best_score) {
			if(verbosity > 3) print("found an equally good solution")
			filled_tree[[1]][[2]][[1]][seq_len+1]<-filled_tree[[1]][[2]][[1]][seq_len+1]+1
			for(j in 1:seq_len) {
				if(anc_dev_seq[j] == -2) {
					filled_tree[[1]][[2]][[2]][j] <- filled_tree[[1]][[2]][[2]][j] + 1
				} else {
					filled_tree[[1]][[2]][[1]][j] <- filled_tree[[1]][[2]][[1]][j]+anc_dev_seq[j]
				}
			}
			for(j in 1:num_child) {
				if(verbosity > 3) cat("found an equally good solution and now going through each child, num.child:",num_child,"j = ",j,"\n")
				sel_evn<-edit.cost(anc_dev_seq,filled_tree[[j+1]][[1]][[1]][ptr.matrix[i,j],(1:seq_len)],2,2,t.pars,seq_len)
				filled_tree[[1]][[2*j+1]]<-filled_tree[[1]][[2*j+1]]+1
				if(length(sel_evn) != 0) {
					for(k in 1:length(sel_evn)) {
						if(sel_evn[k] < 0) {
							filled_tree[[1]][[2*j+2]][[2]][-sel_evn[k]] <- filled_tree[[1]][[2*j+2]][[2]][-sel_evn[k]]+1
						} else {
							filled_tree[[1]][[2*j+2]][[1]][sel_evn[k]] <- filled_tree[[1]][[2*j+2]][[1]][sel_evn[k]]+1
						}
					}
				}
			}
			for(j in 1:num_child) {
				if(filled_tree[[j+1]][[1]][[1]][ptr.matrix[i,j],][seq_len+1] != -1) {
					if(sample.max != 0) {
						# Semi-exhaustive consensus proceedure.
						if(filled_tree[[1]][[2*j+1]] < sample.max) filled_tree[[j+1]]<-.upsweep(filled_tree[[j+1]],ptr.matrix[i,j],edit.cost.func,seq_len,sample.max=sample.max,verbosity=verbosity)
					} else {
						# fully exhaustive proceedure
						filled_tree[[j+1]]<-.upsweep(filled_tree[[j+1]],ptr.matrix[i,j],edit.cost.func,seq_len,sample.max=sample.max,verbosity=verbosity)
					}
				}
			}
		}
	}
 return(filled_tree)
}

.upsweep.simple <-
function (filled_tree,idx,edit.cost.func,seq_len) {
	num_child<-(length(filled_tree[[1]][[1]][idx,]) - seq_len)-1
	ptrs<-filled_tree[[1]][[1]][idx,][(seq_len+2):((seq_len+2)+num_child-1)]
	anc_dev_seq<-filled_tree[[1]][[1]][idx,][1:seq_len]
	dim(anc_dev_seq)<-c(1,seq_len)
	best_score<-0
	child_dims<-c()
	if(edit.cost.func== "event-pair") t.pars<-"parsimov" else t.pars<-edit.cost.func
	for(i in 2:(num_child+1)) {
		best_score<-best_score+edit.cost(anc_dev_seq,filled_tree[[i]][[1]][[1]][ptrs[i-1],][1:seq_len],0,2,edit.cost.func,seq_len)
		child_dims<-c(child_dims,dim(filled_tree[[i]][[1]][[1]])[1])
	}
	filled_tree[[1]][[2]][[1]][seq_len+1]<-filled_tree[[1]][[2]][[1]][seq_len+1]+1
	for(j in 1:seq_len) {
		if(anc_dev_seq[j] == -2) {
			filled_tree[[1]][[2]][[2]][j] <- filled_tree[[1]][[2]][[2]][j] + 1
		} else {
			filled_tree[[1]][[2]][[1]][j] <- filled_tree[[1]][[2]][[1]][j]+anc_dev_seq[j]
		}
	}
	sel_evn<-c()
	for(i in 2:(num_child+1)) {
		sel_evn<-c(sel_evn,list(edit.cost(anc_dev_seq,filled_tree[[i]][[1]][[1]][ptrs[i-1],][1:seq_len],2,2,t.pars,seq_len)))
	}
	for(z in 1:length(sel_evn)) {
		filled_tree[[1]][[2*z+1]]<-filled_tree[[1]][[2*z+1]]+1
		if(length(sel_evn[[z]]) != 0) {
			for(j in 1:length(sel_evn[[z]])) {
				if(sel_evn[[z]][j] < 0) {
					filled_tree[[1]][[2*z+2]][[2]][-sel_evn[[z]][j]] <- filled_tree[[1]][[2*z+2]][[2]][-sel_evn[[z]][j]]+1
				} else {
					filled_tree[[1]][[2*z+2]][[1]][sel_evn[[z]][j]] <- filled_tree[[1]][[2*z+2]][[1]][sel_evn[[z]][j]]+1
				}
			}
		}
	}
	for(i in 2:(num_child+1)) {
		if(filled_tree[[i]][[1]][[1]][ptrs[i-1],][seq_len+1] != -1) {
			filled_tree[[i]]<-.upsweep.simple(filled_tree[[i]],ptrs[i-1],edit.cost.func,seq_len)
		}
	} 
	return(filled_tree)
}

