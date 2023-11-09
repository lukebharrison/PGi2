# PGi R package v2.1 README 
2023/08/07
PGI2: Parsimov-scoared Genetic Inference Method

This package implements the PGi algorithm for the analysis of the evolution of developmental sequences. Many supporting functions and the basic structure of the [vignette](inst/doc/RunningASmallAnalysis.pdf) are inspired or based on functions/vignettes available in the R APE package.

Please refer to the following article for background on the algorithm. This is also the prefered citation for PGi: 

[Harrison, Luke B. and Larsson, Hans C. E. (2008) 'Estimating Evolution of Temporal Sequence Changes: A Practical Approach to Inferring
Ancestral Developmental Sequences and Sequence Heterochrony', Systematic Biology, 57:3, 378 — 387](https://academic.oup.com/sysbio/article/57/3/378/1661823)

Please contact Luke Harrison if you have any comments/questions 
or bug reports etc. luke.harrison@mail.mcgill.ca

## INSTALLATION
To install PGi2, ensure an up-to-date version of R is installed. 
 
Installation from the assemebed R source pacakge:
1. Download the assemled PGi package "pgi2_2.1-1.tar.gz" in this repository
2. Open R, change to the directory containing the package and install the PGi2 package in R:

>install.packages("pgi2_2.1-1.tar.gz")

>libray(pgi2)

Installation from GitHub
1. Open R
>library(devtools)

>install_github("lukebharrison/PGi2")

>libray(pgi2)

## QUICKSTART

Please refer to the R documentation of this package and its functions, each available in the R documentation system after the package is loaded. The key functions are pgi.read.nexus(), pgi(), pgi.supercon(). 

Please also refer to the R vignette: [Running a Small Analysis](inst/doc/RunningASmallAnalysis.pdf) for a worked example of the execution of PGi2 using Velhagen's (1997) data set.

## GETTING DATA INTO PGI2

1. Enter the dataset into a NEXUS file with the tree topology (see vignette for a worked example)
Notes: enter the RANKED data as discrete characters. So for example, if the sequence is 1-2-4-3 (= A-B-D-C), then:  
character 1(A) = 1  
character 2(B) = 2  
character 3(C) = 4  
character 4(D) = 3   
For ranks between 10 and 35 use A,B,C,...,Y, in the nexus file and so on. Enter unscorable/missing data as Z.  
It is possible to have 36+ but the the pgi data structure must be manually edited.  
Make sure the taxon names have no spaces and aren't too long (it'll copy those too).  
Enter a single phylogenetic topology and save it in to the nexus file in the standard parenthetical format (newick).  

3. Load the data into PGi, using the function 
>my_tree<-pgi.read.nexus("filename")  
more details in ?pgi.read.nexus

## RUNNING PGI2

Run the PGi Algorithm in Interative Mode
>pgi(interactive=T)
and follow the instruction prompts

## FREQUENTELY ASKED QUESTIONS

### Q: Can PGi2 be run on a cluster or using multiple cores?
A: The implementation in PGi2 is based on the R parallel/foreach/doMC packages. In PGi2, only the computation of fully independent runs is parallelized. The number of runs is controlled by the nruns variable in the pgi() function: using a vector of two integers nruns=c(nseries,nparallel), where any integer greater than 1 will active parallel processing. The total number of runs is the product of nseries*nparallel. For a clustering-based systems, you would need to execute R/PGi2 on a single computational node only and use multiple cores of that node, not simply submit to the job scheduler.

### Q: What parameters of cycles/replicates/retained ancestral sequences should be used for a given data set?
A: Please refer to the [PGi article above in Systematic Biology](https://academic.oup.com/sysbio/article/57/3/378/1661823) for theoretical details; PGi is a heuristic algorithm, and the solution space of ancestral sequences and heterochronies increases with the number of nodes and sequence length. The objective of the analysis is to sample enough of the solution space for a given tree and set of sequences to arrive at a set of solutions from independent computations that are more or less equally parsimonious in terms of tree length (total number of sequence heterochronies, each run doesn't have to yield necessarily identical lengths, just close). This is to ensure enough sampling of the solution space. This is similar to a bayesian MCMC analysis where convergence of the MCMC chains is desired - although PGi is not a Markov-chain based method). For example, if 4 independent runs with fixed parameters are performed and the tree lengths and consensus solutions are highly divergent, then the analysis is not being run with adequate parameters (in terms of cycles/replicates/retained sequences).  
There are several concerns to also be aware of: if there is a high degree of simultaneity in the data sets, there will be many more possible solutions of similar length - this will require more computation and will also mean that you will discover many different possible solutions that, when summarized in a superconensus of the Indvidual runs, will likely not recover very many heterochronies with a high degree of support.  
Practically speaking: the analysis should be run with the largest values for cycles/replicates/ret.anc.seqs that are computationally feasible to achieve convergences of multiple independent runs. As an example, the Sanchez-Villagra data set (2002) of 24 events and 10 sequences was analyzed with cycles=150, replicates=150, ret.anc.seqs=150 with good results. The [Harington et al. (2013)](https://onlinelibrary.wiley.com/doi/10.1111/ede.12043) analysis of 25 elements and 30 sequences yielded accetable results at 100/100/100. These values do not need to be identical, but in practice, I've typically scaled them together.

## FOR FURTHER DETAILS 

Please consult the R documentation and the [PDF Vignette](inst/doc/RunningASmallAnalysis.pdf) available in the R package and on this repository for further details.

Don't hesitate to contact me (Luke Harrison) - luke.harrison@mail.mcgill.ca with any bugs, comments, suggestions or questions.






