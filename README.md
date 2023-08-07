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

Please also refer to the R vignette: RunningASmallAnalysis.pdf for a worked example of the execution of PGi2 using Velhagen's (1997) data set.

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

## FOR FURTHER DETAILS 

Please consult the R documentation and the [PDF Vignette](inst/doc/RunningASmallAnalysis.pdf) available in the R package and on this repository for further details.

Don't hesitate to contact me (Luke Harrison) - luke.harrison@mail.mcgill.ca with any bugs, comments, suggestions or questions.






