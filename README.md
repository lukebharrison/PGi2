# PGi2
PGi v2.1 README 
2023/08/04
PGI2: Parsimov-scoared Genetic Inference Method

Please cite PGI as:

Harrison, Luke B. and Larsson, Hans C. E. (2008) 'Estimating
Evolution of Temporal Sequence Changes: A Practical Approach to Inferring
Ancestral Developmental Sequences and Sequence Heterochrony', Systematic
Biology, 57:3, 378 — 387


*** Note: The parallel implementaiton is not working on all computers ***


Please contact Luke Harrison if you have any comments/questions 
or bug reports etc.
@ 
luke.harrison@mail.mcgill.ca

To run a PGi analysis, change to the analysis directory and install the pgi package in R:
>install.packages("pgi2_2.1-1.tar.gz")
then load it:
>libray(pgi2)


1. Enter the dataset into a NEXUS file with the tree topology
Notes: enter the RANKED data as discrete characters. So for example, if the sequence is ABDC, then:
character 1 = 1
character 2 = 2
character 3 = 4
character 4 = 3 
and so on.

Enter unscorable/missing data as Z.

For ranks between 10 and 35 use A,B,C,...,Y, in the nexus file and so on.
It is possible to have 36+ but the the pgi data structure must be manually edited.

Make sure the taxon names have no spaces and aren't too long (it'll copy those too).
Enter a single phylogenetic topology and save it in to the nexus file
in the standard parenthetical format (newick).

2. Load the data into PGi, using the function 
my_tree<-pgi.read.nexus("filename")  
more details in ?pgi.read.nexus

-OR- 
Load one of the datasets in datasets using:
> data(velhagen)

Run the PGi Algorithm in Interative Mode
>pgi(interactive=T)
and follow the instructions

Alternatively, please read the manual files ?pgi to execute the algorithm non-interactively.






