# the epiBAT algorithm
This is an implementation of the epiBAT algorithm for SNP epistasis detection. The epiBAT algorithm is based on the bat algorithm and the tabu search. It uses two objective functions: K2-score and Gini score. Objectives are not combined, but are separately used in two populations. The best bats of both populations are then combined together to obtain the candidate set, which is evaluated for significance with modified G-test. The epiBAT algorithm is implemented in Python3 with numpy, pandas and scipy library.

#### Requirements
The epiBAT algorithm requires Python 3.6 with numpy, pandas and scipy library 

### Input parameters
Input parameters must be specified in the head of the epiBAT.py script. There are 6 parameters: 

`gini_population` represents the number of agents in a population of bats with Gini score as its objective    
`k2_population` represents the number of agents in a population of bats with K2 score as its objective    
`iteration_size` denotes the number of iterations for which the epiBAT algorithm will be running     
`alpha_value` denotes p-value threshold (before Bonferroni correction) which must be passed by the SNP combination to be said as significant    
`searching_path` path to the input file     
`output_file` path to the file, where the results of the epiBAT algorithm will be outputted    

Other optional parameters that can be specified by the user are:    
`freq_min` represents the minimum possible frequency of a bat    
`freq_max` represents the maximum possible frequency of a bat     
`alpha` represents the alpha parameter in the bat algorithm (should hold that 0<alpha<1)     
`gamma` represents the gamma parameter in the bat algorithm (should hold that 0<gamma)
`min_loudness_A` represents the minimum initial loudness of a bat     
`max_loudness_A` represents the maximum initial loudness of a bat    
`required_iterations_for_taboo` represents the number of iterations when the best solution in a population is not changed, that is needed to add that best solution to the tabu table   
`zeta_radius` defines the approximity when comparing solutions   

# How to run epiBAT
### Input data file format
Input data must be in a usual comma-delimited file containing the case-control genotype data. The first line in a file denotes the SNP IDs, whereas the last column denotes the class, i.e. case or control. The following rows contain the genotype data and the disease state, while the genotype data should have values 0, 1, 2 (i.e. the homozygous major allele, heterozygous allele, and homozygous minor allele), and the disease state should have values 0 and 1 (i.e. control and case). 
Example of the input data containing 5 SNPs and 3 samples is as follows:

X0,X1,X2,X3,X4,X5,Class      
0,0,0,1,0,0      
0,2,0,2,1,0      
1,0,0,1,1,1      
  
After setting the input parameters and preparing your input data, open the command line, go to the directory, where you have downloaded the epiBAT.py script and call:   
`epiBAT.py`

### Output of the epiBAT algorithm
The epiBAT algorithm outputs in a specified file results of its 1st stage (i.e. candidate set before evaluation by the G-test) and also final results (i.e. candidate set after evaluation by the G-test) with the corresponding p-value.
