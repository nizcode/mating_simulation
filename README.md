# Random Mating Simulation
Starting at grandparents(Generartion1) to probands(Generation3)
The mating is done at random hence all three generations always stay at Hardy Weinberg equilibrium
This can be used as a base for further simulations with different tweeks(ie. dissortative mating)

# Starting Allele frequencies 
0(AA) = 0.25, 1(Aa) = 0.5, 2(aa) = 1
And they stay the same for males and females
Male:Female ratio stays the same for all 3 generations
Generation1 size = N
Gen2 size = N/2
Gen3 size < N/4 
Generation 3 changes slightly beacuse there might not be equal mum to fathers ratio
and because of dissortative mating not all males will find their female


# simulation.R

Input
arg1 = sample population size, arg2 = the number of simulation you want, arg3 = 'yes' if you want disassortative mating or 'no' if you dont
Output
3 files in text format, one for each generation. 
MM MN NN pval chisq heterozygosity, these are the columns and the rows is the simulation number
There is also an .Rdata file which is written
The name of ouput file look like -> Random_sim_gen(arg1)_(arg2)_(arg3).txt/.Rdata

#dissort() function
if args3 = 'yes' then the dissort() function will be executed on 2nd generation only
all males will be paired with females of opposite genotypes 
ie, if MM(0) it can only be paired with a female of genotype Mm(1) or mm(2)
for 1 = 2 or 0
and for 2 = 1 or 0
I took reference from Chapter-3---Systems-of-Mating_2019_Human-Population-Genetics-and-Genomics.pdf




