# microsimr

## What it does?

microsimr provides functions for simulating microsatellite data under
a coalescent using Hudson's ms program. The user is able to analyze
the simulated data directly within R. There is also the option of
outputting the simulated data in a format suitable for Migrate-N.
Functions to output simulated data in other formats can be added
upon request (e.g., IMa2 or LAMARC). Currently, microsimr only simulates 
microsatellite data under a symmetrical strict stepwise mutation model, 
and a two-phase mutation model with symmetrical geometric distribution. Other 
mutation models can be considered if there is interest..

## How to get it?

To install `microsimr` you need Hadley Wickam's `devtools`:

    install.packages("devtools")
  
You can then install `microsimr` by:

    devtools::install_github("andersgs/microsimr")

If you want the vignettes to be installed, you should use:

    devtools::install_github("andersgs/microsimr", build_vignettes = TRUE)
  
But, this will take a little longer to install.

## A quick example

To simulate 100 loci at theta = 5 (per locus mutation rate scaled by the 
effective population size --- 4Nemu) genotyped across 20 individuals, with
10 individuals sampled from each of 2 populations exchanging migrants at the
rate of 15 migrants per generation (4Nem = 15):

    simd_data <- microsimr::sim_microsats(theta = 5,
                                      n_ind = c(10,10),
                                      n_loc = 100,
                                      n_pop = 2,
                                      ms_options = "-I 2 20 20 15")
                                      
The data can then be saved to a Migrate-N formatted input file thus:

    microsimr::save_microsimr(genotypes = simd_data,
                          filename = "infile",
                          n_ind = c(10,10),
                          n_pop = 2,
                          n_loc = 100,
                          format = 'migrate-n')
