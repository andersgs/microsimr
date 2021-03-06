% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_input_files.R
\name{save_microsimr}
\alias{save_microsimr}
\title{Create input files}
\usage{
save_microsimr(genotypes, filename, n_ind, n_pop, n_loc, format = "migrate-n",
  title = "microsimr simulated data")
}
\arguments{
\item{genotypes}{A matrix outputted by sim_microsats function.}

\item{filename}{A character string indicating the filename to be used for saving
the input file.}

\item{n_ind}{An integer indicating the number of sampled diploid individuals
if n_pop = 1, or a vector of integers indicating the number of
diploid individuals sampled in each population in n_pop > 1.}

\item{n_pop}{An integer indicating the number of populations simulated.}

\item{n_loc}{An integer indicating the number of loci simulated.}

\item{format}{A character string indicating the format desired. Currently,
'migrate-n' is the only possible option.}

\item{title}{A character string to be used as the title of the input file,
only used if the format allows it.}
}
\description{
Create input files
}
\examples{
#simulate a dataset with theta = 20, 5 loci genotyped at a total of 10
# individuals, with 5 individuals taken from each of 2 populations. The
# populations are connected by 4 migrants/generation (4Nem = 4).

simd_data <- microsimr::sim_microsats(theta = 20,
                                     n_ind = c(5,5),
                                     n_loc = 5,
                                     n_pop = 2,
                                     ms_options = "-I 2 10 10 4")

# write an input file for migrate-n
\dontrun{
microsimr::save_microsimr(genotypes = simd_data,
                         filename = "infile",
                         n_ind = c(5,5),
                         n_pop = 2,
                         n_loc = 5,
                         format = 'migrate-n')
}

\dontrun{
# changing the format to anything else other than 'migrate-n' currently
generates an error
microsimr::save_microsimr(genotypes = simd_data,
                         filename = "infile",
                         n_ind = c(5,5),
                         n_pop = 2,
                         n_loc = 5,
                         format = 'ima2')
}

}

