#' Create input files
#'
#'@param genotypes A matrix outputted by sim_microsats function.
#'@param filename A character string indicating the filename to be used for saving
#'        the input file.
#'@param n_ind An integer indicating the number of sampled diploid individuals
#'       if n_pop = 1, or a vector of integers indicating the number of
#'       diploid individuals sampled in each population in n_pop > 1.
#'@param n_pop An integer indicating the number of populations simulated.
#'@param n_loc An integer indicating the number of loci simulated.
#'@param format A character string indicating the format desired. Currently,
#'        'migrate-n' is the only possible option.
#'@param title A character string to be used as the title of the input file,
#'        only used if the format allows it.
#'
#'@examples
#'#simulate a dataset with theta = 20, 5 loci genotyped at a total of 10
#'# individuals, with 5 individuals taken from each of 2 populations. The
#'# populations are connected by 4 migrants/generation (4Nem = 4).
#'
#'simd_data <- microsimr::sim_microsats(theta = 20,
#'                                      n_ind = c(5,5),
#'                                      n_loc = 5,
#'                                      n_pop = 2,
#'                                      ms_options = "-I 2 10 10 4")
#'
#'# write an input file for migrate-n
#'\dontrun{
#'microsimr::save_microsimr(genotypes = simd_data,
#'                          filename = "infile",
#'                          n_ind = c(5,5),
#'                          n_pop = 2,
#'                          n_loc = 5,
#'                          format = 'migrate-n')
#'}
#'
#'\dontrun{
#'# changing the format to anything else other than 'migrate-n' currently
#' generates an error
#'microsimr::save_microsimr(genotypes = simd_data,
#'                          filename = "infile",
#'                          n_ind = c(5,5),
#'                          n_pop = 2,
#'                          n_loc = 5,
#'                          format = 'ima2')
#'}
#'
#'@export

save_microsimr <- function(genotypes, filename, n_ind, n_pop, n_loc,
                           format = 'migrate-n',
                           title = "microsimr simulated data") {
  if (format == 'migrate-n') {
    write_migrn(genotypes, filename, n_ind, n_pop, n_loc, title)
  } else {
    stop("microsimr does not yet recognize this format. please use 'migrate-n'
         or request in GitHub issues the desired format.")
  }

}
