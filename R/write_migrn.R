#' Create input file for Migrate-N
#'
#'@description Not a function for general use, please use save_microsimr()
#'              instead

#'@param genotypes A matrix outputted by sim_microsats function.
#'@param filename A character string indicating the filename to be used for saving
#'        the input file.
#'@param n_ind An integer indicating the number of sampled diploid individuals
#'       if n_pop = 1, or a vector of integers indicating the number of
#'       diploid individuals sampled in each population in n_pop > 1.
#'@param n_pop An integer indicating the number of populations simulated.
#'@param n_loc An integer indicating the number of loci simulated.
#'@param title A character string to be used as the title of the input file,
#'        only used if the format allows it.

write_migrn <- function(genotypes, filename, n_ind, n_pop, n_loc, title) {
  pop_refs = genotypes[,1]
  pop_ids = unique(pop_refs)
  total_cols = ncol(genotypes)
  header = paste(" ", n_pop, n_loc, ".", title, "\n")
  cat(header, file = filename)
  for (p in 1:n_pop){
    pop = paste(n_ind[p], "Population", p, "\n")
    cat(pop, file = filename, append = T)
    cat(apply(genotypes[pop_refs == pop_ids[p], 2:total_cols], 1, paste, collapse = " "),
        file = filename, sep = "\n", append = T)
  }
}
