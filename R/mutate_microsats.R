#' A function to mutate microsatellites
#'
#' @param n_mutations A vector indicating the total number of mutations along
#'        each branch of a tree
#' @param mutation_model A character string indicating the mutation model to use.
#'       Currently, only the strict stepwise mutation model of Ohta and Kimura (1973) ('smm'),
#'       and the DiRienzo et al. (1994) two-phase model ('tpm') are implemented.
#'       Default is 'smm'
#' @param sigma2 Variance in allele size to be used in the 'tpm' model
#' @param p_single Probability of a single-step mutation to be used in the 'tmp' model
#'
#' @details
#'
#' A tree is first simulated using 'ms'. Mutations are simulated along the
#' branches of the tree following a Poisson distribution with lambda proportional
#' to branch length times theta (4Nmu). This first part is done in \link{sim_microsats}.
#'
#' Here, the number of mutations at each branch is transformed to either a loss
#' or gain in number of repeats. In the 'smm' model, each mutation represents either
#' a loss or gain of a single repeat unit with equal probability. In the 'tpm'
#' model, with probability \code{p_single} a mutation represents either a
#' gain or loss of a single repeat, and with probability (1 - \code{p_single})
#' the gain/loss is larger following a symmetric geometric distribution.
#'
#' Please refer to the vignette to see a deeper explanation, and test of each
#' model.
#'
#' @references
#' Di Rienzo, A., Peterson, A. C., Garza, J. C., Valdes, A. M., Slatkin, M., & Freimer, N. B. (1994). Mutational processes of simple-sequence repeat loci in human populations. Proceedings of the National Academy of Sciences of the United States of America, 91(8), 3166–3170.
#'
#'Ohta, T., & Kimura, M. (2007). A model of mutation appropriate to estimate the number of electrophoretically detectable alleles in a finite population. Genetical Research, 89(5-6), 367–370. http://doi.org/10.1017/S0016672308009531

mutate_microsats <- function(n_mutations, mutation_model = 'smm', p_single = 0.8, sigma2 = 50){
  if (mutation_model == 'smm') {
    sapply(n_mutations, function(n_mut) {
      #if there are no mutations, return zero
      if (n_mut == 0){
        0
      } else {
        #if n_mut > 0, then generate n_mut random variables (rv) from a Uniform[0,1]
        #distribution. For rv < 0.5, return the loss of one repeat (-1), else
        #return the gain of one repeat (+1).
        #Return the sum of losses and gains of repeats for the absolute change in
        #repeat length at the microsat for a particular branch.
        sum(ifelse(runif(n_mut) < 0.5, -1, 1))
      }
    }
    )
  } else if (mutation_model == 'tpm') {
    sapply(n_mutations, function(n_mut) {
      if (n_mut == 0) {
        0
      } else {
        single_mutations <- rbinom( n = n_mut, size = 1, prob = p_single)
        total_single_step <- sum( ifelse(runif(n_mut) < 0.5, -1, 1) * single_mutations )
        total_multiple_step <- sum( rgeom_dist( n_mut, sigma2 = sigma2 ) * (1 - single_mutations ) )
        sum( total_single_step, total_multiple_step )
      }
    }
    )
  } else {
      stop("Currently, the only mutation model that is implemented is the strict
         stepwise mutation model. Please use 'smm' for this option.")
  }
}
