#' A function to mutate microsatellites
#'
#' @param n_mutations A vector indicating the total number of mutations along
#'        each branch of a tree
#' @param mutation_model A character string indicating the mutation model to use.
#'       Currently, only the stict stepwise mutation model is implemented ('smm').
#'       Default is 'smm'

mutate_microsats <- function(n_mutations, mutation_model = 'smm'){
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
  } else {
    stop("Currently, the only mutation model that is implemented is the strict
         stepwise mutation model. Please use 'smm' for this option.")
  }
}
