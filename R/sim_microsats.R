#' Simulate micrasatellite data using the coalescent
#'
#'@param theta The per locus mutation rate.
#'@param n_ind An integer indicating the number of diploid individuals to be
#'             simulated, or a vector indicating how many individuals to sample
#'             per simulated population
#'@param n_loc An integer indicating the number of loci to be simulated.
#'@param n_pop An integer indicating the number of populations to be simulated.
#'            If NULL, then assume n_pop = 1.
#'@param mutation_model A character string indicating the mutation model to use.
#'       Currently, only the strict stepwise mutation model of Ohta and Kimura (1973) ('smm'),
#'       and the DiRienzo et al. (1994) two-phase model ('tpm') are implemented.
#'       Default is 'smm'
#'@param sigma2 Variance in allele size to be used in the 'tpm' model
#'@param p_single Probability of a single-step mutation to be used in the 'tmp' model
#'@param ancestral_allele_size An integer indicating the ancestral allele
#'       size (i.e, the count of number of repeats of the ancestral allele).
#'       Default is 80 (see Details).
#'@param ms_options A string with additional options to pass on to ms.
#'
#' @details
#'
#' A tree is first simulated using 'ms'. Mutations are simulated along the
#' branches of the tree following a Poisson distribution with lambda proportional
#' to branch length times theta (4Nmu).
#'
#' The number of mutations along the branches are then transformed into gain/loss
#' of repeat units using \code{mutate_microsats} function.
#'
#' Details of each mutation model in the vignette.
#'
#' @references
#' Di Rienzo, A., Peterson, A. C., Garza, J. C., Valdes, A. M., Slatkin, M., & Freimer, N. B. (1994). Mutational processes of simple-sequence repeat loci in human populations. Proceedings of the National Academy of Sciences of the United States of America, 91(8), 3166–3170.
#'
#'Ohta, T., & Kimura, M. (2007). A model of mutation appropriate to estimate the number of electrophoretically detectable alleles in a finite population. Genetical Research, 89(5-6), 367–370. http://doi.org/10.1017/S0016672308009531
#'
#'@examples
#'# generate a simulated dataset for a single population with theta = 5,
#'# from which 100 loci were genotyped for 20 individuals
#'simd_data_1pop <- microsimr::sim_microsats(theta = 5,
#'                                            n_ind = 20,
#'                                            n_loc = 100,
#'                                            n_pop = 1)
#'
#'# generate a simulated dataset for a three populations with theta = 12,
#'# from which 100 loci were genotyped for 20 individuals each sampled from
#'# three populations with 20 migrants per generation
#'
#'simd_data_1pop <- microsimr::sim_microsats(theta = 12,
#'                                            n_ind = c(20, 20, 20),
#'                                            n_loc = 100,
#'                                            n_pop = 3,
#'                                            ms_options = "-I 3 40 40 40 20")
#'# Notice the 'ms_options' parameter, the -I indicates to model an island
#'# population structure model. The three 40's indicate to sample 40 chromosomes
#'# from each of the three modeled populations, and the 20 tells ms to simulate
#'# populations exchanging 20 migrants per generation. Additional options,
#'# including population growth, assymetric migration, isolation, etc. can
#'# be modeled by using the right combination of parameters. These are laid
#'# out in ms's manual. Examples can be added depending on user's needs. Please
#'# request additional examples through the GitHub issues page.
#'
#'@export

sim_microsats <- function(theta,
                          n_ind,
                          n_loc,
                          n_pop = NULL,
                          mutation_model = 'smm',
                          ancestral_allele_size = 80,
                          ms_options = NULL,
                          p_single = NULL,
                          sigma2 = NULL){
  #check if ms_options are ok
  #this section will have to be eventually expanded to allow for more
  #flexible creation of the ms command line, and for error checking.
  if (is.null(ms_options)){
    ms_options = "-T"
  } else {
    ms_options = paste("-T", ms_options)
  }
  #check how many populations should be simulated
  #if n_pop is not null, then leave it as is
  if (is.null(n_pop)) {
    n_pop = 1
  }
  #check if n_ind is an integer or vector, if vector, create total_ind parameter
  #as the sum of individuals to be simulated
  if (is.integer(n_ind)){
    total_ind = n_ind
  } else if (is.vector(n_ind) & length(n_ind) == n_pop) {
    total_ind = sum(n_ind)
  } else {
    stop("n_ind must be an integer or vector. If a vector, it must have length
         equal to n_pop.")
  }

  # checking tpm parameters
  if(mutation_model == 'tmp') {
    if(is.null(p_single)) {
      stop("Specified mutation_model 'tpm', but did not specify a probability of single-step mutations: p_single")
    }
    if(is.null(sigma2)) {
      stop("Specified mutation_model 'tpm', but did not specify variance of multiple-step sizes: sigma2")
    }
  }

  #generate ms trees using the phyclust library
  out_ms = phyclust::ms(nsam = 2 * total_ind, nreps = n_loc, opts = ms_options)
  #pick out only the lines that have trees in them
  trees = out_ms[which(grepl(pattern = "^\\(", out_ms))]
  #create a matrix to store genotypes
  #we need n_loc columns to store genotypes, plus 2 additional columns to
  #store population and individual information
  genotypes = matrix(character(1), ncol = (n_loc + 2), nrow = total_ind)
  pop_ix = rep(1:n_pop, n_ind)
  genotypes[,1] <- paste("pop", pop_ix, sep = "_")
  #generate sequential individual IDs that are exactly 10 characters long
  #this will allow me to generate migrate-n input files, which requires that
  #individual
  inds = sprintf("%-10s", paste("ind", 1:total_ind, sep = "_"))
  genotypes[,2] <- inds
  #a matrix to store estimates of thetas
  #using this for testing at the moment, should be removed later.
  sim_thetas = matrix(0, ncol = 1, nrow = length(trees))
  #cycle through simulated trees and generate genotypes for each locus
  #each tree represents a new locus
  for (t in 1:n_loc){
    #read the tree using the ape function read.tree
    tree_n = ape::read.tree(text = trees[t])
    #find the root of the tree
    root = setdiff(unique(tree_n$edge[,1]), unique(tree_n$edge[,2]))
    #simulate mutation events along the tree branches
    #it follows a Poisson distribution with parameter equal to theta times
    #the length of the branch
    #the output will generate a vector of integers of length equal to the
    #number of branches (i.e., edges) on the tree. Each integer indicates the
    #number of mutations that occurred on that branch
    sim_mutations = rpois(nrow(tree_n$edge), lambda = theta * tree_n$edge.length)
    #once we have the number of mutation events along each branch, we can
    #transform those mutation events into plus or minus one repeat (i.e.,
    #following a strict stepwise muation model), giving us a total number of
    #repeats gained (positive value) or lost (negative value) along a specific
    #branch. The loss or gain of a repeat is randomly assigned for each mutation
    #event
    #NOTE: here is were I will need to modify to incorporate additionaal mutation
    #models
    muts = mutate_microsats(n_mutations = sim_mutations,
                            mutation_model = mutation_model,
                            p_single = p_single,
                            sigma2 = sigma2)
    #create a vector of alleles, primed with the ancestral allele count of repeats.
    #The program will then crawl along the tree starting from the tip, and
    #working towards the root changing the allele size according to the
    #mutation events, by adding or subtracting repeats.
    alleles = rep(ancestral_allele_size, (2 * total_ind))
    #starting from the first tip, crawl across the tree from tip to root
    #adding the absolute change in number of repeats along the branches
    #that connect the tip to the root of the tree
    for (tip in 1: (2 * total_ind)){
      node = tip
      parent = tip
      while(parent != root){
        edge = which(tree_n$edge[,2] == node)
        muts[edge]
        alleles[tip] = alleles[tip] + muts[edge]
        parent = tree_n$edge[edge,1]
        node = parent
      }
    }
    #this is part of testing. might be removed later
    sim_thetas[t, ] <- 2 * var(alleles)
    #order the alleles so that they appear in order of population where they were
    #sampled order rather than in the order they appear on the tree
    alleles <- alleles[as.numeric(gsub(pattern = "s", replacement = "", x = tree_n$tip.label))]
    #use a "." as the allele separator, and collapse alleles that are in sequential
    #order, thus creating diploid genotypes
    als = paste(alleles, c("."," "), sep = "")
    als = apply(matrix(als, ncol = 2, byrow = T),1,paste, collapse = "")
    als = gsub(pattern = "\\s$", replacement = "", x = als)
    #add the genotypes to the genotypes matrix
    genotypes[,(t + 2)] <- als
  }
  return(genotypes)
}
