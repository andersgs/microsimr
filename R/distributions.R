#' DiRienzo Symmetric Geometric distribution
#'
geom_dist <- function(j, p = NULL, sigma2 = NULL ) {
  if (is.null( p ) && !is.null(sigma2)) {
    p = (-1 + sqrt( 1 + 4*sigma2 ))/( 2 * sigma2 )
  }
  return( (1/2) * (p ) * (1-p)^(j - 1) )
}

#' Generate random numbers from the DiRienzo Symmetric Geometric distribution
#'
rgeom_dist <- function( n, p = NULL, sigma2 = NULL ) {
  rv = runif( n = n, min = 0, max = 0.5 )
  out = numeric( n )
  for(i in 1:n) {
    u = rv[i]
    j = 1
    pk = geom_dist( j = j, p = p, sigma2 = sigma2 )
    while( pk < u ) {
      j = j + 1
      pk = pk + geom_dist( j = j, p = p, sigma2 = sigma2 )
    }
    if( rbinom( n = 1, size = 1, prob = 0.5 ) == 0 ) {
      j = -j
    }
    out[i] <- j
}
  return( out )
}
