##' Functional redundancy via Hester's method
##'
##' @export
##' @param trait.contribution #DESCRIBE THIS HERE
##'

hester_fr <- function(trait.contribution) {
  S <- length(trait.contribution)
  o = c() # Again, let's vectorize this
  for (i in 1:(S - 1)) {
    for (j in (i + 1):S) {
      if ((trait.contribution[i] > 0) & (trait.contribution[j] > 0)) {
        o <- c(o, 1)
      } else {
        o <- c(o, 0)
      }
    }
  }
  fr = median(o)
  return(fr)
}

