##' Functional redundancy via Royalty's method
##'
##' @export
##' @description Calculates functional redundancy via the Royalty method, as described in INSERT REFERENCE HERE.
##' @param trait.contribution DESCRIBE THIS
##' @param q.range DESCRIBE THIS
##' @details ANYTHING THAT NEEDS TO GO INTO THE DETAILS SECTION
##' @return EXPLAIN WHAT IT RETURNS, WITH OBJECT TYPE
##' @examples # We should include examples, but I'm kind of confused about how examples work
##' @references ADD MANUSCRIPT REFERENCE HERE

royalty_fr <- function(trait.contribution, q.range = seq(0, 2, by = 0.1)) {

  # ADD SOME CODE TO MAKE SURE THAT PARAMS ARE VALID
  if(!is.vector(trait.contribution, mode = "numeric")) {
    stop("The trait.contribution argument must be a numeric vector")
  }
  if(!is.vector(q.range, mode = "numeric")) {
    stop("The q.range argument must be a numeric vector")
  }


  S <- length(trait.contribution)
    # vanish<-c(0.99999,rep(0.00001/(S-1),S-1))
    #Equation 3 and 4
    # trait.contribution<-abundance*trait/sum(abundance*trait)
    #Equation 6
    even.community <- rep(sum(trait.contribution) / length(trait.contribution),
                          length(trait.contribution))

    n.q <- length(q.range)
    indx <- 1

    fr_out <- data.frame(fr = rep(NaN, n.q),
                         D = rep(NaN, n.q),
                         q = rep(NaN, n.q))

    if (sum(even.community) > 0) {
      for (q in q.range) { # Pretty sure this can be vectorized
        if (q == 1) {
          D.q <-
            exp(-1 * sum(trait.contribution[trait.contribution > 0] * log(trait.contribution[trait.contribution >
                                                                                               0])))
          # delta.q=sum(trait.contribution*log(trait.contribution/even.community))
          # delta.0=sum(vanish*log(vanish/even.community))
          # fr<-1-delta.q/delta.0
          # S<-exp(-1*sum(even.community*log(even.community)))
        } else {
          D.q <- (sum(trait.contribution[trait.contribution > 0] ^ q)) ^ (1 / (1 -
                                                                                 q))
          # D.0<-(sum(vanish^q))^(1/(1-q))
          # delta.q=S-D.q#1/(q-1)*(D.q^(1-q)-S^(1-q))
          # delta.0=S-D.0#1/(q-1)*(D.0^(1-q)-S^(1-q))
          # fr<-1-delta/delta2
          # fr<-(D.q-1)/(S-1)
        }
        fr <- (D.q - 1) / (S - 1)
        fr_out$q[indx] <- q
        fr_out$D[indx] <- D.q
        fr_out$fr[indx] <- fr
        indx <- indx + 1
      }
    } else {
      fr_out$q[indx] <- q
      fr_out$D[indx] <- NaN
      fr_out$fr[indx] <- NaN
    }

    return(fr_out)
  }
















