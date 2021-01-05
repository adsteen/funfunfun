##' Functional redundancy via Ricotta's method
##'
##' @export
##' @param p # WHAT IS THIS
##' @param trait # WHAT IS THIS

ricotta_fr <- function(p, trait) {

  # Check inputs

  ## NOTE THE NAMESPACE DESIGNATION FOR VEGAN
  # Instead of adding vegan to the search path via library(vegan)
  # and then relying on R to figure out where the vegdist function is
  # (Properly "relying on R to identify the correct namespace for vegdist")
  # we specify it explicitly with vegan::vegdist


  t.dist <- as.matrix(vegan::vegdist(as.matrix(trait), method = "euclidean"))
  K.i <- rowSums(sweep(t.dist, MARGIN = 2, as.numeric(p), "*"))
  Q.rao <- sum(p * K.i)
  # Q.rao2<-sum(apply(sweep(t.dist,MARGIN = 2,as.numeric(p),"*"),MARGIN = 1,var)*p)
  U <- Q.rao / sum(p * (1 - p))
  R = 1 - U
  # var.max<-1-(1/length(p))^2*length(p)
  # R<-1-Q.rao/var.max
  d.RQ <- data.frame(R = R, U = U, Q = Q.rao)
  return(d.RQ)
}
