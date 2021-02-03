##' Results of functional redundancy simulations
##' @export
##' @description Performs extinction simulations and calculates functional redundancy using the Royalty method. 
##' @param range.express The number of taxa contributing to a community-aggregated parameter--vector
##' @param range.nonexpress The number of taxa not contributing to a community-aggregated parameter--vector
##' @param q The diversity order to use when calculate functional redundancy--real number
##' @param a.range A parameter modulating how even a lognormal abundance distribution is--vector
##' @param loss.thres A parameter specifying the total loss trait in a community necessary to stop a simulation--number between 0 and 1
##' @param range.trait Possible range of trait values to be assigned to taxa--a vector
##' @param n.community Number of artificial communities to generate--integer greater than 0
##' @param n.extinct.sim Number of extinction simulations to perform on every artificial community--integer greater than 0
##' @details ANYTHING THAT NEEDS TO GO INTO THE DETAILS SECTION
##' @return RReturns a list of fr (functional redundancy), loss.frac (average fraction of a community required to go extinct prior to trait loss), and express.level (the fraction of a community required to go extinct)
##' @examples # We should include examples, but I'm kind of confused about how examples work

extinction_simulation <- function(range.express=c(1:100), range.nonexpress=c(0:100), q = 0.5, a.range = seq(0,0.4,by=0.01), loss.thres = 0.99, range.trait = seq(1,1000,by=1), n.community = 100, n.extinct.sim = 10) {#}, ncore = 1) {
  
  
  community.trait <- data.frame(fr=rep(NaN,n.community),
                                loss.frac=rep(NaN,n.community),
                                express.level=rep(NaN,n.community))
  # cl<-parallel::makeCluster(2) #create cluster object
  # doParallel::registerDoParallel(cores = ncore) #assign number of cores
  
  #parallelize extinction simulations
  
  # community.trait<-foreach::foreach(i=1:n.community, .combine = rbind) %dopar% {
  for (i in 1:n.community){ 
    #define artifical community
    n.e<-sample(range.express,1) #randomly choose number of trait possesing taxa
    n.n<-sample(range.nonexpress,1) #randomly choose taxa not possessing trait
    n.tot=n.e+n.n #calculate community richness
    S0<-1 
    R<-c(0:(n.tot-1)) #octaves
    a<-sample(a.range,1) #choose a parameter for generating lognormal community
    SR<-S0*exp(-a^2*R^2) #generate abundances of lognormal community
    relative.abundance<-matrix(SR/sum(SR), nrow = 1) #convert to proportions
    colnames(relative.abundance) <- 1:n.tot
    rownames(relative.abundance) <- 1
    
    #randomize trait levels among taxa
    trait.rand<-sample(range.trait,n.e,replace = TRUE) #randomly choose, with replacement, trait levels
    trait.rand<-matrix(sample(c(trait.rand,rep(0,n.n)),n.tot,replace = FALSE), ncol = 1) #randomize the order of trait levels; random trait vector is concatenated with a 0's vector of length n.n
    rownames(trait.rand) <- 1:n.tot
    colnames(trait.rand) <- 1
    
    contribution.dist <- t(trait.rand) * relative.abundance
    contribution.dist <- contribution.dist / sum(contribution.dist)
    
    n.member.loss<-rep(NaN,n.extinct.sim) #a dataframe for determining the proportion of a community that goes extinct
    
    #perform extinction simulations,it2 replicates 
    for (k in 1:n.extinct.sim){
      trait.loss <- 0 #intial trait loss from community aggregated parameter
      tmp.contribution <- contribution.dist #create dummy data.frame which updates after each iteration
      member.loss <- 0
      while(trait.loss < loss.thres) { #keep looping until trait level removed exceeds threshold
        indx <- sample(1:length(tmp.contribution),1) #randomly select community member
        trait.loss <- trait.loss + tmp.contribution[indx] #calculate total trait level removed from community aggregated parameter
        tmp.contribution <- tmp.contribution[-indx] #remove member from community
        member.loss <- member.loss+1 #monitor number of extinctions
      }
      n.member.loss[k] <- member.loss #assign total members loss during replicate to data.frame 
    }
    
    fr<-royalty_fr(relative.abundance,trait.rand,q)[1,1] #caluclate functional redundancy using method here
    
    loss.frac=mean(n.member.loss)/n.tot #calculate proportion of community that went extinct
    community.trait[i,]<-data.frame(fr=fr,
                                    loss.frac=loss.frac,
                                    express.level=n.e/n.tot) #create data.frame used for appending replicate simulations together, express level is proportion of community with a trait
    
  }
  
  # parallel::stopCluster(cl)#stop cluster object
  return(community.trait)
}
#write model outputs
# write.csv(community.trait,file =filename,quote = FALSE,row.names = FALSE)

