####Simulate eDNA Detection for given site occupancy and sampling effort (3 heirarchical levels: location/site, sample, subsample/PCR replicate)
### Note - number and definition of levels for eDNA occupancy are somewhat arbitrary - we are following eDNAOccupancy syntax here 

###0. Housekeeping
lapply(c("mvtnorm", "pROC", "mcmcse", "eDNAoccupancy"), require, character.only=T) ## load required packages

###1. User-Defined simulation variables
sim_n <- 10000 ## number of simulations
site_n <- 10 ## number of sites
site_occ_n <- 5 ## how many sites is the species really present at?
samp_n <- 6 ## number of samples taken at each site
pcr_n <- 3 ## number of pcr replicates per sample
samp_p <- 0.2 ## probability that a sample will have target eDNA in it, conditional on target being at site; NOTE - this parameter can vary by environment, species, sample volume, etc.
pcr_p <- 0.8 ## probability that a pcr will detect target eDNA, conditional on target being in sample; NOTE - this parameter can vary by pcr efficiency, inhibition levels, copy number, etc. 

###2. Create empty data.frame for each simulation
col_n <- 3+pcr_n ## need this to create the empty data frame
row_n <- site_n*samp_n ## need this to create the empty data frame
dat_sim <- data.frame(matrix(ncol=col_n, nrow=row_n)) # create the empty dataframe
frame_names <- c("Site", "Sample", "SampleOccupancy", paste("PCR", rep(1:pcr_n), sep="_")) #format
colnames(dat_sim) <- frame_names #format
dat_sim$Site <- rep(paste("Site", rep(1:site_n),sep=""), each=samp_n) #Fill in Site names
dat_sim$Sample <- rep(1:samp_n, site_n) #Fill in Sample number (important that this is an integer)
dat_sim$SampleOccupancy <- 0 #set all sample detections to zero by default
dat_sim[, 4:dim(dat_sim)[2]] = 0 #set all pcr detections to zero by default
site_occ_names <- paste("Site", rep(1:site_occ_n),sep="") #names of occupied sites for matching purposes

###3. Simulate your data
dat_sim[dat_sim$Site %in% site_occ_names, 3] <- rbinom(n=site_occ_n*samp_n, size=1, prob=samp_p) #if site is occupied, simulate sample occurrence with binomial distribution
samp_occ_n <- length(which(dat_sim$SampleOccupancy == 1)) #number of samples occupied
dat_sim[dat_sim$SampleOccupancy == 1, 4:(3+pcr_n)] <- rbinom(n=pcr_n*samp_occ_n, size=1, prob=pcr_p) #if sample is occupied, simulate pcr detections with binomial distribution

###4. Estimate occupancy from simulation results with eDNAoccupancy package (uses a Bayesian MCMC approach to fitting a 3-level occupancy model)
### ...for eDNAOccupancy help, see user manual with >vignette("eDNAoccIntro")
eDNADetections <- occData(dat_sim[, -3], siteColName='Site', sampleColName="Sample") #format input data (need to remove SampleOccupancy column)
mod.fit <- occModel(detectionMats=eDNADetections, niter=5000, niterInterval=1000) #fit the model; NOTE - you may be able to reduce the # of iterations (niter) but check convergence first (see section #6 for that)
posteriorSummary(mod.fit, burnin=1000, mcError=TRUE) # just a summary of results, good for checking parameter estimates, etc.; not necessary
psi <- posteriorSummaryOfSiteOccupancy(mod.fit, burnin=1000) # estimated probability of site occurrence
theta <- posteriorSummaryOfSampleOccupancy(mod.fit, burnin=1000) # estimated probability of conditional sample occurrence
p <- posteriorSummaryOfDetection(mod.fit, burnin=1000) # estimated probability of conditional pcr occurrence

###5. mash the true occupancy rates, mean observed occupancy rates, and estimated occupancy rates into output
observed_PCR_freq <-  sum(dat_sim[which(dat_sim$SampleOccupancy==1), 4:(3+pcr_n)])/(samp_occ_n*pcr_n)
output <- c(trueSite=site_occ_n/site_n, observedSite=site_occ_n/site_n , psi=psi$median[1], #True and observed site occupancy will be identical in this version
            trueSample=samp_p, observedSample=samp_occ_n/(site_occ_n*samp_n), theta=theta$median[1,1],  #true, observed, and estimated sample occupancy
            truePCR=pcr_p, observedPCR= observed_PCR_freq, p=p$median[1,1]) #true, observed, and estimated pcr occupancy
names(output)[3] <-"psi" 
names(output)[6] <-"theta" 
names(output)[9] <-"p" 
par(mar=c(2,8,2,1))
barplot(output, horiz=TRUE, las=1, col=rep(c("blue", "yellow1", "red"), each=3))

###6. Evaluate the model fit, convergence, autocorrelation; note - you don't need to run this everytime, but should check for convergence at given niter
# posteriorSummary(mod.fit, burnin=1000, mcError=TRUE)
# plotTrace(mod.fit, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), burnin=1000)
# plotACF(mod.fit, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), burnin=1000)
