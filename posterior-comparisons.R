## Note: This code assumes that the posterior distributions for
## protea and pelargonium are already available in the corresponding
## objects. It also assumes that the number of posterior samples is the
## same for the two analyses
##
rm(list=ls())
require(ggplot2)

load("results-2014-03-06-09-01-02.Rsave")
protea <- fit$BUGSoutput$sims.list

load("results-2014-03-08-13-25-09.Rsave")
pellie <- fit$BUGSoutput$sims.list

## Posterior Comparison ##

traits <- c("LMA", "AREA", "LWR", "FWC")
betas <- c("cdd", "elev", "inso", "map", "mat", "ratio")

## construct the data frame for plotting
##
genus <- character(0)
covar <- character(0)
resp  <- character(0)
value <- numeric(0)
n.sample <- length(protea$beta.cdd[,1])
for (beta in betas) {
  for (trait in traits) {
    idx <- grep(trait, traits)
    coeff <- paste("beta.", beta, sep="")
    genus <- c(genus, rep("Protea", n.sample))
    covar <- c(covar, rep(beta, n.sample))
    resp  <- c(resp, rep(trait, n.sample))
    value <- c(value, protea[[coeff]][,idx])
    genus <- c(genus, rep("Pelargonium", n.sample))
    covar <- c(covar, rep(beta, n.sample))
    resp  <- c(resp, rep(trait, n.sample))
    value <- c(value, pellie[[coeff]][,idx])
  }
}

combined <- data.frame(genus=genus,
                       covar=covar,
                       resp=resp,
                       value=value)

for (trait in traits) {
  tmp <- subset(combined, resp==trait)
  p <- ggplot(tmp, aes(x=value, fill=genus)) +
    scale_fill_manual(values = c("black", "white") )+
       geom_density(alpha=0.3) +
       facet_wrap(~ covar) +
       geom_vline(yintercept=0, linetype="dashed", alpha=0.4) +
       ggtitle(paste("Posterior distribution of betas for", trait))
  print(p)
}



