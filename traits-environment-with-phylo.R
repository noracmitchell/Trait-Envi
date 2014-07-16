require(R2jags)
require(mvtnorm)
require(ape)

rm(list=ls())

debug <- FALSE
plot <- TRUE
print <- TRUE
report.DIC <- TRUE

protea <- FALSE
schnitzler <- FALSE
pellie <- TRUE



gamma.rate.resid <- 1.0
gamma.shape.resid <- 1.0
gamma.rate.species <- 1.0
gamma.shape.species <-  1.0
beta.par <- 6
max.r <- 0.5


model.file="traits-environment-with-phylo.txt"


if (protea) {
  ## Note: NEXUS.file is used only to get species names if
  ## independent == TRUE
  ##
  if (schnitzler) {
    NEXUS.file <- "schnitzler.nex"
  } else {
    NEXUS.file <- "valentetree.nex"
  }
  CSV.file <- "Protea_CFR_DATA.csv"
} else if (pellie) {
  NEXUS.file <- "FINAL_20112012_ultrametric.nex"
  CSV.file <- "Pel_CFR_2011_2012.csv"
} else {
  stop("Either protea or pellie must be true\n")
}

if (debug) {
  n.chains <- 1
  n.burnin <- 500
  n.iter <- 1000
  n.thin <- 1
  ## to allow replication of results across runs in JAGS
  ##
  set.seed(1)
} else {
  n.chains <- 5
  n.burnin <- 5000
  n.iter <- 25000
  n.thin <- 25
}

standardize <- function(x) {
  if (is.numeric(x)) {
    y <- (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
  } else {
    y <- x
  }
  y
}

get.mean.vector <- function(x) {
  x.mean <- apply(x, 1, mean, na.rm=TRUE)
}

drop.levels <- function(dat) {
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}

likelihood <- function(y, mu, Sigma) {
  llike <- 0.0
  for (i in 1:nrow(y)) {
    llike <- llike + dmvnorm(y[i,], mu[i,], Sigma[,], log=TRUE)
  }
  llike
}

Dbar <- function(y, mu, Sigma) {
  dbar <- numeric(0)
  n.rep <- dim(mu)[1]
  for (i in 1:n.rep) {
    if ((i %% 10) == 0) {
      cat(".", sep="")
      flush.console()
    }
    if ((i %% 500) == 0) {
      cat(i, "\n")
      flush.console()
    }
    dbar[i] <- -2.0*likelihood(y, mu[i,,], Sigma[i,,])
  }
  mean(dbar)
}

Dhat <- function(y, mu.mean, Sigma.mean) {
  -2.0*likelihood(y, mu.mean, Sigma.mean)
}

## read in appropriate data set
##
traits <- read.csv(CSV.file,
                   na.strings=".",
                   header=TRUE)

## add climate data
##
climate <- read.csv("Protea_Pellie_Climate.csv",
                    na.strings=".",
                    header=TRUE)

combined <- merge(traits,
                  climate,
                  by.x="Site_location",
                  by.y="Site_location")
if (debug) {
  print(levels(combined$Site_location))
  cat("nrow(traits):   ", nrow(traits), "\n")
  cat("nrow(climate):  ", nrow(climate), "\n")
  cat("nrow(combined): ", nrow(combined), "\n")
}
rm(traits, climate)

## set up the phylogenetic correlation matrix
##
if (protea) {
  ## The Schnitzler et al. phylogeny includes only Protea.
  ## Nora's phylogeny includes other taxa as outgroups, so
  ## we need to set up a list of taxa to drop
  ##
  if (schnitzler) {
    drop.tips <- FALSE
  } else {
    drop.tips <- TRUE
    drop <- c("Sebarb",
              "Sevill",
              "Sescop",
              "Legand",
              "Farubr",
              "Fasali",
              "Faroc",
              "Fagalp",
              "Famacn")
  }
} else {
  drop.tips <- FALSE
}
## read the data file
##
phylo <- read.nexus(NEXUS.file)
## and drop tips if necessary
##
if (drop.tips) {
  phylo <- drop.tip(phylo, drop)
}
G <- vcv(phylo, corr=TRUE)
Ginv <- solve(G)
## no longer need phylo, discard it
##
rm(phylo)

## prepare the data for JAGS
##
## 1. pull the relevant data into columns
##
site <- combined$Site_location
if (protea) {
  if (schnitzler) {
    species <- combined$Schnitzler_name
  } else {
    species <- combined$Valente_name
  }
} else {
  if (pellie) {
    species <- combined$Species_name
  }
}

lma <- combined$LMA
area <- combined$Canopy_area
lwr <- combined$LWratio
fwc <- combined$FWC
map <- combined$MAP
mat <- combined$MAT
ratio <- combined$ratio
cdd <- combined$CDD
inso <- combined$Insolation
elev <- combined$Elevation

## 2. put them back in a data frame and include only complete
## cases

tmp <- data.frame(species=factor(species, rownames(Ginv)),
                    site=site,
                    lma=lma,
                    area=area,
                    lwr=lwr,
                    fwc=fwc,
                    map=map,
                    mat=mat,
                    ratio=ratio,
                    cdd=cdd,
		                inso=inso,
                    elev=elev)

if (debug) {
  print(rownames(G))
  check <- sample(nrow(tmp), 15, replace=FALSE)
  print(as.numeric(tmp$species[check]))
  print(tmp[check,])
  flush.console()
}

ok <- complete.cases(tmp)
tmp <- tmp[ok,]
tmp <- drop.levels(tmp)


## 3. pull the vectors back out and standardize for JAGS
##
species <- as.numeric(tmp$species)
site <- tmp$site
lma <- standardize(tmp$lma)
area <- standardize(tmp$area)
lwr <- standardize(tmp$lwr)
fwc <- standardize(tmp$fwc)
map <- standardize(tmp$map)
mat <- standardize(tmp$mat)
ratio <- standardize(tmp$ratio)
cdd <- standardize(tmp$cdd)
inso <- standardize(tmp$inso)
elev <- standardize(tmp$elev)


n.samp <- nrow(tmp)
n.species <- max(species)
n.dim <- n.species*4

## construct response matrix
##
y <- as.matrix(data.frame(lma,
                          area,
                          lwr,
                          fwc))

## parameters for Wishart prior
##
## nrow = ncol = # of parameters, i.e., ncol(y) == 4
## nu <- nrow + 2 makes it as vague as possible
## Note: nu > nrow + 1 required for distribution to be
##       non-degenerate
##
Omega <- diag(x=1.0, nrow=ncol(y), ncol=ncol(y))
nu <- nrow(Omega) + 2

## prior precision on regression coefficients
##
tau <- 0.1


jags.data <- c("species",
                 "y",
                 "map",
                 "mat",
                 "ratio",
                 "cdd",
		             "inso",
		             "elev",
                 "n.samp",
                 "n.species",
                 "n.dim",
                 "Ginv",
                 "tau",
                 "gamma.rate.resid",
                 "gamma.shape.resid",
                 "gamma.rate.species",
                 "gamma.shape.species",
                 "beta.par",
                 "max.r")
jags.par <- c("beta.map",
                "beta.mat",
                "beta.ratio",
                "beta.cdd",
		            "beta.inso",
		            "beta.elev",
                "mu",
                "rho.resid",
                "Sigma.resid",
                "rho.species",
                "Sigma.species")

fit <- jags(data=jags.data,
            inits=NULL,
            parameters=jags.par,
            model.file=model.file,
            n.chains=n.chains,
            n.burnin=n.burnin,
            n.iter=n.iter,
            n.thin=n.thin,
            DIC=TRUE,
            working.directory=".")

mu.mean <- fit$BUGSoutput$mean$mu
Sigma.mean <- fit$BUGSoutput$mean$Sigma.resid
mu <- fit$BUGSoutput$sims.list$mu
Sigma <- fit$BUGSoutput$sims.list$Sigma.resid
summary<-fit$BUGSoutput$summary

if (print) {
  opt.old <- options(width=120)
  filename <- paste("results-",
                    gsub(":", "-",
                         gsub(" ", "-", Sys.time())),
                    ".txt",
                    sep="")
  if (!debug) {
    sink(filename, split=TRUE)
  }
  
    cat("Using ", NEXUS.file, "\n\n", sep="")
  
  cat("gamma.rate.resid:    ", gamma.rate.resid, "\n")
  cat("gamma.shape.resid:   ", gamma.shape.resid, "\n")
  cat("gamma.rate.species:  ", gamma.rate.species, "\n")
  cat("gamma.shape.species: ", gamma.shape.species, "\n")
  cat("beta.par:            ", beta.par, "\n")
  cat("max.r:               ", max.r, "\n\n")
  print(fit, digits.summary=3)
  if (!debug) {
    sink()
  }
  if (report.DIC) {
    dbar <- Dbar(y, mu, Sigma)
    dhat <- Dhat(y, mu.mean, Sigma.mean)
    pD <- dbar - dhat
    DIC <- dbar + pD
    if (!debug) {
      sink(filename, append=TRUE, split=TRUE)
    }
    cat("\n",
        "Dbar: ", dbar, "\n",
        "Dhat: ", dhat, "\n",
        "pD:   ", pD, "\n",
        "DIC:  ", DIC, "\n")
    if (!debug) {
      sink()
    }
  }
  options(opt.old)
}

if (plot) {
  dev.new()
  old.par <- par(mfrow=c(2,2))
  for (i in 1:4) {
    x.plot <- tapply(y[,i], site, mean)
    y.plot <- tapply(mu.mean[,i], site, mean)
    plot(x.plot, y.plot, xlab="Observed", ylab="Predicted",
         pch=16, cex=0.75,
         main=c("LMA", "Area", "LWR", "FWC")[i])
  }
  par(old.par)
  dev.new()
  old.par <- par(mfrow=c(2,2))
  for (i in 1:4) {
    x.plot <- tapply(y[,i], site, mean)
    y.plot <- tapply(y[,i]-mu.mean[,i], site, mean)
    plot(x.plot, y.plot, xlab="Observed", ylab="Residual",
         pch=16, cex=0.75,
         main=c("LMA", "Area", "LWR", "FWC")[i])
  }
  par(old.par)
}

filename <- paste("results-",
                  gsub(":", "-",
                       gsub(" ", "-", Sys.time())),
                  ".Rsave",
                  sep="")
save(fit, file=filename)

