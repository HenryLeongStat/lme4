### suppressPackageStartupMessages(...)  as we have an *.Rout.save to Rdiff against
stopifnot(suppressPackageStartupMessages(require(lme4)))
options(show.signif.stars = FALSE, useFancyQuotes=FALSE)

source(system.file("test-tools-1.R", package = "Matrix"))# identical3() etc
all.EQ <- function(u,v, ...) all.equal.X(u, v, except = c("call", "frame"), ...)
S4_2list <- function(obj) {   # no longer used
   sn <- slotNames(obj)
   structure(lapply(sn, slot, object = obj), .Names = sn)
}

if (lme4:::testLevel() <= 1)
    quit("no")
## otherwise *print* normally:

oldOpts <- options(digits=2)
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm1a <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, REML = FALSE))
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1, fm2)

## Now works for glmer
fm1. <- suppressWarnings(glmer(Reaction ~ Days + (Days|Subject), sleepstudy))
## default family=gaussian/identity link -> automatically calls  lmer()  (but with a warning)
## hack call -- comes out unimportantly different
fm1.@call[[1]] <- quote(lmer)
stopifnot(all.equal(fm1, fm1.))
## Test against previous version in lmer1 (using bobyqa for consistency)
#(fm1. <- lmer1(Reaction ~ Days + (Days|Subject), sleepstudy, opt = "bobyqa"))
#stopifnot(all.equal(fm1@devcomp$cmp['REML'], fm1.@devcomp$cmp['REML']),
#          all.equal(fixef(fm1), fixef(fm1.)),
#          all.equal(fm1@re@theta, fm1.@theta, tolerance = 1.e-7),
#          all.equal(ranef(fm1), ranef(fm1.)))

## compDev = FALSE no longer applies to lmer
## Test 'compDev = FALSE' (vs TRUE)
## fm1. <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy,
##              compDev = FALSE)#--> use R code (not C++) for deviance computation
## stopifnot(all.equal(fm1@devcomp$cmp['REML'], fm1.@devcomp$cmp['REML']),
##           all.equal(fixef(fm1), fixef(fm1.)),
##           all.equal(fm1@re@theta, fm1.@re@theta, tolerance = 1.e-7),
##           all.equal(ranef(fm1), ranef(fm1.), tolerance = 1.e-7))

stopifnot(
    all.equal(fixef(fm1), fixef(fm2), tolerance = 1.e-13)
   ,
    all.equal(unname(fixef(fm1)),
              c(251.405104848485, 10.467285959595), tolerance = 1e-13)
   ,
    all.equal(Matrix::cov2cor(vcov(fm1))["(Intercept)", "Days"],
              -0.1375, tolerance = 4e-4)
)

fm1ML <- refitML(fm1)
fm2ML <- refitML(fm2)
(cbind(AIC= c(m1= AIC(fm1ML), m2= AIC(fm2ML)),
       BIC= c(    BIC(fm1ML),     BIC(fm2ML))) -> ICm)
stopifnot(all.equal(c(ICm), c(1763.94, 1762, 1783.1, 1777.97),
                    tolerance = 1e-5))# see 1.2e-6

(fm3 <- lmer(Yield ~ 1|Batch, Dyestuff2))
stopifnot(all.equal(coef(summary(fm3)),
		    array(c(5.6656, 0.67838803150, 8.3515624346),
			  c(1,3), dimnames = list("(Intercept)",
				  c("Estimate", "Std. Error", "t value")))))
showProc.time() #

### {from ../man/lmer.Rd } --- compare lmer & lmer1 ---------------
(fmX1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm.1 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))

#(fmX2 <- lmer2(Reaction ~ Days + (Days|Subject), sleepstudy))
#(fm.2 <- lmer2(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
## check update(<mer>, <formula>):
fm.3 <- update(fmX1, . ~ Days + (1|Subject) + (0+Days|Subject))
stopifnot(all.equal(fm.1, fm.3))

fmX1s <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy )# no longer:, sparseX=TRUE)
#fmX2s <- lmer2(Reaction ~ Days + (Days|Subject), sleepstudy, sparseX=TRUE)
options(oldOpts)  ## restore digits

showProc.time() #

for(nm in c("coef", "fixef", "ranef", "sigma",
	     "model.matrix", "model.frame" , "terms")) {
    cat(sprintf("%15s : ", nm))
    FUN <- get(nm)
    F.fmX1s <- FUN(fmX1s)
#    F.fmX2s <- FUN(fmX2s)
#    if(nm == "model.matrix") {
#        F.fmX1s <- as(F.fmX1s, "denseMatrix")
#        F.fmX2s <- as(F.fmX2s, "denseMatrix")
#	FF <- function(.) {r <- FUN(.); row.names(r) <- NULL
#			   as(r, "generalMatrix") }
#    } # else
    FF <- FUN
    stopifnot(
	      all.equal( FF(fmX1), F.fmX1s, tolerance =  1e-6)
#	      ,
#	      all.equal( FF(fmX2), F.fmX2s, tolerance = 1e-5)
#              ,
#	      all.equal( FF(fm.1), F.fmX2s, tolerance = 9e-6) ## these are different models
#              ,
#              all.equal(F.fmX2s,   F.fmX1s, tolerance = 6e-6)
#              ,
#              all.equal(FUN(fm.1), FUN(fm.2), tolerance = 6e-6)
              ,
              TRUE)
    cat("[Ok]\n")
}


## transformed vars should work[even if non-sensical as here;failed in 0.995-1]
fm2l <- lmer(log(Reaction) ~ log(Days+1) + (log(Days+1)|Subject),
             data = sleepstudy, REML = FALSE)
## no need for an expand method now : xfm2 <- expand(fm2)

stopifnot(dim(ranef(fm2l)[[1]]) == c(18, 2),
          is((c3 <- coef(fm3)), "coef.mer"),
          all(fixef(fm3) == c3$Batch),## <-- IFF  \hat{\sigma^2} == 0
          TRUE)



## Simple example by Andrew Gelman (2006-01-10) ----
n.groups <- 10 ; n.reps <- 2
n <- length(group.id <- gl(n.groups, n.reps))
## simulate the varying parameters and the data:
set.seed(0)
a.group <- rnorm(n.groups, 1, 2)
y <- rnorm (n, a.group[group.id], 1)
## fit and summarize the model
fit.1 <- lmer (y ~ 1 + (1 | group.id))
oldOpts <- options(digits=3)
coef (fit.1)
options(oldOpts)
## check show( <"summary.mer"> ):
(sf1 <- summary(fit.1)) # --> now looks as for fit.1

stopifnot(all.equal(fixef(fit.1), c("(Intercept)" = 1.571312129)),
	  all.equal(unname(ranef(fit.1, drop=TRUE)[["group.id"]]), structure(
		   c(1.8046888, -1.8097665, 1.6146451, 1.5408268, -0.1331995,
                     -3.3306655, -1.8259277, -0.8735145, -0.3591311,  3.3720441),
                   postVar = rep.int(0.311091076, 10)),
		    tolerance = 1e-5)
	  )


## ranef and coef
rr <- ranef(fm1)
stopifnot(is.list(rr), length(rr) == 1, is.data.frame(rr[[1]]))
print(plot(rr))
stopifnot(is(cc <- coef(fm1), "coef.mer"),
	  is.list(cc), length(cc) == 1, is.data.frame(cc[[1]]))
print(plot(cc))
rr <- ranef(fm2)
stopifnot(is.list(rr), length(rr) == 1, is.data.frame(rr[[1]]))
print(plot(rr))
stopifnot(is(cc <- coef(fm2), "coef.mer"),
	  is.list(cc), length(cc) == 1, is.data.frame(cc[[1]]))
print(plot(cc))

showProc.time() #

## Invalid factor specification -- used to seg.fault:
set.seed(1)
dat <- within(data.frame(lagoon = factor(rep(1:4,each = 25)),
                         habitat = factor(rep(1:20, each = 5))),
          {
              y <- round(10*rnorm(100, m = 10*as.numeric(lagoon)))
          })

tt <- suppressWarnings(try(reg <- lmer(y ~ habitat + (1|habitat*lagoon), data = dat)
                                        )
                                         # did seg.fault)
    ) # now gives error                 ^- should be ":"
## suppress warning that uses different quoting conventions on
## R-release vs. R-devel

## ignore singular fits as well as hess/grad problems
## (Windows gets singular fits, other platforms don't ...)
ctrl0 <- lmerControl(
    check.conv.singular="ignore",
    check.conv.hess="ignore",
    check.conv.grad="ignore")
r1  <- lmer(y ~ 0+habitat + (1|habitat:lagoon), data = dat,
            control=ctrl0) # ok, but senseless
r1b <- lmer(y ~ 0+habitat + (1|habitat), data = dat,
            control=ctrl0) # same model, clearly unidentifiable
## "TODO" :  summary(r1)  should ideally warn the user
stopifnot(all.equal(fixef(r1), fixef(r1b), tolerance= 1e-15),
          all.equal(ranef(r1), ranef(r1b), tolerance= 1e-15, check.attributes=FALSE))

## Use a more sensible model:
r2.0 <- lmer(y ~ 0+lagoon + (1|habitat:lagoon), data = dat) # ok
r2   <- lmer(y ~ 0+lagoon + (1|habitat), data = dat) # ok, and more clear
stopifnot(all.equal(fixef(r2), fixef(r2.0), tolerance= 1e-15),
          all.equal(ranef(r2), ranef(r2.0), tolerance= 1e-15, check.attributes=FALSE))
V2 <- vcov(r2)
assert.EQ.mat(V2, diag(x = 9.9833/3, nr = 4))
stopifnot(all.equal(unname(fixef(r2)) - (1:4)*100,
		    c(1.72, 0.28, 1.76, 0.8), tolerance = 1e-13))

## sparseX version should give same numbers:
## (only gives a warning now -- sparseX disregarded)
if(FALSE) { ## no longer
r2.  <- lmer(y ~ 0+lagoon + (1|habitat), data = dat,
             sparseX = TRUE)

## the summary() components we do want to compare 'dense X' vs 'sparse X':
nmsSumm <- c("methTitle", "devcomp", "logLik", "ngrps", "coefficients",
             "sigma", "REmat", "AICtab")
sr2  <- summary(r2)
sr2. <- summary(r2.)
sr2.$devcomp$dims['spFe'] <- 0L       # to allow for comparisons below
stopifnot(all.equal(sr2[nmsSumm], sr2.[nmsSumm], tolerance= 1e-14)
          , all.equal(ranef(r2), ranef(r2.), tolerance= 1e-14)
          , Matrix:::isDiagonal(vcov(r2.)) # ok
          , all.equal(Matrix::diag(vcov(r2.)), rep.int(V2[1,1], 4), tolerance= 1e-13)
#          , all(vcov(r2.)@factors$correlation == diag(4))  # not sure why this fails
          , TRUE)
r2.
}

### mcmcsamp() :
## From: Andrew Gelman <gelman@stat.columbia.edu>
## Date: Wed, 18 Jan 2006 22:00:53 -0500

if (FALSE) {  # mcmcsamp still needs work
    ## NB: Need to restore coda to the Suggests: field of DESCRIPTION
    ## file if this code block is reinstated.
    ## has.coda <- require(coda)
    ## if(!has.coda)
    ##     cat("'coda' package not available; some outputs will look suboptimal\n")

    ## Very simple example
    y <- 1:10
    group <- gl(2,5)
    (M1 <- lmer (y ~ 1 + (1 | group))) # works fine
    (r1 <- mcmcsamp (M1))              # dito
    r2 <- mcmcsamp (M1, saveb = TRUE)  # gave error in 0.99-* and 0.995-[12]
    (r10 <- mcmcsamp (M1, n = 10, saveb = TRUE))

    ## another one, still simple
    y <- (1:20)*pi
    x <- (1:20)^2
    group <- gl(2,10)
    M1 <- lmer (y ~ 1 | group)
    mcmcsamp (M1, n = 2, saveb=TRUE) # fine

    M2 <- lmer (y ~ 1 + x + (1 + x | group)) # false convergence
    ## should be identical (and is)
    M2 <- lmer (y ~ x + ( x | group))#  false convergence -> simulation doesn't work:
    if(FALSE) ## try(..) fails here (in R CMD check) [[why ??]]
        mcmcsamp (M2, saveb=TRUE)
    ## Error: inconsistent degrees of freedom and dimension ...

    ## mcmc for glmer:
    rG1k <- mcmcsamp(m1, n = 1000)
    summary(rG1k)
    rG2 <- mcmcsamp(m1, n = 3, verbose = TRUE)
}

## Spencer Graves' example (from a post to S-news, 2006-08-03) ----------------
## it should give an error, rather than silent non-sense:
tstDF <- data.frame(group = letters[1:5], y = 1:5)
assertError(## Now throws an error, as desired :
            lmer(y ~ 1 + (1|group), data = tstDF)
            )

showProc.time() #

## Wrong formula gave a seg.fault at times:
set.seed(2)# !
D <-  data.frame(y= rnorm(12,10), ff = gl(3,2,12),
                 x1=round(rnorm(12,3),1), x2=round(rnorm(12,7),1))
## NB: The first two are the same, having a length-3 R.E. with 3 x 3 vcov-matrix:
## --> do need CPU
## suppressWarnings() for warning about too-few random effects levels
tmpf <- function(form) lmer(form, data = D , control=lmerControl(check.conv.singular="ignore",
                                                                 check.nobs.vs.nRE="ignore",
                                                                 calc.derivs=FALSE))
m0 <- tmpf(y ~ (x1 + x2)|ff)
m1 <- tmpf(y ~ x1 + x2|ff)
m2 <- tmpf(y ~ x1 + (x2|ff))
m3 <- tmpf(y ~ (x2|ff) + x1)
suppressWarnings(stopifnot(all.equal(ranef(m0), ranef(m1), tolerance = 1e-5),
          all.equal(ranef(m2), ranef(m3), tolerance = 1e-5),
          inherits(tryCatch(lmer(y ~ x2|ff + x1, data = D), error = function(e)e),
                   "error")))

showProc.time() #

## Reordering of grouping factors should not change the internal structure
#Pm1  <- lmer1(strength ~ (1|batch) + (1|sample), Pastes, doFit = FALSE)
#Pm2  <- lmer1(strength ~ (1|sample) + (1|batch), Pastes, doFit = FALSE)
#P2.1 <- lmer (strength ~ (1|batch) + (1|sample), Pastes, devFunOnly = TRUE)
#P2.2 <- lmer (strength ~ (1|sample) + (1|batch), Pastes, devFunOnly = TRUE)

## The environments of Pm1 and Pm2 should be identical except for
## "call" and "frame":
#stopifnot(## all.EQ(env(Pm1), env(Pm2)),
#	  all.EQ(S4_2list(P2.1),
#		 S4_2list(P2.2)))


## example from Kevin Thorpe: synthesized equivalent
## http://thread.gmane.org/gmane.comp.lang.r.lme4.devel/9835

## NA issue: simpler example
d <- data.frame(y=1:60,f=factor(rep(1:6,each=10)))
d$y[2] <- NA
d$f[3:4] <- NA
lmer(y~(1|f),data=d)
glmer(y~(1|f),data=d,family=poisson)

## we originally thought that these examples should be
## estimating non-zero variances, but they shouldn't ...
## number of levels with each level of replication
levs <- c(800,300,150,100,50,50,50,20,20,5,2,2,2,2)
n <- seq_along(levs)
flevels <- seq(sum(levs))
set.seed(101)
fakedat <- data.frame(DA = factor(rep(flevels,rep(n,levs))),
                      zbmi=rnorm(sum(n*levs)))
## add NA values
fakedat[sample(nrow(fakedat),100),"zbmi"] <- NA
fakedat[sample(nrow(fakedat),100),"DA"] <- NA

m5 <- lmer(zbmi ~ (1|DA) , data = fakedat,
	   control=lmerControl(check.nobs.vs.rankZ="ignore"))
m6 <- update(m5, data=na.omit(fakedat))
stopifnot(VarCorr(m5)[["DA"]] == 0,
	  VarCorr(m6)[["DA"]] == 0)

showProc.time()
