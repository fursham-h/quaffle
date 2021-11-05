
.posterior <- function(prob, a=1, b=1){

    successes <- ceiling(prob * 1000)
    total <- 1000
    # Adapted from triplot() in the LearnBayes package
    # Plot the prior, likelihood and posterior:
    likelihood_a = successes + 1
    likelihood_b = total - successes + 1
    posterior_a = a + successes
    posterior_b = b + total - successes
    theta = seq(0.005, 0.995, length = 5000)
    prior = stats::dbeta(theta, a, b)
    likelihood = stats::dbeta(theta, likelihood_a, likelihood_b)
    posterior  = stats::dbeta(theta, posterior_a, posterior_b)
    m = max(c(prior, likelihood, posterior))

    # Print out summary statistics for the prior, likelihood and posterior:
    calcBetaMode <- function(aa, bb) { BetaMode <- (aa - 1)/(aa + bb - 2); return(BetaMode); }
    calcBetaMean <- function(aa, bb) { BetaMean <- (aa)/(aa + bb); return(BetaMean); }
    calcBetaSd   <- function(aa, bb) { BetaSd <- sqrt((aa * bb)/(((aa + bb)^2) * (aa + bb + 1))); return(BetaSd); }
    prior_mode      <- calcBetaMode(a, b)
    likelihood_mode <- calcBetaMode(likelihood_a, likelihood_b)
    posterior_mode  <- calcBetaMode(posterior_a, posterior_b)
    prior_mean      <- calcBetaMean(a, b)
    likelihood_mean <- calcBetaMean(likelihood_a, likelihood_b)
    posterior_mean  <- calcBetaMean(posterior_a, posterior_b)
    prior_sd        <- calcBetaSd(a, b)
    likelihood_sd   <- calcBetaSd(likelihood_a, likelihood_b)
    posterior_sd    <- calcBetaSd(posterior_a, posterior_b)
    return(posterior_mean)
}




#
#
# .optimalBeta <-  function(x) {
#
#
#     x <- broom::tidy(summary(x))
#     # find the beta prior using quantile1 and quantile2
#     q.min = list(p=0.25, x=x$q1)
#     q.med = list(p=0.5, x=x$median)
#     q.max = list(p=0.75, x=x$q3)
#
#     # prior parameters using median and max, and median and min
#     prior.A = LearnBayes::beta.select(q.med,q.max)
#
#     prior.B = LearnBayes::beta.select(q.med,q.min)
#     prior.Aa = prior.A[1]
#     prior.Ab = prior.A[2]
#     prior.Ba = prior.B[1]
#     prior.Bb = prior.B[2]
#     ## find the best possible beta prior
#     ## Set a start and stop point range to find the best parameters
#     if (prior.Aa < prior.Ba) {
#         start.a = prior.Aa
#         stop.a = prior.Ba
#     } else {
#         start.a = prior.Ba
#         stop.a = prior.Aa
#     }
#     if (prior.Ab < prior.Bb) {
#         start.b = prior.Ab
#         stop.b = prior.Bb
#     } else {
#         start.b = prior.Bb
#         stop.b = prior.Ab
#     }
#     seq.a = seq(from=start.a, to=stop.a, length.out=1000)
#     seq.b = seq(from=start.b, to=stop.b, length.out=1000)
#     seq.grid = expand.grid(seq.a, seq.b)
#     prior.C.q1 = qbeta(x$q1, seq.grid[,1], seq.grid[,2])
#     prior.C.q2 = qbeta(x$median, seq.grid[,1], seq.grid[,2])
#     prior.C.q3 = qbeta(x$q3, seq.grid[,1], seq.grid[,2])
#     ## Different distance measurements, manhattan, euclidean, or otherwise.
#     ## It would be interesting to run a simulation to measure a variety of distance measurements.
#     prior.C.delta = abs(prior.C.q1 - x$q1) + abs(prior.C.q2 - x$median) + abs(prior.C.q3 - x$q3)
#     ## prior.C.delta = sqrt( (prior.C.q1 - q1p)^2 + (prior.C.q2 - q2p)^2 + (prior.C.q3 - q3p)^2 )
#     optimize.seq = cbind(seq.grid, prior.C.q1, prior.C.q2, prior.C.q3, prior.C.delta)
#     ## Minimize the delta, if the min-delta is not unique then choose the first occurence
#     best.a = optimize.seq[,1][ optimize.seq[,6]==min(optimize.seq[,6])][1]
#     best.b = optimize.seq[,2][ optimize.seq[,6]==min(optimize.seq[,6])][1]
#     return(best.a)
# }










