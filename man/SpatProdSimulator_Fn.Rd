\name{SpatProdSimulator_Fn}
\alias{SpatProdSimulator_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate data from a population following a spatial surplus production process
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SpatProdSimulator_Fn(MoveMat, SD_omega = 1, SD_epsilon = 1, SD_effort = 1, CV_obs = 1, effort_par = c(0.2, 0.5), sizepar = c(1, 0.5), Scale, Dynamical_Model, n_s, n_t, r_s, n_r, loc_r, logmeanu0, alpha, beta, km2_r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MoveMat}{
%%     ~~Describe \code{MoveMat} here~~
}
  \item{SD_omega}{
%%     ~~Describe \code{SD_omega} here~~
}
  \item{SD_epsilon}{
%%     ~~Describe \code{SD_epsilon} here~~
}
  \item{SD_effort}{
%%     ~~Describe \code{SD_effort} here~~
}
  \item{CV_obs}{
%%     ~~Describe \code{CV_obs} here~~
}
  \item{effort_par}{
%%     ~~Describe \code{effort_par} here~~
}
  \item{sizepar}{
%%     ~~Describe \code{sizepar} here~~
}
  \item{Scale}{
%%     ~~Describe \code{Scale} here~~
}
  \item{Dynamical_Model}{
%%     ~~Describe \code{Dynamical_Model} here~~
}
  \item{n_s}{
%%     ~~Describe \code{n_s} here~~
}
  \item{n_t}{
%%     ~~Describe \code{n_t} here~~
}
  \item{r_s}{
%%     ~~Describe \code{r_s} here~~
}
  \item{n_r}{
%%     ~~Describe \code{n_r} here~~
}
  \item{loc_r}{
%%     ~~Describe \code{loc_r} here~~
}
  \item{logmeanu0}{
%%     ~~Describe \code{logmeanu0} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{km2_r}{
%%     ~~Describe \code{km2_r} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (MoveMat, SD_omega = 1, SD_epsilon = 1, SD_effort = 1, 
    CV_obs = 1, effort_par = c(0.2, 0.5), sizepar = c(1, 0.5), 
    Scale, Dynamical_Model, n_s, n_t, r_s, n_r, loc_r, logmeanu0, 
    alpha, beta, km2_r) 
{
    require(RandomFields)
    RF_omega = RMgauss(var = SD_omega^2, scale = Scale)
    RF_epsilon = RMgauss(var = SD_epsilon^2, scale = Scale)
    RF_effort = RMgauss(var = SD_effort^2, scale = Scale)
    effortdens_t = rlnorm(n_t, meanlog = log(effort_par[1]) - 
        effort_par[2]^2/2, sdlog = effort_par[2])
    effortdens_r = exp(RFsimulate(model = RF_effort, x = loc_r[, 
        1], y = loc_r[, 2])@data[, 1] - SD_effort^2/2)
    effortdens_rt = outer(effortdens_r, rep(1, n_t)) * outer(rep(1, 
        n_r), effortdens_t)
    catch_rt = upred_rt = u_rt = Epsilon_rt = matrix(NA, ncol = n_t + 
        n_years_burnin, nrow = n_r)
    Omega_r = RFsimulate(model = RF_omega, x = loc_r[, 1], y = loc_r[, 
        2])@data[, 1] - SD_omega^2/2
    for (t in 1:n_t) {
        Epsilon_rt[, t] = RFsimulate(model = RF_epsilon, x = loc_r[, 
            1], y = loc_r[, 2])@data[, 1]
        if (t == 1) {
            if (Dynamical_Model == "Gompertz") 
                u_rt[, t] = km2_r * exp(logmeanu0 + Omega_r)
            if (Dynamical_Model == "Ricker") 
                u_rt[, t] = km2_r * exp(logmeanu0) + Omega_r
            catch_rt[, t] = (1 - exp(-effortdens_rt[, t])) * 
                u_rt[, t]
            u_rt[, t] = exp(-effortdens_rt[, t]) * u_rt[, t]
            u_rt[, t] = as.vector(MoveMat \%*\% u_rt[, t])
            u_rt[, t] = u_rt[, t] * exp(Epsilon_rt[, t])
        }
        if (t >= 2) {
            catch_rt[, t] = (1 - exp(-effortdens_rt[, t])) * 
                u_rt[, t - 1]
            upred_rt[, t] = exp(-effortdens_rt[, t]) * u_rt[, 
                t - 1]
            upred_rt[, t] = as.vector(MoveMat \%*\% upred_rt[, 
                t])
            if (Dynamical_Model == "Gompertz") 
                u_rt[, t] = upred_rt[, t] * exp(alpha + Omega_r - 
                  beta * log(upred_rt[, t]/km2_r) + Epsilon_rt[, 
                  t])
            if (Dynamical_Model == "Ricker") 
                u_rt[, t] = upred_rt[, t] * exp(alpha + Omega_r - 
                  beta * (upred_rt[, t]/km2_r) + Epsilon_rt[, 
                  t])
        }
    }
    DF = expand.grid(s_i = 1:n_s, t_i = 1:n_t)
    DF = cbind(DF, r_i = r_s[DF[, "s_i"]], km2_i = 1)
    DF = cbind(DF, cexp_i = u_rt[as.matrix(DF[, c("r_i", "t_i")])] * 
        DF[, "km2_i"]/km2_r[DF[, "r_i"]])
    DF = cbind(DF, c_i = rpois(n_s * n_t, lambda = DF[, "cexp_i"]))
    DF = cbind(DF, encounterprob_i = 1 - exp(-DF[, "cexp_i"]))
    logSD_obs = sqrt(log(CV_obs^2 + 1))
    DF = cbind(DF, zinflognorm_i = ifelse(DF[, "c_i"] > 0, 1, 
        0) * rlnorm(nrow(DF), meanlog = log(DF[, "cexp_i"]/DF[, 
        "encounterprob_i"]), sdlog = logSD_obs))
    DF = cbind(DF, zinfgamma_i = ifelse(DF[, "c_i"] > 0, 1, 0) * 
        rgamma(nrow(DF), shape = CV_obs^(-2), scale = DF[, "cexp_i"]/DF[, 
            "encounterprob_i"] * CV_obs^2))
    catch_t = colSums(catch_rt)
    SettingsList = list(MoveMat = MoveMat, SD_omega = SD_omega, 
        SD_epsilon = SD_epsilon, SD_effort = SD_effort, CV_obs = CV_obs, 
        effort_par = effort_par, sizepar = sizepar, Scale = Scale, 
        Dynamical_Model = Dynamical_Model, n_s = n_s, n_t = n_t, 
        r_s = r_s, n_r = n_r, loc_r = loc_r, logmeanu0 = logmeanu0, 
        alpha = alpha, beta = beta, km2_r = km2_r)
    Return = list(SettingsList = SettingsList, DF = DF, catch_t = catch_t, 
        catch_rt = catch_rt, effortdens_rt = effortdens_rt, MoveMat = MoveMat, 
        upred_rt = upred_rt, u_rt = u_rt, Epsilon_rt = Epsilon_rt, 
        Omega_r = Omega_r)
    return(Return)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
