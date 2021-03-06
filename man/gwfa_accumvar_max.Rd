\name{gwfa.Accumvar_max.calc}
\alias{gwfa.Accumvar_max.calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate Detection of the maximum of the minimum accumulated variation for a specified bandwidth for GWFA

}
\description{
This function finds the maximum of the minimum accumulated variation for a specified bandwidth for GWFA. It can be used to construct the bandwidth function across all possible bandwidths and compared to that found automatically via bw_gwfa function.


}
\usage{
gwfa.Accumvar_max.calc(bw, x, dp.locat, k, robust, kernel, adaptive,
p, theta, longlat, dMat, vars, n.obs, fm, rotate, scores, oblique.scores, timeout, foreach)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bw}{
bandwidth used in the weighting function;fixed (distance) or adaptive bandwidth(number of nearest neighbours)

Description from GWmodel::gwpca.cv
}
  \item{x}{
Same as GWmodel::gwpca. A Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp.

}
  \item{dp.locat}{
Same as GWmodel::gwpca.cv. A two-column numeric array of observation coordinates

}
  \item{k}{
Same as GWmodel::gwpca.cv. The number of retained components; k must be less than the number of variables

}
  \item{robust}{
Same as GWmodel::gwpca.cv. If TRUE, robust GWPCA will be applied; otherwise basic GWPCA will be applied

}
  \item{kernel}{
Same as GWmodel::gwpca.cv. Function chosen as follows:
gaussian: wgt = exp(-.5*(vdist/bw)^2);
exponential: wgt = exp(-vdist/bw);
bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
boxcar: wgt=1 if dist < bw, wgt=0 otherwise

}
  \item{adaptive}{
Same as GWmodel::gwpca.cv. If TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance)

}
  \item{p}{
Same as GWmodel::gwpca. The power of the Minkowski distance, default is 2, i.e. the Euclidean distance.

}
  \item{theta}{
Same as GWmodel::gwpca. An angle in radians to rotate the coordinate system, default is 0.

}
  \item{longlat}{
Same as GWmodel::gwpca. If TRUE, great circle distances will be calculated.

}
  \item{dMat}{
Same as GWmodel::gwpca. A pre-specified distance matrix, it can be calculated by the function gw.dist .

}
  \item{vars}{
Same as GWmodel::gwpca. The number of retained components; k must be less than the number of variables.

}
  \item{n.obs}{
Same as psych::fa. Number of observations used to find the correlation matrix if using a correlation matrix. Used for finding the goodness of fit statistics. Must be specified if using a correlaton matrix and finding confidence intervals.

}
  \item{fm}{
Same as psych::fa. Factoring method fm="minres" will do a minimum residual as will fm="uls". Both of these use a first derivative. fm="ols" differs very slightly from "minres" in that it minimizes the entire residual matrix using an OLS procedure but uses the empirical first derivative. This will be slower. fm="wls" will do a weighted least squares (WLS) solution, fm="gls" does a generalized weighted least squares (GLS), fm="pa" will do the principal factor solution, fm="ml" will do a maximum likelihood factor analysis. fm="minchi" will minimize the sample size weighted chi square when treating pairwise correlations with different number of subjects per pair. fm ="minrank" will do a minimum rank factor analysis. "old.min" will do minimal residual the way it was done prior to April, 2017 (see discussion below). fm="alpha" will do alpha factor analysis as described in Kaiser and Coffey (1965).

}
  \item{rotate}{
Same as psych::fa. "none", "varimax", "quartimax", "bentlerT", "equamax", "varimin", "geominT" and "bifactor" are orthogonal rotations. "Promax", "promax", "oblimin", "simplimax", "bentlerQ, "geominQ" and "biquartimin" and "cluster" are possible oblique transformations of the solution. The default is to do a oblimin transformation, although versions prior to 2009 defaulted to varimax. SPSS seems to do a Kaiser normalization before doing Promax, this is done here by the call to "promax" which does the normalization before calling Promax in GPArotation.

}

  \item{scores}{
the default="regression" finds factor scores using regression. Alternatives for estimating factor scores include simple regression ("Thurstone"), correlaton preserving ("tenBerge") as well as "Anderson" and "Bartlett" using the appropriate algorithms ( factor.scores). Although scores="tenBerge" is probably preferred for most solutions, it will lead to problems with some improper correlation matrices.
}
\item{oblique.scores}{
 When factor scores are found, should they be based on the structure matrix (default) or the pattern matrix (oblique.scores=TRUE). Now it is always false. If you want oblique factor scores, use tenBerge. (See ?psych::fa)
 }
  \item{timeout}{
A numeric specifying the maximum number of seconds the expression is allowed to run before being interrupted by the timeout. (See ?R.utils::wituTimeout)
}
  \item{foreach}{
default:FALSE. If TRUE, foreach function works to implement calculation using multicores.

}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the maximum of the minimum accumulated variation.
}
\references{
  Isabella Gollini, Binbin Lu, Martin Charlton, Christopher Brunsdon, Paul Harris (2015). GWmodel:
  An R Package for Exploring Spatial Heterogeneity Using Geographically Weighted Models. Journal
  of Statistical Software, 63(17), 1-50. URL http://www.jstatsoft.org/v63/i17/.

  Binbin Lu, Paul Harris, Martin Charlton, Christopher Brunsdon (2014). The GWmodel R package:
  further topics for exploring spatial heterogeneity using geographically weighted models.
  Geo-spatial Information Science, 17(2), 85-101. URL
  http://dx.doi.org/10.1080/10095020.2014.917453

  Revelle, W. (2017) psych: Procedures for Personality and Psychological Research, Northwestern
  University, Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.7.8.

  Tsutsumida N., Harris P., Comber A. (2017) The application of a geographically weighted principal components analysis for exploring 23 years of goat population change across Mongolia, Annals of the American Association of Geographers, 107(5) 1060-1074; doi:  http://dx.doi.org/10.1080/24694452.2017.1309968
}
\author{
N. Tsutsumida,...
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%% TBD
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "SPATIAL" }% use one of  RShowDoc("KEYWORDS")
\keyword{ "MULTIVALIATE" }% __ONLY ONE__ keyword per line
