\name{gwfa}
\alias{gwfa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geographically Weighted Factor Analysis
}
\description{
This function implements GWFA. most of description in this manual from Description from GWmodel::gwpca or psych::fa.
}
\usage{
gwfa(data, elocat, vars, bw, k=2, kernel, adaptive=TRUE, p=2, theta=0, longlat=FALSE, dMat=NULL,
                 n.obs = NA, n.iter=1, rotate="oblimin", scores="regression",
                  residuals=FALSE, SMC=TRUE, covar=FALSE,missing=FALSE,impute="median",
                  min.err = 0.001,  max.iter = 50,symmetric=TRUE, warnings=TRUE, fm="minres",
                  alpha=.1,pr=.05,oblique.scores=FALSE,np.obs=NULL,use="pairwise",cor="cor",
                  correct=.5,weight=NULL,timeout,foreach, organize=FALSE...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a Spatial*DataFrame either SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp.

}
  \item{elocat}{
Same as GWmodel::gwpcaa. two-column numeric array or Spatial*DataFrame object for providing evaluation locations, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp.

}
  \item{vars}{
the number of retained components; k must be less than the number of variables.


}
  \item{kernel}{
Same as GWmodel::gwpca. Function chosen as follows:
gaussian: wgt = exp(-.5*(vdist/bw)^2);
exponential: wgt = exp(-vdist/bw);
bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
boxcar: wgt=1 if dist < bw, wgt=0 otherwise

see help(GWmodel::gw.weight) more detail.
}
 \item{adaptive}{
 if TRUE calculate an adaptive kernel where the bandwidth corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance).

}
  \item{bw}{
bandwidth used in the weighting function, possibly calculated by bw.gwpca;fixed (distance) or adaptive bandwidth(number of nearest neighbours).

}
  \item{k}{
the number of retained components; k must be less than the number of variables

}

  \item{p}{
the power of the Minkowski distance, default is 2, i.e. the Euclidean distance.

}
  \item{theta}{
an angle in radians to rotate the coordinate system, default is 0.

}
  \item{longlat}{
if TRUE, great circle distances will be calculated.

}
  \item{dMat}{
a pre-specified distance matrix, it can be calculated by the function gw.dist .

}
  \item{n.obs}{
Number of observations used to find the correlation matrix if using a correlation matrix. Used for finding the goodness of fit statistics. Must be specified if using a correlaton matrix and finding confidence intervals.

}
  \item{fm}{
Same as psych::fa.Factoring method fm="minres" will do a minimum residual as will fm="uls". Both of these use a first derivative. fm="ols" differs very slightly from "minres" in that it minimizes the entire residual matrix using an OLS procedure but uses the empirical first derivative. This will be slower. fm="wls" will do a weighted least squares (WLS) solution, fm="gls" does a generalized weighted least squares (GLS), fm="pa" will do the principal factor solution, fm="ml" will do a maximum likelihood factor analysis. fm="minchi" will minimize the sample size weighted chi square when treating pairwise correlations with different number of subjects per pair. fm ="minrank" will do a minimum rank factor analysis. "old.min" will do minimal residual the way it was done prior to April, 2017 (see discussion below). fm="alpha" will do alpha factor analysis as described in Kaiser and Coffey (1965).


}
  \item{rotate}{
Same as psych::fa. "none", "varimax", "quartimax", "bentlerT", "equamax", "varimin", "geominT" and "bifactor" are orthogonal rotations. "Promax", "promax", "oblimin", "simplimax", "bentlerQ, "geominQ" and "biquartimin" and "cluster" are possible oblique transformations of the solution. The default is to do a oblimin transformation, although versions prior to 2009 defaulted to varimax. SPSS seems to do a Kaiser normalization before doing Promax, this is done here by the call to "promax" which does the normalization before calling Promax in GPArotation.

}
  \item{n.iter}{
Same as psych::fa. Number of bootstrap interations to do in fa or fa.poly

}
  \item{residuals}{
Same as psych::fa. Should the residual matrix be shown

}

  \item{scores}{
Same as psych::fa. The default="regression" finds factor scores using regression. Alternatives for estimating factor scores include simple regression ("Thurstone"), correlaton preserving ("tenBerge") as well as "Anderson" and "Bartlett" using the appropriate algorithms ( factor.scores). Although scores="tenBerge" is probably preferred for most solutions, it will lead to problems with some improper correlation matrices.

}
  \item{SMC}{
Same as psych::fa. Use squared multiple correlations (SMC=TRUE) or use 1 as initial communality estimate. Try using 1 if imaginary eigen values are reported. If SMC is a vector of length the number of variables, then these values are used as starting values in the case of fm='pa'.

}
  \item{covar}{
Same as psych::fa. If covar is TRUE, factor the covariance matrix, otherwise factor the correlation matrix

}
  \item{missing}{
Same as psych::fa. If scores are TRUE, and missing=TRUE, then impute missing values using either the median or the mean

}
  \item{impute}{
Same as psych::fa."median" or "mean" values are used to replace missing values

}
  \item{min.err}{
Same as psych::fa. Iterate until the change in communalities is less than min.err

}
  \item{max.iter}{
Same as psych::fa. Maximum number of iterations for convergence

}

  \item{symmetric}{
Same as psych::fa. Symmetric=TRUE forces symmetry by just looking at the lower off diagonal values.

}
  \item{warnings}{
Same as psych::fa. warnings=TRUE => warn if number of factors is too many

}
  \item{alpha}{
Same as psych::fa. Alpha level for the confidence intervals for RMSEA

}
  \item{pr}{
Same as psych::fa. If doing iterations to find confidence intervals, what probability values should be found for the confidence intervals

}
  \item{oblique.scores}{
Same as psych::fa. When factor scores are found, should they be based on the structure matrix (default) or the pattern matrix (oblique.scores=TRUE). Now it is always false. If you want oblique factor scores, use tenBerge.

}
  \item{weight}{
Same as psych::fa. If not NULL, a vector of length n.obs that contains weights for each observation. The NULL case is equivalent to all cases being weighted 1.

}
  \item{np.obs}{
Same as psych::fa. The pairwise number of observations. Used if using a correlation matrix and asking for a minchi solution.

}
  \item{use}{
Same as psych::fa. How to treat missing data, use="pairwise" is the default". See cor for other options.

}
  \item{cor}{
Same as psych::fa. How to find the correlations: "cor" is Pearson", "cov" is covariance, "tet" is tetrachoric, "poly" is polychoric, "mixed" uses mixed cor for a mixture of tetrachorics, polychorics, Pearsons, biserials, and polyserials, Yuleb is Yulebonett, Yuleq and YuleY are the obvious Yule coefficients as appropriate

}
  \item{correct}{
Same as psych::fa. When doing tetrachoric, polycoric, or mixed cor, how should we treat empty cells. (See the discussion in the help for tetrachoric.)

}
  \item{timeout}{
A numeric specifying the maximum number of seconds the expression is allowed to run before being interrupted by the timeout. (See ?R.utils::wituTimeout)

}
  \item{foreach}{
default:FALSE. If TRUE, foreach function works to implement calculation using multicores.

}

  \item{organize}{
default:FALSE.if TRUE, ordering factor loadings and scores are based on the similarity of conventional FA. mean absolute biases (MAB) between loadings from local fa and conventional fa are calculated and the minimum MAB is used for the matching. if FALSE, pre-rotated order is used for ordering factors and associated outputs are reordered as factor names. see: https://stats.stackexchange.com/a/169991

}

  \item{...}{
Same as psych::fa. Additional parameters, specifically, keys may be passed if using the target rotation, or delta if using geominQ, or whether to normalize if using Varimax

}

}
\details{
This function implements a explanatory factor analysis, via Geographically weighted approach.

The following detail is from psych::fa
Factor analysis is an attempt to approximate a correlation or covariance matrix with one of lesser rank. The basic model is that nRn = nFk kFn' + U2 where k is much less than n. There are many ways to do factor analysis, and maximum likelihood procedures are probably the most commonly preferred (see factanal ). The existence of uniquenesses is what distinguishes factor analysis from principal components analysis (e.g., principal). If variables are thought to represent a “true" or latent part then factor analysis provides an estimate of the correlations with the latent factor(s) representing the data. If variables are thought to be measured without error, then principal components provides the most parsimonious description of the data.

Factor loadings will be smaller than component loadings for the later reflect unique error in each variable. The off diagonal residuals for a factor solution will be superior (smaller) that of a component model. Factor loadings can be thought of as the asymptotic component loadings as the number of variables loading on each factor increases.

The fa function will do factor analyses using one of six different algorithms: minimum residual (minres, aka ols, uls), principal axes, alpha factoring, weighted least squares, minimum rank, or maximum likelihood.

Principal axes factor analysis has a long history in exploratory analysis and is a straightforward procedure. Successive eigen value decompositions are done on a correlation matrix with the diagonal replaced with diag (FF') until sum(diag(FF')) does not change (very much). The current limit of max.iter =50 seems to work for most problems, but the Holzinger-Harmon 24 variable problem needs about 203 iterations to converge for a 5 factor solution.

Not all factor programs that do principal axes do iterative solutions. The example from the SAS manual (Chapter 33) is such a case. To achieve that solution, it is necessary to specify that the max.iter = 1. Comparing that solution to an iterated one (the default) shows that iterations improve the solution.

In addition, fm="mle" produces even better solutions for this example. Both the RMSEA and the root mean square of the residuals are smaller than the fm="pa" solution.

However, simulations of multiple problem sets suggest that fm="pa" tends to produce slightly smaller residuals while having slightly larger RMSEAs than does fm="minres" or fm="mle". That is, the sum of squared residuals for fm="pa" seem to be slightly smaller than those found using fm="minres" but the RMSEAs are slightly worse when using fm="pa". That is to say, the "true" minimal residual is probably found by fm="pa".

Following extensive correspondence with Hao Wu and Mikko Ronkko, in April, 2017 the derivative of the minres and uls) fitting was modified. This leads to slightly smaller residuals (appropriately enough for a method claiming to minimize them) than the prior procedure. For consistency with prior analyses, "old.min" was added to give these slightly larger residuals. The differences between old.min and the newer "minres" and "ols" solutions are at the third to fourth decimal, but none the less, are worth noting. For comparison purposes, the fm="ols" uses empirical first derivatives, while uls and minres use equation based first derivatives. The results seem to be identical, but the minres and uls solutions require fewer iterations for larger problems and are faster. Thanks to Hao Wu for some very thoughtful help.

Although usually these various algorithms produce equivalent results, there are several data sets included that show large differences between the methods. Schutz produces Heywood and super Heywood cases, blant leads to very different solutions.

Principal axes may be used in cases when maximum likelihood solutions fail to converge, although fm="minres" will also do that and tends to produce better (smaller RMSEA) solutions.

The fm="minchi" option is a variation on the "minres" (ols) solution and minimizes the sample size weighted residuals rather than just the residuals. This was developed to handle the problem of data that Massively Missing Completely at Random (MMCAR) which a condition that happens in the SAPA project.

A problem in factor analysis is to find the best estimate of the original communalities. Using the Squared Multiple Correlation (SMC) for each variable will underestimate the communalities, using 1s will over estimate. By default, the SMC estimate is used. In either case, iterative techniques will tend to converge on a stable solution. If, however, a solution fails to be achieved, it is useful to try again using ones (SMC =FALSE). Alternatively, a vector of starting values for the communalities may be specified by the SMC option.

The iterated principal axes algorithm does not attempt to find the best (as defined by a maximum likelihood criterion) solution, but rather one that converges rapidly using successive eigen value decompositions. The maximum likelihood criterion of fit and the associated chi square value are reported, and will be (slightly) worse than that found using maximum likelihood procedures.

The minimum residual (minres) solution is an unweighted least squares solution that takes a slightly different approach. It uses the optim function and adjusts the diagonal elements of the correlation matrix to mimimize the squared residual when the factor model is the eigen value decomposition of the reduced matrix. MINRES and PA will both work when ML will not, for they can be used when the matrix is singular. Although before the change in the derivative, the MINRES solution was slightly more similar to the ML solution than is the PA solution. With the change in the derivative of the minres fit, the minres, pa and uls solutions are practically identical. To a great extent, the minres and wls solutions follow ideas in the factanal function with the change in the derivative.

The weighted least squares (wls) solution weights the residual matrix by 1/ diagonal of the inverse of the correlation matrix. This has the effect of weighting items with low communalities more than those with high communalities.

The generalized least squares (gls) solution weights the residual matrix by the inverse of the correlation matrix. This has the effect of weighting those variables with low communalities even more than those with high communalities.

The maximum likelihood solution takes yet another approach and finds those communality values that minimize the chi square goodness of fit test. The fm="ml" option provides a maximum likelihood solution following the procedures used in factanal but does not provide all the extra features of that function. It does, however, produce more expansive output.

The minimum rank factor model (MRFA) roughly follows ideas by Shapiro and Ten Berge (2002) and Ten Berge and Kiers (1991). It makes use of the glb.algebraic procedure contributed by Andreas Moltner. MRFA attempts to extract factors such that the residual matrix is still positive semi-definite. This version is still being tested and feedback is most welcome.

Alpha factor analysis finds solutions based upon a correlation matrix corrected for communalities and then rescales these to the original correlation matrix. This procedure is described by Kaiser and Coffey, 1965.

Test cases comparing the output to SPSS suggest that the PA algorithm matches what SPSS calls uls, and that the wls solutions are equivalent in their fits. The wls and gls solutions have slightly larger eigen values, but slightly worse fits of the off diagonal residuals than do the minres or maximum likelihood solutions. Comparing the results to the examples in Harman (1976), the PA solution with no iterations matches what Harman calls Principal Axes (as does SAS), while the iterated PA solution matches his minres solution. The minres solution found in psych tends to have slightly smaller off diagonal residuals (as it should) than does the iterated PA solution.

Although for items, it is typical to find factor scores by scoring the salient items (using, e.g., scoreItems) factor scores can be estimated by regression as well as several other means. There are multiple approaches that are possible (see Grice, 2001) and one taken here was developed by tenBerge et al. (see factor.scores). The alternative, which will match factanal is to find the scores using regression – Thurstone's least squares regression where the weights are found by W = inverse(R)S where R is the correlation matrix of the variables ans S is the structure matrix. Then, factor scores are just Fs = X W.

In the oblique case, the factor loadings are referred to as Pattern coefficients and are related to the Structure coefficients by S = P Phi and thus P = S Phi^-1. When estimating factor scores, fa and factanal differ in that fa finds the factors from the Structure matrix while factanal seems to do it from the Pattern matrix. Thus, although in the orthogonal case, fa and factanal agree perfectly in their factor score estimates, they do not agree in the case of oblique factors. Setting oblique.scores = TRUE will produce factor score estimate that match those of factanal.

It is sometimes useful to extend the factor solution to variables that were not factored. This may be done using fa.extension. Factor extension is typically done in the case where some variables were not appropriate to factor, but factor loadings on the original factors are still desired.

For dichotomous items or polytomous items, it is recommended to analyze the tetrachoric or polychoric correlations rather than the Pearson correlations. This may be done by specifying cor="poly" or cor="tet" or cor="mixed" if the data have a mixture of dichotomous, polytomous, and continous variables.

Analysis of dichotomous or polytomous data may also be done by using irt.fa or simply setting the cor="poly" option. In the first case, the factor analysis results are reported in Item Response Theory (IRT) terms, although the original factor solution is returned in the results. In the later case, a typical factor loadings matrix is returned, but the tetrachoric/polychoric correlation matrix and item statistics are saved for reanalysis by irt.fa. (See also the mixed.cor function to find correlations from a mixture of continuous, dichotomous, and polytomous items.)

Of the various rotation/transformation options, varimax, Varimax, quartimax, bentlerT, geominT, and bifactor do orthogonal rotations. Promax transforms obliquely with a target matix equal to the varimax solution. oblimin, quartimin, simplimax, bentlerQ, geominQ and biquartimin are oblique transformations. Most of these are just calls to the GPArotation package. The “cluster” option does a targeted rotation to a structure defined by the cluster representation of a varimax solution. With the optional "keys" parameter, the "target" option will rotate to a target supplied as a keys matrix. (See target.rot.)

Two additional target rotation options are available through calls to GPArotation. These are the targetQ (oblique) and targetT (orthogonal) target rotations of Michael Browne. See target.rot for more documentation.

The "bifactor" rotation implements the Jennrich and Bentler (2011) bifactor rotation by calling the GPForth function in the GPArotation package and using two functions adapted from the MatLab code of Jennrich and Bentler. This seems to have a problem with local minima and multiple starting values should be used.

There are two varimax rotation functions. One, Varimax, in the GPArotation package does not by default apply Kaiser normalization. The other, varimax, in the stats package, does. It appears that the two rotation functions produce slightly different results even when normalization is set. For consistency with the other rotation functions, Varimax is probably preferred.

The rotation matrix (rot.mat) is returned from all of these options. This is the inverse of the Th (theta?) object returned by the GPArotation package. The correlations of the factors may be found by Phi = Th' Th

There are two ways to handle dichotomous or polytomous responses: fa with the cor="poly" option which will return the tetrachoric or polychoric correlation matrix, as well as the normal factor analysis output, and irt.fa which returns a two parameter irt analysis as well as the normal fa output.

When factor analyzing items with dichotomous or polytomous responses, the irt.fa function provides an Item Response Theory representation of the factor output. The factor analysis results are available, however, as an object in the irt.fa output.

fa.poly is deprecated, for its functioning is matched by setting cor="poly". It will produce normal factor analysis output but also will save the polychoric matrix (rho) and items difficulties (tau) for subsequent irt analyses. fa.poly will, by default, find factor scores if the data are available. The correlations are found using either tetrachoric or polychoric and then this matrix is factored. Weights from the factors are then applied to the original data to estimate factor scores.

The function fa will repeat the analysis n.iter times on a bootstrapped sample of the data (if they exist) or of a simulated data set based upon the observed correlation matrix. The mean estimate and standard deviation of the estimate are returned and will print the original factor analysis as well as the alpha level confidence intervals for the estimated coefficients. The bootstrapped solutions are rotated towards the original solution using target.rot. The factor loadings are z-transformed, averaged and then back transformed. This leads to an error in the case of Heywood cases. The probably better alternative is to just find the mean bootstrapped value and find the confidence intervals based upon the observed range of values. The default is to have n.iter =1 and thus not do bootstrapping.

If using polytomous or dichotomous items, it is perhaps more useful to find the Item Response Theory parameters equivalent to the factor loadings reported in fa.poly by using the irt.fa function.

Some correlation matrices that arise from using pairwise deletion or from tetrachoric or polychoric matrices will not be proper. That is, they will not be positive semi-definite (all eigen values >= 0). The cor.smooth function will adjust correlation matrices (smooth them) by making all negative eigen values slightly greater than 0, rescaling the other eigen values to sum to the number of variables, and then recreating the correlation matrix. See cor.smooth for an example of this problem using the burt data set.

One reason for this problem when using tetrachorics or polychorics seems to be the adjustment for continuity. Setting correct=0 turns this off and seems to produce more proper matrices.

For those who like SPSS type output, the measure of factoring adequacy known as the Kaiser-Meyer-Olkin KMO test may be found from the correlation matrix or data matrix using the KMO function. Similarly, the Bartlett's test of Sphericity may be found using the cortest.bartlett function.

For those who want to have an object of the variances accounted for, this is returned invisibly by the print function. (e.g., p <- print(fa(ability))$Vaccounted ). This is now returned by the fa function as well (e.g. p <- fa(ability)$Vaccounted ).

The output from the print.psych.fa function displays the factor loadings (from the pattern matrix, the h2 (communalities) the u2 (the uniquenesses), com (the complexity of the factor loadings for that variable (see below). In the case of an orthogonal solution, h2 is merely the row sum of the squared factor loadings. But for an oblique solution, it is the row sum of the orthogonal factor loadings (remember, that rotations or transformations do not change the communality).

fa.sapa will do iterative solutions for successive random samples of fractions (frac) of the data set. This allows us to find the stability of solutions for various sample sizes.
}
\value{
  \item{loadings }{Arrays of factor loadings at each location.}
  \item{score }{A matrix of an estimated factor score at each location.}
  \item{uniquenesses }{A matrix of uniquenesses at each location.}
  \item{loading_score }{A matrix of loading scores at each location.}
  \item{ss }{A matrix of proportion variances at each location.}
  \item{cor.mt }{A matrix of correlation coefficients at each location.}
  \item{rmsea }{A vector of RMSEA.}
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
