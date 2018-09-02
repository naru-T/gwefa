\name{bw_gwnsprcomp}
\alias{bw_gwnsprcomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bandwidth selection for GW non-negative PCA
}
\description{
A function for automatic bandwidth selection to calibrate a  GW non-negative PCA
}
\usage{
bw.gwnsprcomp(data, vars, k = 2, kernel = "bisquare",
          adaptive = FALSE, p = 2, theta = 0, longlat = F, dMat,nneg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Same as GWmodel::gwpca. A Spatial*DataFrame, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp.

}

  \item{vars}{
Same as GWmodel::gwpca.cv. The number of retained components; k must be less than the number of variables.

}

  \item{k}{
Same as GWmodel::gwpca.cv. The number of retained components; k must be less than the number of variables

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

  \item{nneg}{

a logical value indicating whether the loadings should be non-negative, i.e. the PAs should be constrained to the non-negative orthant.

}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the adaptive or fixed distance bandwidth.
}
\references{
  Isabella Gollini, Binbin Lu, Martin Charlton, Christopher Brunsdon, Paul Harris (2015). GWmodel:
  An R Package for Exploring Spatial Heterogeneity Using Geographically Weighted Models. Journal
  of Statistical Software, 63(17), 1-50. URL http://www.jstatsoft.org/v63/i17/.

  Binbin Lu, Paul Harris, Martin Charlton, Christopher Brunsdon (2014). The GWmodel R package:
  further topics for exploring spatial heterogeneity using geographically weighted models.
  Geo-spatial Information Science, 17(2), 85-101. URL
  http://dx.doi.org/10.1080/10095020.2014.917453


}
\author{
N.Tsutsumida,...
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "SPATIAL" }% use one of  RShowDoc("KEYWORDS")
\keyword{ "MULTIVALIATE" }% __ONLY ONE__ keyword per line
