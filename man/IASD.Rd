\name{IASD}
\alias{IASD}
\title{
Model selection for index of asymmetry distribution
}
\description{
Calculate AIC's and AICc's of unimodal model (one normal distribution) and bimodal model(a mixture of two normal distributions) which fit the distribution of indices of asymmetry (IAS), and plot their density, to help determine IAS distribution is unimodel or bimodal.
}
\usage{
IASD(df, cols = NA, fixSignApproximation = FALSE, plotGraph = TRUE, 
     plotToScreen = FALSE, filePrefix = "Graph", xlimMin = NA, xlimMax = NA,  
     ylimMin = 0, ylimMax = NA, dHist = NA, dFunc = NA, meanStartSymmetric = NA, 
     sdStartSymmetric = NA, meanStartAsymmetric = NA, sdStartAsymmetric = NA, 
     positiveRatioStartAsymmetric = NA)
}
\arguments{
  \item{df}{
data frame containg the data to be investigated.
}
  \item{cols}{
Columns in df to be processed.  If NA, they are from the second to the last columns for multi-column data frame and the sole column for single-column data frame.
}
  \item{fixSignApproximation}{
If TRUE, parameters of normal distributions are determined only by using absolute values, rather than MLE.  Each column can be separately controled by using vector value.
}
  \item{plotGraph}{
If TRUE, histograms and density plots are plotted and saved to PDF, if FALSE, they are not plotted.
}
  \item{plotToScreen}{
If TRUE, plotted graphs are also shown in the screen.
}
  \item{filePrefix}{
File names of each saved plot start with this value.
}
  \item{xlimMin}{
Minimum of plot range.  If NA, it is determined from the data.  Each column can be separately controled by using vector value.
}
  \item{xlimMax}{
Maximum of plot range.  If NA, it is determined from the data.  Each column can be separately controled by using vector value.
}
  \item{dHist}{
Width of histogram bars.  If NA, it is one 20th of the plot range.  Each column can be separately controled by using vector value.
}
  \item{dFunc}{
Broken line step for the plot of density functions.  If NA, it is one 200th of the plot range.  Each column can be separately controled by using vector value.
}
  \item{ylimMin}{
Minimum of vertical axis of plots.  If NA, it is determined by "hist()" function of R.  Each column can be separately controled by using vector value.
}
  \item{ylimMax}{
Maximum of vertical axis of plots.  If NA, it is determined by "hist()" function of R.  Each column can be separately controled by using vector value.
}
  \item{meanStartSymmetric}{
Start value of mean for mle() in bimodal symmetric model.  If NA, it is calculated by using ablolute values of the data.  Each column can be separately controled by using vector value.
}
  \item{sdStartSymmetric}{
Start value of sd for mle() in bimodal symmetric model.  If NA, it is calculated by using ablolute values of the data.  Each column can be separately controled by using vector value.
}
  \item{meanStartAsymmetric}{
Start value of mean for mle() in bimodal asymmetric model.  If NA, it is calculated by using ablolute values of the data.  Each column can be separately controled by using vector value.
}
  \item{sdStartAsymmetric}{
Start value of sd for mle() in bimodal asymmetric model.  If NA, it is calculated by using ablolute values of the data.  Each column can be separately controled by using vector value.
}
  \item{positiveRatioStartAsymmetric}{
Start value of positiveRatio for mle() in bimodal asymmetric model.  If NA, it is ratio of positive data.  Each column can be separately controled by using vector value.
}
}
\details{
Calculate AIC and AICc for the following four models and plot their densities.
\enumerate{
\item unimodel symmetric distribution (normal distribution with mean = 0)\cr
    N(0, sd^2) 
\item unimodel asymmetric distribution (normal distribution)\cr
    N(mean, sd^2)
\item bimodel symmetric distribution (mixture of two normal distributions with opposite sign of mean but same absolute values and weights)\cr
    0.5*N(mean, sd^2) + 0.5*N(- mean, sd^2)
\item bimodel asymmetric distribution (weighted mixture of two normal distributions with opposite sign of mean and the same absolute values)\cr
    positiveRatio*N(mean, sd^2) + (1 - positiveRatio)*N(- mean, sd^2)
 	}
 	
 	If the start vlues for mle() (meanStartSymmetric, sdStartSymmetric, meanStartAsymmetric, sdStartAsymmetric, positiveRatioStartAsymmetric) are inappropriate values, mle() does not work properly.  If they are not assigned (NA), mean and sd are those of absolute values of the data, and positiveRatio is the ratio of positive data.
}
\value{
  \item{colname}{list for the column with colname in the argument df}
  \item{UnimodalSymmetric}{list for the unimodal symmetric model}
  \item{UnimodalAsymmetric}{list for the unimodal asymmetric model}
  \item{BimodalSymmetric}{list for the bimodal symmetric model}
  \item{BimodalAsymmetric}{list for the bimodal asymmetric model}
  \item{AIC}{AIC (Akaike's information riterion)}
  \item{AICc}{AICc (AIC with a correction for finite sample sizes)}
  \item{mean}{estimated value of mean}
  \item{sd}{estimated value of sd}
  \item{positiveRatio}{estimated value of positiveRatio}
  \item{f}{density function}
}
%\references{
%}
\author{
Satoshi Takahashi
}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
\examples{
	df = data.frame(ID = c(1:5), IAS = c(8.3, 12.7, -12.7, -7.3, -8.1),
	   IAS2 = c(14.2, 8.8, -12.7, -8.6, -10.5),
	   IAS3 = c(1.04, 1.28, -0.78, -0.84, -0.85))
	result = IASD(df) # calculate AIC's and AICc's
	result = IASD(df, cols = c(2,4), plotGraph = FALSE) 
		# use data in the second and fourth columns, do not plot graphs
	result = IASD(df, filePrefix="P.microlepis", xlimMin = -15, 
	xlimMax = 15, dHist = c(1, 1, 0.1))  
		# file name of each plot starts with "P.microlepis", plot range 
		# and width of histgram bar is changed
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ AIC }
\keyword{ IAS }
