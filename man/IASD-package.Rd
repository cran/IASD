\name{IASD-package}
\alias{IASD-package}
\docType{package}
\title{
Model Selection for Index of Asymmetry Distribution
}
\description{
Calculate AIC's and AICc's of unimodal model (one normal distribution) and bimodal model(a mixture of two normal distributions) which fit the distribution of indices of asymmetry (IAS), and plot their density, to help determine IAS distribution is unimodal or bimodal.
}
\details{
\tabular{ll}{
Package: \tab IASD\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2023-09-01\cr
License: \tab  GPL (>= 2)\cr
}
IASD(df) calculates AIC's and AICc's of unimodal model and bimodal model for the distribution of indices of asymmetry in the data frame df, and plots density functions.
}
\author{
Satoshi Takahashi

Maintainer: Satoshi Takahashi <takahasi@lisboa.ics.nara-wu.ac.jp>
}
\keyword{ package }
\examples{
	df = data.frame(ID = c(1:5), IAS = c(8.3, 12.7, -12.7, -7.3, -8.1),
	   IAS2 = c(14.2, 8.8, -12.7, -8.6, -10.5),
	   IAS3 = c(1.04, 1.28, -0.78, -0.84, -0.85))
		# Do not use the option 'generateFiles = FALSE', in the following IASD commands.
	result = IASD(df, generateFiles = FALSE) # calculate AIC's and AICc's
	result = IASD(df, cols = c(2,4), plotGraph = FALSE, generateFiles = FALSE) 
		# use data in the second and fourth columns, do not plot graphs
	result = IASD(df, filePrefix="P.microlepis", xlimMin = -15, 
		xlimMax = 15, dHist = c(1, 1, 0.1), generateFiles = FALSE)  
		# file name of each plot starts with "P.microlepis", plot range 
		# and width of histgram bar is changed
}
