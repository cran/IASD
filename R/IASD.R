IASD <-
function(df, cols = NA, fixSignApproximation = FALSE, plotGraph = TRUE, 
		plotToScreen = FALSE, filePrefix = "Graph", xlimMin = NA, xlimMax = NA, 
		ylimMin = 0, ylimMax = NA, dHist = NA, dFunc = NA, meanStartSymmetric = NA, 
		sdStartSymmetric = NA, meanStartAsymmetric = NA, sdStartAsymmetric = NA, 
		positiveRatioStartAsymmetric = NA){  
	# df: data frame
	# range: range of the distribution graph, 
	# dHist: division width of histgram
	# dFunc: division of function graph
	
	#Sample Variance. var in R is unbiased variance
	BVar <- function(x){
		sum((x-mean(x))*(x-mean(x)))/length(x)
	}
	
	#used to calculated AICc, the small-sample-size corrected version of AIC
	AICC <- function(p, n){
		p*n/(n-p-1)
	}

	PlotHist <- function(ylimMin, ylimMax){
		if (is.na(ylimMin) || is.na(ylimMax))
			hist(ias, freq=F, breaks=breaks, xlab=colname, 
				main=sprintf("Histogram of %s", colname))
		else
			hist(ias, freq=F, breaks=breaks, ylim = c(ylimMin, ylimMax), xlab=colname, 
				main=sprintf("Histogram of %s", colname))
	}
		
	if (is.na(cols)){
		ncdf <- ncol(df)
		if (ncdf >1)
			cols = c(2:ncdf)
		else
			cols <- 1
	}
	nc <- length(cols)
	sd1 <- rep(0, nc)
	sd2 <- rep(0, nc)
	sd3 <- rep(0, nc)
	sd4 <- rep(0, nc)
	mean2 <- rep(0, nc)
	mean3 <- rep(0,nc)
	mean4 <- rep(0, nc)
	positiveRatio4 <- rep(0, nc)
	UniSymAIC <- rep(0, nc)
	UniSymAICC <- rep(0, nc)
	UniAsymAIC <- rep(0, nc)
	UniAsymAICC <- rep(0, nc)
	BiSymAIC <- rep(0, nc)
	BiSymAICC <- rep(0, nc)
	BiAsymAIC <- rep(0, nc)
	BiAsymAICC <- rep(0, nc)
	UnimodalSymmetric <- as.list(NA, nc)
	UnimodalAsymmetric <- as.list(NA, nc)
	BimodalSymmetric <- as.list(NA, nc)
	BimodalAsymmetric <- as.list(NA, nc)
	
	RepVar <- function(x){
		if (length(x) >= nc)
			return(x)
		else
			return(rep(x, nc))
	}
	if (nc > 1){
		fixSignApproximation <- RepVar(fixSignApproximation)
		xlimMin <- RepVar(xlimMin)
		xlimMax <- RepVar(xlimMax)
		dHist <- RepVar(dHist)
		dFunc <- RepVar(dFunc)
		ylimMin <- RepVar(ylimMin)
		ylimMax <- RepVar(ylimMax)
		meanStartSymmetric <- RepVar(meanStartSymmetric)
		sdStartSymmetric <- RepVar(sdStartSymmetric)
		meanStartAsymmetric <- RepVar(meanStartAsymmetric)
		sdStartAsymmetric <- RepVar(sdStartAsymmetric)
		positiveRatioStartAsymmetric <- RepVar(positiveRatioStartAsymmetric)
	}
	i <- 0
	for(col in cols){
		i <- i+1
		colname <- colnames(df)[col]
		ias <- df[[col]]
		n <- length(ias)
		absMax = max(abs(ias))
		if (absMax > 6)
			absMax <- ceiling(absMax/5)*5
		else
			absMax <- ceiling(absMax)
		if (is.na(xlimMin[i]))
			xlimMin[i] <- - absMax
		if (is.na(xlimMax[i]))
			xlimMax[i] <- absMax
		if (is.na(dHist[i]))
			dHist[i] <- (xlimMax[i] - xlimMin[i])/20
		if (is.na(dFunc[i]))
			dFunc[i] <- (xlimMax[i] - xlimMin[i])/200
		if (plotGraph){
			breaks <- seq(xlimMin[i], xlimMax[i], dHist[i])
			iv <- seq(xlimMin[i], xlimMax[i], dFunc[i])
		}
		sd <- sqrt(sum(ias*ias)/n)
		df$PUniSym <- dnorm(ias, sd = sd)
		UniSymAIC[i] <- -2*sum(log(df$PUniSym)) + 2*1
		UniSymAICC[i] <- -2*sum(log(df$PUniSym)) + 2*AICC(1,n)
		sd1[i] <- sd
		switch(i,
			f <- function(x){dnorm(x, sd = sd1[1])},
			f <- function(x){dnorm(x, sd = sd1[2])},
			f <- function(x){dnorm(x, sd = sd1[3])},
			f <- function(x){dnorm(x, sd = sd1[4])},
			f <- function(x){dnorm(x, sd = sd1[5])},
			f <- function(x){dnorm(x, sd = sd1[6])},
			f <- function(x){dnorm(x, sd = sd1[7])},
			f <- function(x){dnorm(x, sd = sd1[8])},
			f <- function(x){dnorm(x, sd = sd1[9])},
			f <- function(x){dnorm(x, sd = sd1[10])},
		)
		UnimodalSymmetric[[i]] <- list(AIC = UniSymAIC[i], AICc = UniSymAICC[i], 
			sd = sd, f = f)
		names(UnimodalSymmetric)[i] = colname
		if (plotGraph){
			fileName = sprintf("%s-%s-Unimodal-Symmetric.pdf", 
				filePrefix, colname)
			if (!plotToScreen)
				pdf(file = fileName)
			PlotHist(ylimMin[i], ylimMax[i])
			lines(iv, dnorm(iv, sd = sd))
			if (plotToScreen)			
				dev.print(pdf,file = fileName)
			else
				dev.off()
		}
			
		mean <- mean(ias)
		sd <- sqrt(BVar(ias))
		df$PUniAsymB <- dnorm(ias, mean = mean, sd = sd)
		UniAsymAIC[i] <- -2*sum(log(df$PUniAsymB)) + 2*2
		UniAsymAICC[i] <- -2*sum(log(df$PUniAsymB)) + 2*AICC(2,n)	
		mean2[i] <- mean
		sd2[i] <- sd
		switch(i,
			f <- function(x){dnorm(x, mean = mean2[1], sd = sd2[1])},
			f <- function(x){dnorm(x, mean = mean2[2], sd = sd2[2])},
			f <- function(x){dnorm(x, mean = mean2[3], sd = sd2[3])},
			f <- function(x){dnorm(x, mean = mean2[4], sd = sd2[4])},
			f <- function(x){dnorm(x, mean = mean2[5], sd = sd2[5])},
			f <- function(x){dnorm(x, mean = mean2[6], sd = sd2[6])},
			f <- function(x){dnorm(x, mean = mean2[7], sd = sd2[7])},
			f <- function(x){dnorm(x, mean = mean2[8], sd = sd2[8])},
			f <- function(x){dnorm(x, mean = mean2[9], sd = sd2[9])},
			f <- function(x){dnorm(x, mean = mean2[10], sd = sd2[10])},
		)
		UnimodalAsymmetric[[i]] = list(AIC = UniAsymAIC[i], AICc = UniAsymAICC[i], 
			mean = mean, sd = sd, f = f)
		names(UnimodalAsymmetric)[i] = colname
		if (plotGraph){
			fileName = sprintf("%s-%s-Unimodal-Asymmetric.pdf", 
				filePrefix, colname)
			if (!plotToScreen)
				pdf(file = fileName)
			PlotHist(ylimMin[i], ylimMax[i])
			lines(iv, dnorm(iv, mean = mean, sd = sd))
			if (plotToScreen)			
				dev.print(pdf,file = fileName)
			else
				dev.off()
		}
		
		DBiNorm <- function(ias, mean = 0, sd = 1, positiveRatio = 0.5) {
			positiveRatio*dnorm(ias, mean = mean, sd = sd) +
				 (1-positiveRatio)*dnorm(ias, mean = -mean, sd = sd)
		}
		mean <- mean(abs(ias))
		sd <- sqrt(BVar(abs(ias)))
		if (fixSignApproximation[i]){
			df$PBiB <- DBiNorm(ias, mean = mean, sd = sd)
			BiSymAIC[i] <- -2*sum(log(df$PBiB)) + 2*2
			BiSymAICC[i] <- -2*sum(log(df$PBiB)) + 2*AICC(2,n)		
		}else{
			if (is.na(meanStartSymmetric[i]))
				meanStartSymmetric[i] <- mean
			if (is.na(sdStartSymmetric[i]))
				sdStartSymmetric[i] <- sd
			BiSymFunc <- function(mean = 2, sd = 2)
					-sum(log(DBiNorm(ias, mean = mean, sd = sd)))			
			BiSymFit <- mle(BiSymFunc, start = list(mean = meanStartSymmetric[i], 
				sd = sdStartSymmetric[i]))						
			BiSymAIC[i] <- -2*as.numeric(logLik(BiSymFit)) + 2*2
			BiSymAICC[i] <- -2*as.numeric(logLik(BiSymFit)) + 2*AICC(2,n)			
			mean = coef(BiSymFit)["mean"]
			sd = coef(BiSymFit)["sd"]
		}
		mean3[i] <- mean
		sd3[i] <- sd
		switch(i,
			f <- function(x){DBiNorm(x, mean = mean3[1], sd = sd3[1])},
			f <- function(x){DBiNorm(x, mean = mean3[2], sd = sd3[2])},
			f <- function(x){DBiNorm(x, mean = mean3[3], sd = sd3[3])},
			f <- function(x){DBiNorm(x, mean = mean3[4], sd = sd3[4])},
			f <- function(x){DBiNorm(x, mean = mean3[5], sd = sd3[5])},
			f <- function(x){DBiNorm(x, mean = mean3[6], sd = sd3[6])},
			f <- function(x){DBiNorm(x, mean = mean3[7], sd = sd3[7])},
			f <- function(x){DBiNorm(x, mean = mean3[8], sd = sd3[8])},
			f <- function(x){DBiNorm(x, mean = mean3[9], sd = sd3[9])},
			f <- function(x){DBiNorm(x, mean = mean3[10], sd = sd3[10])},
		)
		BimodalSymmetric[[i]] <- list(AIC = BiSymAIC[i], AICc = BiSymAICC[i], 
			mean = mean, sd = sd, f = f)
		names(BimodalSymmetric)[i] = colname
		if (plotGraph){
			fileName = sprintf("%s-%s-Bimodal-Symmetric.pdf", 
				filePrefix, colname)
			if (!plotToScreen)
				pdf(file = fileName)
			PlotHist(ylimMin[i], ylimMax[i])
			lines(iv, DBiNorm(iv, mean = mean, sd = sd))
			if (plotToScreen)			
				dev.print(pdf,file = fileName)
			else
				dev.off()
		}	
		
		mean <- mean(abs(ias))
		sd <- sqrt(BVar(abs(ias)))
		positiveRatio <- sum(ias > 0)/n
		if (fixSignApproximation[i]){
			df$PBiAsymB <- DBiNorm(ias, mean = mean, sd = sd, positiveRatio = positiveRatio)
			BiAsymAIC[i] <- -2*sum(log(df$PBiAsymB)) + 2*3
			BiAsymAICC[i] <- -2*sum(log(df$PBiAsymB)) + 2*AICC(3,n)
		}else{
			if (is.na(meanStartAsymmetric[i]))
				meanStartAsymmetric[i] <- mean
			if (is.na(sdStartAsymmetric[i]))
				sdStartAsymmetric[i] <- sd
			if (is.na(positiveRatioStartAsymmetric[i]))
				positiveRatioStartAsymmetric[i] <- positiveRatio
			BiAsymFunc <- function(mean = 2, sd =1 , positiveRatio = 0.5)
					-sum(log(DBiNorm(ias, mean = mean, sd = sd, positiveRatio = positiveRatio)))
			#BiAsymFit <- mle(BiAsymFunc,lower=c(-100,0,0), upper=c(100,50,1)) 
				#bounds can only be used with method L-BFGS-B or Brent
			BiAsymFit <- mle(BiAsymFunc, start = list(mean = meanStartAsymmetric[i], 
				sd = sdStartAsymmetric[i], positiveRatio = positiveRatioStartAsymmetric[i]))
			BiAsymAIC[i] <- -2*as.numeric(logLik(BiAsymFit)) + 2*3
			BiAsymAICC[i] <- -2*as.numeric(logLik(BiAsymFit)) + 2*AICC(3,n)
			mean <- coef(BiAsymFit)["mean"]
			sd <- coef(BiAsymFit)["sd"]
			positiveRatio <- coef(BiAsymFit)["positiveRatio"]
		}
		mean4[i] <- mean
		sd4[i] <- sd
		positiveRatio4[i] <- positiveRatio
		switch(i,
			f <- function(x){DBiNorm(x, mean = mean4[1], sd = sd4[1], positiveRatio4[1])},
			f <- function(x){DBiNorm(x, mean = mean4[2], sd = sd4[2], positiveRatio4[2])},
			f <- function(x){DBiNorm(x, mean = mean4[3], sd = sd4[3], positiveRatio4[3])},
			f <- function(x){DBiNorm(x, mean = mean4[4], sd = sd4[4], positiveRatio4[4])},
			f <- function(x){DBiNorm(x, mean = mean4[5], sd = sd4[5], positiveRatio4[5])},
			f <- function(x){DBiNorm(x, mean = mean4[6], sd = sd4[6], positiveRatio4[6])},
			f <- function(x){DBiNorm(x, mean = mean4[7], sd = sd4[7], positiveRatio4[7])},
			f <- function(x){DBiNorm(x, mean = mean4[8], sd = sd4[8], positiveRatio4[8])},
			f <- function(x){DBiNorm(x, mean = mean4[9], sd = sd4[9], positiveRatio4[9])},
			f <- function(x){DBiNorm(x, mean = mean4[10], sd = sd4[10], positiveRatio4[10])},
		)
		BimodalAsymmetric[[i]] <- list(AIC = BiAsymAIC[i], AICc = BiAsymAICC[i], 
			mean = mean, sd = sd, positiveRatio = positiveRatio, f = f)
		names(BimodalAsymmetric)[i] = colname
		if (plotGraph){
			fileName = sprintf("%s-%s-Bimodal-Asymmetric.pdf", 
				filePrefix, colname)
			if (!plotToScreen)
				pdf(file = fileName)
			PlotHist(ylimMin[i], ylimMax[i])
			lines(iv, DBiNorm(iv, mean = mean, sd = sd, positiveRatio = positiveRatio))
			if (plotToScreen)			
				dev.print(pdf,file = fileName)
			else
				dev.off()
		}
	}
	
	PrintNameValue <- function(s, x){
		cat(sprintf("   %-20s", s), sprintf("%9.3f", x), "\n")
	}
	cat(sprintf("%-20s   ", "AIC:"))	
	for (col in cols)
		cat(sprintf("%10s", colnames(df)[col]))
	cat("\n")		
	PrintNameValue("Unimodal Symmetric", UniSymAIC)
	PrintNameValue("Unimodal Asymmetric", UniAsymAIC)
	PrintNameValue("Bimodal Symmetric", BiSymAIC)
	PrintNameValue("Bimodal Asymmetric", BiAsymAIC)	
	cat(sprintf("%-20s   ", "AICc:"))	
	#for (col in cols)
	#	cat(sprintf("%10s", colnames(df)[col]))
	cat("\n")		
	PrintNameValue("Unimodal Symmetric", UniSymAICC)
	PrintNameValue("Unimodal Asymmetric", UniAsymAICC)
	PrintNameValue("Bimodal Symmetric", BiSymAICC)
	PrintNameValue("Bimodal Asymmetric", BiAsymAICC)
	
	return(list(UnimodalSymmetric = UnimodalSymmetric, 
		UnimodalAsymmetric = UnimodalAsymmetric,
		BimodalSymmetric = BimodalSymmetric, 
		BimodalAsymmetric = BimodalAsymmetric))
}
