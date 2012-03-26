#-------------------------------------- Function to calculate cross spectrum from correlation function
corr2spec <- function( corr ){
	nspec <- length(corr)/2
	spec <- fft(corr)[1:nspec]
	return(spec)
}
#-------------------------------------- Function to calculate correlation function from cross spectrum
spec2corr <- function(spec){
	nspec <- length(spec)
	tmpspec <- c(spec, 0, Conj(spec[nspec:2]))
	corr <- Re(fft(tmpspec, inverse=TRUE) / nspec)
	return(c(corr[(nspec+1):(2*nspec)], corr[1:nspec]))
}

#-------------------------------------- Delay Search
delay_search <- function( spec ){
	nspec <- length( spec )
	#-------- Search for delay
	delay <- as.numeric(which.max(Mod(spec2corr(spec))) - nspec )	# Coarse Delay
	delay_step <- 1													# Initial step (in unit of sampling period)
	while( delay_step > 0.01){
		trial_delay <- seq(-delay_step, delay_step, by=delay_step); trial_amp <- numeric(3)	# Trials with a 3-lag window
		for(index in 1:3){ trial_amp[index] <- Mod(sum(delay_cal(spec, delay + trial_delay[index]))) } # Amplitudes for the 3-lag window
		temp_df <- data.frame( x=trial_delay, y=trial_amp )    # Store them into a data frame
		fit <- lm( data = temp_df, formula = y ~ x + I(x^2) )  # Fit for a quadratic function
		
		delay <- delay - 0.5*fit[[1]][2]/fit[[1]][3]	# Correction for the best delay value
		delay_step <- delay_step / 4					# Narrower steps
	}
	return(list(delay = delay, phase = Arg(sum(delay_cal(spec, delay)))))
}

#-------------------------------------- Delay calibration
delay_cal <- function( spec, delay ){
	# spec : input spectrum (complex)
	# delay : delay[1] = initial phase, delay[2] = delay
	# delay_cal() returns delay-calibrated spectrum
	#
	nspec <- length( spec )
	twiddle <- complex(modulus=rep(1, nspec), argument = delay* pi* seq((-nspec/2), (nspec/2 - 1), by=1) / nspec)
	return( spec * twiddle )
}
