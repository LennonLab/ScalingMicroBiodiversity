##Do microbes and macrobes differ in how *N* relates to aspects of local diversity (dominance, evenness, rarity)?

## MODEL: Dominance ~ log(N) + microORmacro                            
	
	R-squared = 0.926
	Adj. R-squared = 0.926
	F = 8505;  p = 0.00
	Log-Likelihood = -1347.1
	AIC = 2700 ; BIC = 2716
	
	No. Observations = 1365   
	Df Residuals = 1362   
	Df Model = 2                                         

                 coef   std err      t       P>|t|      95.0% CI
	-------------------------------------------------------------------
	Intercept  -1.5149   0.065    -23.179    0.000   -1.643  -1.387
	T.micro    -0.2889   0.047     -6.121    0.000   -0.381  -0.196
	N           1.0065   0.010    101.813    0.000    0.987   1.026

	Durbin-Watson = 1.315 
	Omnibus = 102.741; p = 0.000
	Cond. No.  32.8
	*Reject normality*
	
	Skew = -0.645   
	Kurtosis = 3.843
	Jarque-Bera = 134.941; p = 4.99e-30
	*Reject normality*
	
====

## MODEL: Evenness ~ log(N) + microORmacro

	R-squared = 0.638
	Adj. R-squared = 0.638
	F = 1201;  p = 1.98e-301
	Log-Likelihood = -503.80
	AIC = 1014 ; BIC = 1029
	
	No. Observations = 1365   
	Df Residuals = 1362   
	Df Model = 2                                         

                 coef   std err      t       P>|t|      95.0% CI
	-------------------------------------------------------------------
	Intercept  0.3642    0.035     10.462    0.000     0.296   0.432
	T.micro    0.2995    0.025     12.046    0.000     0.251   0.348
	N         -0.2292    0.005    -44.204    0.000    -0.239  -0.219

	Durbin-Watson = 0.788 
	Omnibus = 12.901; p = 0.002
	Cond. No.  32.2
	*Reject normality*

	
	Skew = -0.234   
	Kurtosis = 3.093
	Jarque-Bera = 12.984; p = 0.00152
	*Reject normality*
	
====

## MODEL: Rarity ~ log(N) + microORmacro

	R-squared = 0.455
	Adj. R-squared = 0.454
	F = 568.7;  p = 2.81e-180
	Log-Likelihood = -3871.4
	AIC = 7749 ; BIC = 7764
	
	No. Observations = 1365   
	Df Residuals = 1362   
	Df Model = 2                                         

                 coef   std err      t       P>|t|      95.0% CI
	-------------------------------------------------------------------
	Intercept  -3.8229  0.406     -9.412      0.000    -4.620  -3.026
	T.micro    3.1206   0.292     10.682      0.000     2.548   3.694
	N         1.0810    0.061     17.836      0.000     0.962   1.200

	Durbin-Watson = 1.237 
	Omnibus = 760.585; p = 0.000
	Cond. No.  31.8
	*Reject normality*
	
	Skew = 2.223
	Kurtosis = 17.220
	Jarque-Bera = 12625.152; p = 0.000
	*Reject normality*