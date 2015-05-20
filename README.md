Estimating Heterogenous Treatments and Effects using Super Learning

The 'HetEffects' package uses the method outlined in Grimmer, Messing and Westwood (2014) to estimate heterogenous treatments and heterogenous treatment effects.  An ensemble of constituent models are fit based on user-supplied formulas, weighted according to cross-validated model performance, then used to predict the expected marginal conditional averages.  Then the treatment and control are differenced, conditional on key user-specified covariates to produce Marginal Conditional Average Treatment Effects (MCATEs).

Additional details can be found at http://stanford.edu/~jgrimmer/het.pdf. 
