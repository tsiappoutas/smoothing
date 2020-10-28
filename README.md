# smoothing
Jump Smoothing Algorithm                        
A geospatial territory rating algorithm for smoothing relativities (actuarial factors) at a zip code level is demonstrated with code in R. The inputs are the relativities and exposure at a zip code level; the output is a smoothed relativity at a zip code level. The method is algorithmic, not model-based, which makes it easily implementable, interpretable, and explainable to the Departments of Insurance. It is particularly suitable small to medium sized insurance carriers who have thin, noncredible data at a zip code level. This is because the smoothing is not based on Pure Premium (Loss Cost) or Loss Ratio predictions, which is unreliable with very thin data especially at small coverages, but rather, each zip code relativity is weighted by the exposure-weightedaverage relativity of all its neighboring zip codes. R code and data available on GitHub.
