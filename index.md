Generalized linear mixed models ([GLMM](https://en.wikipedia.org/wiki/Generalized_linear_mixed_model)) are addressed for inference and
prediction in a wide range of different applications providing a powerful
scientific tool for the researchers and analysts coming from different
fields. In most of these fields more and more sources of data are becoming
available introducing a variety of hypothetical explanatory variables for
these models to be considered. Selection of an optimal combination of these
variables is thus becoming crucial in a Bayesian setting. The posterior
distribution of the models can be viewed as a relevant measure for the model
evidence, based on the observed data. The number of models to select from is
exponential in the number of candidate variables, moreover the search space
in this context is often extremely non-concave and has numerous local
extrema or statistically speaking modes. Hence efficient search algorithms
have to be adopted for evaluating the posterior distribution within a
reasonable amount of time. In this paper we introduce and implement
efficient mode jumping MCMC algorithms for calculating posterior
probabilities of the models for generalized linear models with a random
effect. Marginal likelihoods of models, given the specific choice of priors
and any choice of covariates, can be efficiently calculated using the
integrated nested Laplace approximations approach ([INLA](http://www.r-inla.org/)) for the class of
models addressed, however for some particular cases exact results are also
available. We further apply the suggested algorithm to some simulated [data1](https://github.com/aliaksah/EMJMCMC2016/tree/master/examples/Simulated%20Data%20%28Example%201%29) and [data2](https://github.com/aliaksah/EMJMCMC2016/tree/master/examples/Simulated%20Logistic%20Data%20With%20Multiple%20Modes%20%28Example%203%29),
the famous [U.S. crime data](https://github.com/aliaksah/EMJMCMC2016/tree/master/examples/US%20Data), [protein activity data](https://github.com/aliaksah/EMJMCMC2016/tree/master/examples/Protein%20Activity%20Data), [real epigenetic data](https://github.com/aliaksah/EMJMCMC2016/tree/master/examples/Epigenetic%20Data), and [NEO-nonNEO asteroids data](https://github.com/aliaksah/EMJMCMC2016/tree/master/examples/asteroid%20data) and compare its performance to some of the existing approaches like [BAS](https://cran.r-project.org/web/packages/BAS/index.html), RS
or MC3. 

* Full text of the paper is available on [arXiv](http://arxiv.org/abs/1604.06398)

* Short version of the paper is available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/blob/master/paper/paper_short.pdf)

* Proofs and pseudo-codes are available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/blob/master/paper/appendix.pdf)

* Presentations of the talks are available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/tree/master/presentations)

* Latest issues of the package are available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/issues)

* EMJMCMC installation procedure is available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/blob/master/examples/install/install.R)

* Examples of the use of EMJMCMC are available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/tree/master/examples/)  

***

![Concept](https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/illustrations/opt_symmetric.png)

## _** Developed by [Aliaksandr Hubin](http://www.mn.uio.no/math/english/people/aca/aliaksah/) and [Geir Storvik](http://www.mn.uio.no/math/personer/vit/geirs/)**_

***
**N/B To install the latest version run `install.packages("https://github.com/aliaksah/EMJMCMC2016/files/270429/EMJMCMC_1.2.tar.gz", repos = NULL, type="source")` in R.**

* Notice that some dependencies might be required. The full installation would need some additional [packages](https://github.com/aliaksah/EMJMCMC2016/blob/master/examples/install/install.R)

* Usage example can be found [here](https://github.com/aliaksah/EMJMCMC2016/blob/master/examples/US%20Data/mode_jumping_package_class_crime_bas_data_3237.r)

***

The research was presented during several talks

1. [The first annual conference of NORBIS, Rosendal, Norway, October 2015](http://norbis.no/files/2015/03/Full-program-NORBIS-Annual-Meeting.pdf)
2. Statistics in Genomics Discussion Group Meeting, Oslo, Norway, November 2015
3. [8th Conference on Computational and Methodological Statistics, London, UK, December 2015](http://cmstatistics.org/RegistrationsV2/CFE2015/viewSubmission.php?id=1533&token=044snso7ns3q3041q240qr64s2o38p20)
4. [Colloquium Talk at Wiener Biometrische Sektion, Vienna, Austria, March 2016](http://www.meduniwien.ac.at/wbs/kolloquien.html)
5. [Colloquium Talk at Belarus State University, Minsk, Belarus, March 2016](http://www.fpmi.bsu.by/ImgFpmi/Cache/Page/15303.pdf)
6. [Game of Epigenetics Conference, Dubrovnik, Croatia, April 2016](http://goe.irb.hr/Programme/Variable-selection-in-binomial-regression-with-latent-Gaussian-field-models-for-analysis-of-epigenetic-data)
7. Astrophysics talk (Department of Astrophysics, University of Oslo), Oslo, Norway, May 2016
8. [NordStat 2016 conference in mathematical statistics, Copenhagen, Denmark, June 2016](http://nordstat2016.dk/posterabstracts.php#1) 
9. [11th International Conference
COMPUTER DATA ANALYSIS & MODELING 2016
Theoretical & Applied Stochastics
, Minsk, Belarus, September 2016](http://www.cdam.bsu.by/en/sm.aspx?guid=3033)

***
