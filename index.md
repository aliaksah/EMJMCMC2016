## Bayesian model configuration, selection and averaging in complex regression contexts


In this R package problems of Bayesian model selection and model averaging are addressed in various complex regression contexts. The approaches developed within the package are based on the idea of marginalizing out parameters from the likelihood. This allows to work on the marginal space of models, which simplifies the search algorithms significantly. For the generalized linear mixed models an efficient mode jumping Monte Carlo Markov chain (MJMCMC) algorithm is implemented. The approach performs very well on simulated and real data. Further, the algorithm is extended to work with logic regressions, where one has a feature space consisting of various complicated logical expressions, which makes enumeration of all features computationally and memory infeasible in most of the cases. The genetically modified MJMCMC (GMJMCMC) algorithm is simplemented suggested to tackle this issue. The algorithm combines the idea of keeping and updating the populations of highly predictive logical expressions combined with MJMCMC for the efficient exploration of the model space. Several simulation and real data studies show that logical expressions of high orders can be recovered with large power and low false discovery rate. Moreover, the GMJMCMC approach is adapted to make inference within the class of deep Bayesian regression models (which is a suggested in the thesis extension of various machine and statistical learning models like artificial neural networks, classification and regression trees, logic regressions and linear models). The reversible GMJMCMC, named RGMJMCMC, is also suggested. It makes transitions between the populations of variables in a way that satisfies the detailed balance equation. Based on several examples, it is shown that the DBRM approach can be efficient for both inference and prediction in various applications. In particular, two ground physical laws (planetary mass law and third Kepler’s law) can be recovered from the data with large power and low false discovery rate. Three classification examples are also studied, where the comparison to other popular machine and statistical learning approaches is performed.

* Full text of the paper introducing MJMCMC for Bayesian variable selection: [arXiv](http://arxiv.org/abs/1604.06398)
* Full text of the paper introducing GMJMCMC for inference on Bayesian logic regressions: [arXiv](https://arxiv.org/abs/1705.07616)
* Full text of the paper introducing DBRM and GMJMCMC, RGMJMCMC algorithms for DBRM: [arXiv](https://arxiv.org/abs/1806.02160)
* Presentations of the talks are available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/tree/master/presentations)
* Latest issues of the package are available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/)
* Some applications of the package are available on [GitHub](https://github.com/aliaksah/EMJMCMC2016/tree/master/examples/)  

***

![Concept](https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/master/illustrations/opt_symmetric.png)

## Developed by [Aliaksandr Hubin](https://scholar.google.com/citations?user=Lx-G8ckAAAAJ&hl=en/), [Geir Storvik](https://scholar.google.no/citations?user=0xDw_sQAAAAJ&hl=en) and [Florian Frommlet](https://scholar.google.com/citations?user=Nmh2LqgAAAAJ&hl=en)

***
**N/B To install the latest version run `install_github("aliaksah/EMJMCMC2016")` or choose the version of interest on https://github.com/aliaksah/EMJMCMC2016/  and install it directly `install.packages("https://github.com/aliaksah/EMJMCMC2016/raw/master/EMJMCMC_1.4_bin.tar.gz", repos = NULL, type="source")`.**

* Notice that some dependencies might be required. The full installation would need some additional [packages](https://github.com/aliaksah/EMJMCMC2016/blob/master/examples/install/install.R)

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
