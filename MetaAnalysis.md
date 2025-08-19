---
name: MetaAnalysis
topic: Meta-Analysis
maintainer: Michael Dewey, Wolfgang Viechtbauer
email: lists@dewey.myzen.co.uk
version: 2025-05-27
source: https://github.com/cran-task-views/MetaAnalysis/
---

This task view covers packages which include facilities for meta-analysis of summary statistics from primary studies. The task view does not consider the meta-analysis of individual participant data (IPD) which can be handled by any of the standard linear modelling functions but it does include some packages which offer special facilities for IPD.

The standard meta-analysis model is a form of weighted least squares and so any of the wide range of R packages providing weighted least squares would in principle be able to fit the model. The advantage of using a specialised package is that (a) it takes care of the small tweaks necessary (b) it provides a range of ancillary functions for displaying and investigating the model. Where the model is referred to below, it is this model which is meant.

Where summary statistics are not available, a meta-analysis of significance levels is possible. This is not completely unconnected with the problem of adjustment for multiple comparisons but the packages below which offer this, chiefly in the context of genetic data, also offer additional functionality.

### Univariate meta-analysis

#### Preparing for meta-analysis

- The primary studies often use a range of statistics to present their results. Convenience functions to convert these onto a common metric are presented by: `r pkg("compute.es")` which converts from various statistics to d, g, r, z, and the log odds ratio, `r pkg("MAd")` which converts to mean differences, and `r pkg("metafor", priority = "core")` which converts to an extensive set of effect size measures for comparative studies (such as binary data, person years, mean differences and ratios, and so on), for studies of association (a wide range of correlation types), and for non-comparative studies (proportions, incidence rates, and mean change). It also provides for a measure used in psychometrics (Cronbach's alpha). `r pkg("esc")` provides a range of effect size calculations with partial overlap with `r pkg("metafor")` but with some extras, noticeably for converting test statistics, and also includes a convenience function for collating its output for input to another package like `r pkg("metafor")` or producing a CSV file. `r pkg("estimraw")` estimates the cell frequencies from odds ratios, risk ratios, or risk differences. `r pkg("effsize")` contains functions to compute mean difference effect sizes (Cohen's d and Hedges' g) and measures of dominance (Cliff's delta) and stochastic superiority (Vargha-Delaney A). `r pkg("effectsize")` provides a large number of different effect sizes and converts between them. `r pkg("psychmeta")` provides extensive facilities for converting effect sizes and correcting for various restrictions and measurement artifacts. `r pkg("metansue")` provides some methods for converting to effect sizes, while `r pkg("es.dif")` computes Cohen's d, Hedges' g, biased/unbiased c (an effect size between a mean and a constant), and e (an effect size between means without assuming variance equality) from raw data. `r pkg("MOTE")` provides a variety of conversions based on Cohen's d. `r pkg("estmeansd")` converts between quantiles and means and standard deviations. `r pkg("metaBLUE")` estimates means and standard deviations from various order statistics. `r pkg("SingleCaseES")` provides basic effect sizes for single-case designs including parametric and non-overlap measures. `r pkg("smd")` computes standardised mean differences. `r pkg("metaHelper")` calculates values commonly used in meta-analysis with a focus on standardized mean differences and related statistics. `r pkg("metaConvert")` estimates 11 effect size measures.
- `r pkg("meta", priority = "core")` provides functions to read and work with files output by RevMan 4 and 5.
- `r pkg("metagear")` provides many tools for the systematic review process including screening articles, downloading the articles, generating a PRISMA diagram, and some tools for effect sizes. `r pkg("revtools")` provides tools for downloading from bibliographic databases and uses machine learning methods to process them. `r pkg("citationchaser")` assists in the process of chasing citations. `r pkg("MetaNLP")` pre--processes titles and abstracts for use in screening before a meta-analysis.
- `r pkg("esci")` also provides a large range of effect sizes.
- `r pkg("metacor")` calculates effect sizes in pre-post studies and imputes missing variances and correlations.
- `r pkg("metavcov")` computes the variance-covariance matrix for multivariate meta-analysis when correlations between outcomes can be provided but not between treatment effects.
- `r pkg("clubSandwich")` and `r pkg("metafor")` provide functions to impute variance-covariance matrix for multivariate meta-analysis.
- `r pkg("metafuse")` uses a fused lasso to merge covariate estimates across a number of independent datasets.
- `r pkg("metapower")` provides power analysis for meta-analysis and meta-regression. `r pkg("POMADE")` does the same for the overall average effect size in a meta-analysis of dependent effect sizes.
- `r pkg("PRISMA2020")` produces an interactive flow diagram that conforms to PRISMA 2020 version, `r pkg("PRISMAstatement")` also generates flowcharts conforming to the PRISMA statement.
- `r pkg("reappraised")` provides tools for checking the integrity of groups of trials.
- `r pkg("RTSA")` provides Trial Sequential Analysis.
- `r pkg("AIscreenR")` uses a large language model to screen titles and abstracts.
- Several packages provide assistance in digitising data from published figures: `r pkg("metaDigitise")` and `r pkg("juicr")` provide graphical interfaces and accept various input formats. `r pkg("digitize")` has a more limited range of facilities.

#### Fitting the model

- Four packages provide the inverse variance weighted, Mantel-Haenszel, and Peto methods: `r pkg("meta")`, `r pkg("metafor")`, `r pkg("rmeta")`, and `r pkg("epiR")`.
- For binary and time-to-event data, `r pkg("metafor")` also provides binomial-normal and the Poisson-normal models.
- Packages which work with specific effect sizes may be more congenial to workers in some areas of science. `r pkg("MAd")` provides functions for the meta-analysis of standardized mean differences and provides a range of graphics.
- `r pkg("psychmeta")` implements the Hunter-Schmidt method (also known as psychometric meta-analysis) including corrections for reliability and other artifacts.
- Bayesian approaches are contained in various packages. `r pkg("bspmma")` provides two different models: a non-parametric and a semi-parametric. Graphical display of the results is provided. `r pkg("bayesmeta")` includes shrinkage estimates, meta-regression, posterior predictive p-values, and forest plots via either `r pkg("metafor")` or `r pkg("forestplot")`. Diagnostic graphical output is available. `r pkg("metaBMA")` provides a Bayesian approach using model averaging; a variety of priors are provided and it is possible for the user to define new ones. `r pkg("MetaStan")` includes binomial-normal hierarchical models and can use weakly informative priors for the heterogeneity and treatment effect parameters. `r pkg("baggr")` provides facilities using Stan for hierarchical Bayesian models; graphical facilities are provided. `r pkg("brms")` can also fit Bayesian meta-analytic models using Stan as the backend. `r pkg("BayesCombo")` provides facilities using a Bayesian approach and has graphical facilities. `r pkg("RBesT")` uses Bayesian synthesis to generate priors from various sources. `r pkg("metamisc")` provides a method with priors suggested by Higgins. `r pkg("RoBMA")` provides a framework for estimating ensembles of meta-analytic models and using Bayesian model averaging to combine them it also provides meta-regression. `r pkg("ra4bayesmeta")` provides principled reference analysis within the Bayesian normal-normal model. `r pkg("metabup")` provides a Bayesian approach using basic uncertainty pooling.
- Some packages concentrate on providing a specialised version of the core meta-analysis function without providing a full range of ancillary functions. These are: `r pkg("metaLik")` which uses a more sophisticated approach to the likelihood, `r pkg("metatest")` which provides another improved method of obtaining confidence intervals, and `r pkg("CoTiMA")` which performs meta-analyses of correlation matrices of repeatedly measured variables for studies with different time lags using a SEM framework with OpenMx as the engine.
- `r pkg("metaplus")` fits random effects models relaxing the usual assumption that the random effects have a normal distribution by providing t or a mixture of normals.
- `r pkg("ratesci")` fits random effects models to binary data using a variety of methods for confidence intervals.
- `r pkg("RandMeta")` estimates exact confidence intervals in random effects models using an efficient algorithm.
- `r pkg("rma.exact")` estimates exact confidence intervals in random effects normal-normal models and also provides plots of them.
- `r pkg("pimeta")` provides a range of methods for prediction interval estimation from random effects models and has graphical facilities.
- `r pkg("metamedian")` implements several methods to meta-analyze one-group or two-group studies that report the median as the outcome. These methods estimate the pooled median in the one-group context and the pooled raw difference of medians across groups in the two-group context and also analyses median survival time studies. `r pkg("meta")` and `r pkg("metafor")` also provide methods for medians.
- `r pkg("MetaUtility")` proposes a metric for estimating the proportion of effects above a cut-off of scientific importance.
- `r pkg("metasens")` provides imputation methods for missing binary data.
- `r pkg("metagam")` provides a framework for meta-analysis of generalised additive models including the case where individual participant data cannot be shared across locations and `r pkg("EvidenceSynthesis")` also combines across sites where individual participant data cannot be shared.
- `r pkg("metawho")` implements a method for combining within study interactions.
- `r pkg("metarep")` provides replicability analyses after a conventional analysis.
- `r pkg("rema")` uses a permutation approach to handle meta-analyses of rare event data.
- `r pkg("meta.shrinkage")` uses shrinkage methods to provide better estimates of individual means in meta-analysis.
- `r pkg("metaumbrella")` provides facilities for umbrella reviews.
- `r pkg("vcmeta")` provides functions for varying-coefficient meta-analysis as an alternative to the usual fixed- or random-effects methods.
- `r pkg("robustmeta")` provides methods for meta-analysis for cases where primary studies may have influential outlying values.
- `r pkg("coefa")` provides a method for conducting a meta-analysis of factor analyses based on co-occurrence matrices.
- `r pkg("CausalMetaR")` provides robust and efficient methods for estimating causal effects in a target population using a multi-source dataset.
- `r pkg("metainc")` assesses inconsistency in meta-analyses by calculating the Decision Inconsistency index (DI) and the Across-Studies Inconsistency (ASI) index. Allows input from a variety of packages.
- `r pkg("mars")` performs univariate and multivariate meta-analysis with estimation of the asymptotic variance-covariance matrix.
- `r pkg("twotrials")` deals with the special case of exactly two trials.
- `r pkg("mmeta")` provides for multivariate analysis of two by two tables.
- The Doi model (IVhet) is covered in `r pkg("meta")` and `r pkg("metafor")`.

#### Graphical methods

An extensive range of graphical procedures is available.

- Forest plots are provided in `r pkg("forplo")`, `r pkg("forestly")` (interactive plots), `r pkg("forestmodel")` (using ggplot2), `r pkg("forestplot")`, `r pkg("forestploter")`, `r pkg("meta")`, `r pkg("metafor")`, `r pkg("metansue")`, `r pkg("psychmeta")`, and `r pkg("rmeta")`. Although the most basic plot can be produced by any of them, they each provide their own choice of enhancements. `r pkg("metaviz")` provides a range of enhancements.
- Funnel plots are provided in `r pkg("meta")`, `r pkg("metafor")`, `r pkg("metansue")`, `r pkg("rmeta")` and `r pkg("weightr")`. In addition to standard funnel plots, a funnel plot for limit meta-analysis is provided in `r pkg("metasens")`, and `r pkg("metaviz")` provides an extensive range of enhanced funnel plots and also facilities for their use in the context of visual inference.
- Radial (Galbraith) plots are provided in `r pkg("meta")` and `r pkg("metafor")`.
- L'Abbe plots are provided in `r pkg("meta")` and `r pkg("metafor")`.
- Baujat plots are provided in `r pkg("meta")` and `r pkg("metafor")`.
- `r pkg("meta")` provides drapery plots.
- `r pkg("MetaAnalyser")` provides an interactive visualisation of the results of a meta-analysis.
- `r pkg("metaviz")` provides rainforest plots, an enhanced version of forest plots. It accepts input from `r pkg("metafor")`.
- `r pkg("DTAplots")` produces various plots for diagnostic studies including forest and SROC plots.
- GOSH plots are provided in `r pkg("metafor")`.
- `r pkg("robvis")` can be used to visualize the results of risk-of-bias (rob) assessments.
- `r pkg("xmeta")` provides galaxy plots, an analogue of funnel plots for multivariate meta-analysis.

#### Investigating heterogeneity

- Confidence intervals for the heterogeneity parameter are provided in `r pkg("metafor")`.
- `r pkg("altmeta")` presents a variety of alternative methods for measuring and testing heterogeneity with a focus on robustness to outlying studies.
- `r pkg("mc.heterogeneity")` implements a Monte Carlo based test for heterogeneity.
- `r pkg("boot.heterogeneity")` provides a bootstrap test for heterogeneity for mean differences, correlations, and odds ratios.
- `r pkg("heterometa")` converts between various summary measures of heterogeneity.

#### Model criticism

- An extensive series of plots of diagnostic statistics for detecting outliers and influential studies is provided in `r pkg("metafor")`.
- `r pkg("metaplus")` provides outlier diagnostics.
- `r pkg("EValue")` provides a sensitivity analysis of the effect of unmeasured confounders.
- `r pkg("boutliers")` provides bootstrap distributions for outlier detection and influence diagnostics.
- `r pkg("metaconfoundr")` provides a number of ways to visualise confounding relationships in meta-analysis.
- `r pkg("RoBMA")` allows comparison of different meta-analytic models.

#### Small study bias / Unobserved studies / Publication bias

The issue of whether small studies give different results from large studies can be addressed by visual examination of the funnel plots mentioned above. In addition:

- `r pkg("meta")` and `r pkg("metafor")` provide both the non-parametric rank correlation test suggested by Begg and Mazumdar and a range of regression tests modelled after the approach of Egger.
- `r pkg("metamisc")` provides funnel plots and tests for asymmetry.
- `r pkg("xmeta")` provides methods for small study effects in multivariate meta-analysis.

A related issue in meta-analysis is the problem of unobserved studies and publication bias.

- Rosenthal's fail safe n is provided by `r pkg("MAd")`. `r pkg("metafor")` provides it as well as two more recent methods by Orwin and Rosenberg and a generalization to random-effects models.
- `r pkg("fsn")` computes the fail-safe number with confidence interval.
- Duval's trim and fill method is provided by `r pkg("meta")` and `r pkg("metafor")`.
- `r pkg("metasens")` provides the Copas selection model and also the method of limit meta-analysis (a regression based approach for dealing with small study effects) due to Rücker and colleagues.
- `r pkg("RobustBayesianCopas")` fits a robust version of the Copas selection model.
- `r pkg("selectMeta")` provides various selection models: the parametric model of Iyengar and Greenhouse, the non-parametric model of Dear and Begg, and proposes a new non-parametric method imposing a monotonicity constraint.
- `r pkg("weightr")` provides facilities for using the weight function model of Vevea and Hedges. `r pkg("metafor")` also provides an implementation thereof and a variety of other selection models.
- `r pkg("RoBMA")` includes Bayesian versions of selection models.
- `r pkg("phacking")` models only non-affirmative studies to allow for selection of studies and selection within studies.
- `r pkg("puniform")` provides methods using only the statistically significant studies, methods for the special case of replication studies, and sample size determinations.
- `r pkg("PublicationBias")` performs sensitivity analysis of the number of unpublished studies needed to have a specified influence.
- The `r pkg("metansue")` package allows the inclusion by multiple imputation of studies known only to have a non-significant result.
- `r pkg("publipha")` estimates models accounting for publication bias or p-hacking using a Bayesian framework.
- `r pkg("metafor")` provides the test of excess significance.
- `r pkg("multibiasmeta")` conducts sensitivity analyses for the joint effects of internal and publication biases.
- `r pkg("metabias")` provides common components (classes, methods, documentation) for several other packages to investigate within- and across-study biases in meta-analysis.

#### Other study designs

- `r pkg("SCMA")` provides single case meta-analysis. It is part of a suite of packages dedicated to single-case designs.
- `r pkg("joint.Cox")` provides facilities for the meta-analysis of studies of joint time-to-event and disease progression.
- `r pkg("metamisc")` provides for meta-analysis of prognostic studies.
- `r pkg("metamicrobiomeR")` provides meta-analysis of zero-inflated beta microbiome data fitted via GAMLSS models.
- `r pkg("metaSurvival")` estimates the survival curves from data extracted from primary study survival curves.

#### Meta-analysis of significance values

- Fisher's method and Lancaster's are available in `r pkg("aggregation")`, `r pkg("metap")`, `r pkg("metapro")`, and `r pkg("poolr")`.
- Stouffer's, Tippett's, and Wilkinson's method are available in `r pkg("metap")` and `r pkg("poolr")`.
- Edgington's method, inverse-t, logit, mean of p, and mean of z are all available in `r pkg("metap")`.

In all cases `r pkg("poolr")` considers correlated p-values in addition to independent. The others above do not.

- `r pkg("TFisher")` provides Fisher's method using both hard and soft thresholding for the p-values. There is a wrapper in `r pkg("metap")` for the hard threshold case.
- `r pkg("harmonicmeanp")` uses a method based on the harmonic mean of p-values which is robust to correlation between the p-values.
- `r pkg("amanida")` provides meta-analysis of metabolite data using p-values and fold change.
- `r pkg("metap")` provides simple graphics including albatross plots.

Some methods are also provided in some of the genetics packages mentioned below.

### Multivariate meta-analysis

Standard methods outlined above assume that the effect sizes are independent. This assumption may be violated in a number of ways: within each primary study multiple treatments may be compared to the same control, each primary study may report multiple endpoints or multiple assessments, or primary studies may be clustered for instance because they come from the same country or the same research team.

- `r pkg("mvmeta")` assumes the within study covariances are known and provides a variety of options for fitting random effects. `r pkg("metafor")` provides fixed effects and likelihood based random effects model fitting procedures. Both these packages include meta-regression, `r pkg("metafor")` also provides for clustered and hierarchical models.
- `r pkg("mixmeta")` provides an integrated interface to standard meta-analysis and extensions like multivariate and dose-response models.
- `r pkg("mvtmeta")` provides multivariate meta-analysis using the method of moments for random effects although not meta-regression.
- `r pkg("clubSandwich")`, `r pkg("robumeta")`, and `r pkg("metafor")` provide cluster-robust variance estimates for clustered and hierarchical estimates.
- `r pkg("wildmeta")` conducts single coefficient tests and multiple-contrast hypothesis tests of meta-regression models using cluster wild bootstrapping.
- `r pkg("metaSEM")` provides multivariate (and univariate) meta-analysis and meta-regression by embedding it in the structural equation framework and using OpenMx for the structural equation modelling. It can provide a three-level meta-analysis taking account of clustering and allowing for level 2 and level 3 heterogeneity. It also provides via a two-stage approach meta-analysis of correlation or covariance matrices.
- `r pkg("dosresmeta")` concentrates on the situation where individual studies provide information about a dose-response relationship. `r pkg("MBNMAdose")` provides a Bayesian analysis using network meta-analysis of dose-response studies.
- `r pkg("CIAAWconsensus")` has a function for multivariate meta-analysis in the context of atomic weights and estimating isotope ratios.
- `r pkg("BayesMultMeta")` provides Bayesian inference for the parameters of a multivariate random-effects model with application to multivariate meta-analysis.
- `r pkg("remaCor")` provides multivariate meta-analysis when the correlation between effect sizes is known.
- `r pkg("xmeta")` provides a variety of methods for multivariate meta-analysis.
- `r pkg("crwbmetareg")` uses the WLS estimator of Stanley and Doucouliagos both with and without covariates followed by computing the relevant p-values using the cluster robust wild bootstrap methodology.

### Meta-analysis of studies of diagnostic tests

A special case of multivariate meta-analysis is the case of summarising studies of diagnostic tests. This gives rise to a bivariate, binary meta-analysis with the within-study correlation assumed zero although the between-study correlation is estimated. This is an active area of research and a variety of methods are available including what is referred to here as Reitsma's method and the hierarchical summary receiver operating characteristic (HSROC) method. In many situations these are equivalent.

- `r pkg("mada")` provides various descriptive statistics and univariate methods (diagnostic odds ratio and Lehman model) as well as the bivariate method due to Reitsma. Meta-regression is provided.
- `r pkg("bamdit")` provides Bayesian meta-analysis with a bivariate random effects model (using JAGS to implement the MCMC method).
- `r pkg("meta4diag")` provides Bayesian inference analysis for bivariate meta-analysis of diagnostic test studies and an extensive range of graphical methods.
- `r pkg("diagmeta")` considers the case where the primary studies use an analysis using multiple cut-offs.
- `r pkg("NMADiagT")` provides network meta-analysis of diagnostic tests in a Bayesian framework using Stan as the engine.
- `r pkg("DTAplots")` produces various plots for diagnostic studies including forest and SROC plots. The packages above also provide various graphical methods.
- `r pkg("MVPBT")` provides tests for small study effects in diagnostic test meta-analysis.
- `r pkg("dmetatools")` provides confidence intervals for the AUC of the summary ROC curve and related methods.
- `r pkg("CopulaREMADA")` provides the bivariate copula mixed model for meta-analysis of diagnostic test accuracy studies.

### Meta-regression

Where suitable moderator variables are available they may be included using meta-regression. All these packages are mentioned above, this just draws that information together.

- `r pkg("metafor")` provides meta-regression (multiple moderators are catered for). Various packages rely on `r pkg("metafor")` to provide meta-regression (`r pkg("meta")`, and `r pkg("MAd")`) and all of these provide bubble plots. `r pkg("psychmeta")` also uses `r pkg("metafor")`.
- `r pkg("metaLik")`, `r pkg("metansue")`, `r pkg("metaSEM")`, and `r pkg("metatest")` also provide meta-regression.
- `r pkg("mvmeta")` provides meta-regression for multivariate meta-analysis as do `r pkg("metafor")` and `r pkg("metaSEM")`.
- `r pkg("mada")` provides for the meta-regression of diagnostic test studies.
- `r pkg("GENMETA")` uses generalised meta-analysis to handle the situation where the studies do not all use the same regressors.
- `r pkg("jarbes")` uses a Bayesian approach of hierarchical meta-regression.
- `r pkg("metacart")` uses classification and regression trees to identify interactions between moderators.
- `r pkg("metaforest")` investigates heterogeneity using random forests. Note that it has nothing to do with forest plots.
- `r pkg("pema")` provides a penalised approach to meta-regression useful for situations with a large number of moderators relative to observations.
- `r pkg("bayesmeta")` includes meta-regression in a Bayesian framework.
- `r pkg("CAMAN")` offers the possibility of using finite semiparametric mixtures as an alternative to the random effects model where there is heterogeneity. Covariates can be included to provide meta-regression.
- `r pkg("crwbmetareg")` fits meta-regression models using weighted least squares and then uses cluster robust wild bootstrap methodology for inferences.

### Individual participant data (IPD)

Where all studies provide individual participant data, software for the analysis of multi-centre trials or multi-centre cohort studies should prove adequate but is outside the scope of this task view (see the `r view("MixedModels")` task view for relevant packages). Other packages which provide facilities related to IPD are:

- `r pkg("ecoreg")` which is designed for ecological studies enables estimation of an individual level logistic regression from aggregate data or individual data.
- `r pkg("multinma")` provides network meta-analysis and network meta-regression models for aggregate data, individual patient data, and mixtures thereof.
- `r pkg("MetaIntegration")` combines IPD data with external models.
- `r pkg("bipd")` uses a Bayesian approach for IPD. It includes facilities for multiple imputation using mice.

### Network meta-analysis

Also known as multiple treatment comparison, this is a very active area of research and development. Note that some of the packages mentioned above under multivariate meta-analysis can also be used for network meta-analysis with appropriate setup.

- `r pkg("netmeta")` works in a frequentist framework. It provides an extensive range of graphical and other displays including network graphs and a heatmap for displaying inconsistency and heterogeneity. A frequentist analogue of SUCRA is also available.
- A Bayesian framework is provided by `r pkg("pcnetmeta")`, which uses JAGS. It provides a number of data-sets. `r pkg("nmaINLA")` uses integrated nested Laplace approximations as an alternative to MCMC. It provides a number of data-sets. `r pkg("NMADiagT")` provides network meta-analysis of diagnostic tests in a Bayesian framework using Stan as the engine; graphical output is provided. `r pkg("gemtc")` acts as a front-end to BUGS or JAGS, `r pkg("bnma")` provides arm-based methods using JAGS as the engine, and `r pkg("metapack")` provides methods using built-in MCMC code.
- `r pkg("multinma")` provides network meta-analysis and network meta-regression models for aggregate data, individual patient data, and mixtures thereof.
- `r pkg("netdose")` provides NMA of does-response studies in a frequentist way.
- `r pkg("nmathresh")` provides decision-invariant bias adjustment thresholds and intervals the smallest changes to the data that would result in a change of decision. `r pkg("NMAoutlier")` detects outliers in NMA using a forward search.
- `r pkg("pcnetmeta")` provides network graphs. `r pkg("nmaplateplot")` displays the results from an NMA using a heatplot style and also displays SUCRA.
- `r pkg("nmarank")` evaluates hierarchies of evidence in network meta-analysis.
- `r pkg("rnmamod")` uses a Bayesian approach to perform NMA while addressing (aggregate) missing participant outcome data and provides graphics.
- `r pkg("viscomp")` provides several visualization tools for exploring the behavior of the components in a network meta-analysis of multi-component interventions.
- `r pkg("OssaNMA")` calculates minimum total sample size needed for prespecified power or optimal allocation for each treatment group with a fixed total sample size.
- `r pkg("MBNMAtime")` provides for analysis of multiple time points from studies using a Bayesian framework.
- `r pkg("crossnma")` can be used to conduct a network meta-analysis for individual participant data, aggregate data, and mixtures thereof in a Bayesian framework.
- `r pkg("PINMA")` provides improved methods to construct prediction intervals for network meta-analysis.
- `r pkg("rankinma")` provides treatment ranking in NMA.
- `r pkg("ssifs")` evaluates consistency in NMA.
- `r pkg("NMA")` uses a contrast-based approach and includes standard diagnostic and graphical methods. It uses an improved REML estimation procedure.
- `r pkg("CBnetworkMA")` performs contrast-based NMA in a Bayesian framework.
- `r pkg("closeloop")` calculate distance between single-arm observational studies using covariate information to remove heterogeneity in NMA.

### Genetics

There are a number of packages specialising in genetic data: `r pkg("catmap")` combines case-control and family study data, graphical facilities are provided, `r pkg("CPBayes")` uses a Bayesian approach to study cross-phenotype genetic associations, `r pkg("corrmeta")` performs correlated meta-analysis across multiple scans, `r pkg("gap")` combines p-values, `r pkg("getmstatistic")` quantifies systematic heterogeneity, `r pkg("getspres")` uses standardised predictive random effects to explore heterogeneity in genetic association meta-analyses, `r pkg("GMCM")` uses a Gaussian mixture copula model for high-throughput experiments, `r pkg("GSEMA")` performs the different steps of gene set enrichment meta-analysis, `r pkg("MendelianRandomization")` provides several methods for performing Mendelian randomisation analyses with summarised data, `r pkg("metaGE")` provides functions for a meta-analysis of genome-wide association studies for studying genotype x environment interactions. `r pkg("MetaHD")` performs multivariate meta-analysis for high-dimensional metabolomics data, `r pkg("MetaIntegrator")` provides meta-analysis of gene expression data, `r pkg("metaMA")` provides meta-analysis of p-values or moderated effect sizes to find differentially expressed genes, `r pkg("metaRNASeq")` does meta-analysis from multiple RNA sequencing experiments, `r pkg("MetaSKAT")` provides for meta-analysis of the SKAT, `r pkg("MetaSubtract")` uses leave-one-out methods to validate meta-GWAS results, `r pkg("ofGEM")` provides a method for identifying gene-environment interactions using meta-filtering, `r pkg("RobustRankAggreg")` provides methods for aggregating lists of genes, and `r pkg("SPAtest")` combines association results.

### Data-sets

- `r pkg("metadat")` provides a large number of data-sets used in meta-analysis.
- `r pkg("psymetadata")` provides more data-sets from meta-analyses in psychology research.
- `r pkg("nmadb")` provides access to a database of network meta-analyses.
- `r pkg("tracenma")` provides a database of network systematic reviews for use in checking transitivity.
- `r pkg("KenSyn")` provides data-sets to accompany a French language book on meta-analysis in the agricultural sciences.
- `r pkg("metabolic")` provides data and code for a published meta-analysis.

### Interfaces

- Plug-ins for Rcmdr are provided by: `r pkg("RcmdrPlugin.EZR")` which uses `r pkg("meta")` and `r pkg("metatest")`, `r pkg("RcmdrPlugin.MA")` which uses `r pkg("MAd")` and `r pkg("metafor")`, and `r pkg("RcmdrPlugin.RMTCJags")` for network meta-analysis using BUGS code.
- Both [JASP](https://jasp-stats.org) and [jamovi](https://www.jamovi.org) provide meta-analysis with R as the backend (the latter via the [MAJOR](https://github.com/kylehamilton/MAJOR) module).
- `r pkg("miniMeta")` provides a shiny interface to `r pkg("meta")`.

### Links

[R mailing list for meta-analysis](https://stat.ethz.ch/mailman/listinfo/r-sig-meta-analysis/)
