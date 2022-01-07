---
name: MetaAnalysis
topic: Meta-Analysis
maintainer: Michael Dewey
email: lists@dewey.myzen.co.uk
version: 2022-01-07
source: https://github.com/cran-task-views/MetaAnalysis/
---

This task view covers packages which include facilities for
meta-analysis of summary statistics from primary studies. The task view
does not consider the meta-analysis of individual participant data (IPD)
which can be handled by any of the standard linear modelling functions
but it does include some packages which offer special facilities for
IPD.

The standard meta-analysis model is a form of weighted least squares and
so any of the wide range of R packages providing weighted least squares
would in principle be able to fit the model. The advantage of using a
specialised package is that (a) it takes care of the small tweaks
necessary (b) it provides a range of ancillary functions for displaying
and investigating the model. Where the model is referred to below it is
this model which is meant.

Where summary statistics are not available a meta-analysis of
significance levels is possible. This is not completely unconnected with
the problem of adjustment for multiple comparisons but the packages
below which offer this, chiefly in the context of genetic data, also
offer additional functionality.

#### Univariate meta-analysis

*Preparing for meta-analysis*

-   The primary studies often use a range of statistics to present their
    results. Convenience functions to convert these onto a common metric
    are presented by: `r pkg("compute.es")` which converts
    from various statistics to d, g, r, z and the log odds ratio,
    `r pkg("MAc")` which converts to correlation
    coefficients, `r pkg("MAd")` which converts to mean
    differences, and `r pkg("metafor", priority = "core")`
    which converts to effect sizes an extensive set of measures for
    comparative studies (such as binary data, person years, mean
    differences and ratios and so on), for studies of association (a
    wide range of correlation types), for non-comparative studies
    (proportions, incidence rates, and mean change). It also provides
    for a measure used in psychometrics (Cronbach's alpha).
    `r pkg("esc")` provides a range of effect size
    calculations with partial overlap with
    `r pkg("metafor")` but with some extras, noticeably for
    converting test statistics, also includes a convenience function for
    collating its output for input to another package like
    `r pkg("metafor")` or producing a CSV file.
    `r pkg("estimraw")` estimates the cell frequencies from
    one of odds ratio, risk ratio, or risk difference.
    `r pkg("effsize")` contains functions to compute effect
    sizes mean difference (Cohen's d and Hedges g), dominance matrices
    (Cliff's Delta) and stochastic superiority (Vargha-Delaney A).
    `r pkg("effectsize")` provides a large number of
    different effect sizes and converts between them.
    `r pkg("psychmeta")` provides extensive facilties for
    converting effect sizes and for correcting for a variety of
    restrictions and measurement errors. `r pkg("metansue")`
    provides some methods for converting to effect sizes
    `r pkg("es.dif")` from raw data computes Cohen's d,
    Hedges' d, biased/unbiased c (an effect size between a mean and a
    constant) and e (an effect size between means without assuming the
    variance equality). `r pkg("MOTE")` provides a variety
    of conversions based on Cohen's d. `r pkg("estmeansd")`
    converts between quantiles and means and standard deviations.
    `r pkg("metaBLUE")` estimates means and standard
    deviations from various order statistics.
    `r pkg("SingleCaseES")` provides basic effect sizes for
    single-case designs, both parametric and non-overlap.
    `r pkg("smd")` computes standardised mean differences
-   `r pkg("meta", priority = "core")` provides functions to
    read and work with files output by RevMan 4 and 5.
-   `r pkg("metagear")` provides many tools for the
    systematic review process including screening articles, downloading
    the articles, generating a PRISMA diagram, and some tools for effect
    sizes. `r pkg("revtools")` provides tools for
    downloading from bibliographic databases and uses machine learning
    methods to process them.
-   `r pkg("metavcov")` computes the variance-covariance
    matrix for multivariate meta-analysis when correlations between
    outcomes can be provided but not between treatment effects
-   `r pkg("clubSandwich")` imputes variance-covariance
    matrix for multivariate meta-analysis
-   `r pkg("metafuse")` uses a fused lasso to merge
    covariate estimates across a number of independent datasets.
-   `r pkg("metapower")` provides power analysis for
    meta-analysis and meta-regression

*Fitting the model*

-   Four packages provide the inverse variance weighted,
    Mantel-Haenszel, and Peto methods: `r pkg("epiR")`,
    `r pkg("meta")`, `r pkg("metafor")`, and
    `r pkg("rmeta")`.
-   For binary data `r pkg("metafor")` provides the
    binomial-normal model.
-   For sparse binary data `r pkg("exactmeta")` provides an
    exact method which does not involve continuity corrections.
-   Packages which work with specific effect sizes may be more congenial
    to workers in some areas of science and include
    `r pkg("MAc")` and `r pkg("metacor")` which
    provide meta-analysis of correlation coefficients and
    `r pkg("MAd")` which provides meta-analysis of mean
    differences. `r pkg("MAc")` and
    `r pkg("MAd")` provide a range of graphics.
    `r pkg("psychometric")` provides an extensive range of
    functions for the meta-analysis of psychometric studies.
    `r pkg("mixmeta")` provides an integrated interface to
    standard meta-analysis and extensions like multivariate and
    dose-response.
-   `r pkg("psychmeta")` implements the Hunter-Schmidt
    method including corrections for reliability and range-restriction
    issues
-   `r pkg("concurve")` provides consonance curves relying
    on `r pkg("metafor")`
-   `r pkg("clubSandwich")` gives cluster-robust variance
    estimates.
-   Bayesian approaches are contained in various packages.
    `r pkg("bspmma")` which provides two different models: a
    non-parametric and a semi-parametric. Graphical display of the
    results is provided. `r pkg("mmeta")` provides
    meta-analysis using beta-binomial prior distributions.
    `r pkg("bayesmeta")` includes shrinkage estimates,
    posterior predictive p-values and forest plots via either
    `r pkg("metafor")` or `r pkg("forestplot")`.
    Diagnostic graphical output is available.
    `r pkg("MetaStan")` includes binomial-normal
    hierarchical models and can use weakly informative priors for the
    heterogeneity and treatment effect parameters.
    `r pkg("baggr")` provides facilities using Stan for
    hierarchical Bayesian models, graphical facilities are provided.
    `r pkg("BayesCombo")` provides facilities using a
    Bayesian approach and has graphical facilities.
    `r pkg("RBesT")` uses Bayesian synthesis to generate
    priors from various sources. `r pkg("metamisc")`
    provides a method with priors suggested by Higgins.
    `r pkg("RoBMA")` provides a framework for estimating
    ensembles of meta-analytic models and using Bayesian model averaging
    to combine them. `r pkg("ra4bayesmeta")` provides
    principled reference analysis within the Bayesian normal-normal
    model.
-   Some packages concentrate on providing a specialised version of the
    core meta-analysis function without providing the range of ancillary
    functions. These are: `r pkg("metaLik")` which uses a
    more sophisticated approach to the likelihood, and
    `r pkg("metatest")` which provides another improved
    method of obtaining confidence intervals.
    `r pkg("metaBMA")` has a Bayesian approach using model
    averaging, a variety of priors are provided and it is possible for
    the user to define new ones. `r pkg("gmeta")` which
    subsumes a very wide variety of models under the method of
    confidence distributions and also provides a graphical display,
    `r pkg("EvidenceSynthesis")` combines causal effect
    estimates and study diagnostics across studies.
    `r pkg("CoTiMA")` performs meta-analyses of correlation
    matrices of repeatedly measured variables for studies with different
    time lags using a SEM framework with OpenMx as the engine
-   `r pkg("metaplus")` fits random effects models relaxing
    the usual assumption that the random effects have a normal
    distribution by providing t or a mixture of normals.
-   `r pkg("ratesci")` fits random effects models to binary
    data using a variety of methods for confidence intervals.
-   `r pkg("RandMeta")` estimates exact confidence intervals
    in random effects models using an efficient algorithm.
-   `r pkg("rma.exact")` estimates exact confidence
    intervals in random effects normal-normal models and also provides
    plots of them.
-   `r pkg("pimeta")` provides a range of methods for
    prediction interval estimation from random effects models and ahs
    graphical facilities.
-   `r pkg("metamedian")` implements several methods to
    meta-analyze one-group or two-group studies that report the median
    of the outcome. These methods estimate the pooled median in the
    one-group context and the pooled raw difference of medians across
    groups in the two-group context `r pkg("meta")` also
    provides methods for medians
-   `r pkg("MetaUtility")` proposes a metric for estimating
    the proportion of effects above a cut-off of scientific importance
-   `r pkg("metasens")` provides imputation methods for
    missing binary data.
-   `r pkg("metagam")` provides a framework for
    meta-analysis of generalised additive models including the case
    where individual paticipant data cannot be shared across locations.
-   `r pkg("metawho")` implements a method for combining
    within study interactions
-   `r pkg("metarep")` provides replicability analyses after
    a conventional analysis

*Graphical methods*

An extensive range of graphical procedures is available.

-   Forest plots are provided in `r pkg("forplo")`,
    `r pkg("forestmodel")` (using ggplot2),
    `r pkg("forestplot")`, `r pkg("meta")`,
    `r pkg("metafor")`, `r pkg("metansue")`,
    `r pkg("psychmeta")`, and `r pkg("rmeta")`.
    Although the most basic plot can be produced by any of them they
    each provide their own choice of enhancements.
    `r pkg("metaviz")` provides a range of enhancements.
-   Funnel plots are provided in `r pkg("meta")`,
    `r pkg("metafor")`, `r pkg("metansue")`,
    `r pkg("psychometric")` `r pkg("rmeta")` and
    `r pkg("weightr")`. In addition to the standard funnel
    plots an enhanced funnel plot to assess the impact of extra evidence
    is available in `r pkg("extfunnel")`, a funnel plot for
    limit meta-analysis in `r pkg("metasens")`, and
    `r pkg("metaviz")` provides an extensive range of
    enhanced funnel plots and also facilities for their use in the
    context of visual inference.
-   Radial (Galbraith) plots are provided in `r pkg("meta")`
    and `r pkg("metafor")`.
-   L'Abbe plots are provided in `r pkg("meta")` and
    `r pkg("metafor")`.
-   Baujat plots are provided in `r pkg("meta")` and
    `r pkg("metafor")`.
-   `r pkg("meta")` provides drapery plots
-   `r pkg("metaplotr")` provides a crosshair plot
-   `r pkg("MetaAnalyser")` provides an interactive
    visualisation of the results of a meta-analysis.
-   `r pkg("metaviz")` provides rainforestplots, an enhanced
    version of forest plots. It accepts input from
    `r pkg("metafor")`.
-   `r pkg("DTAplots")` produces various plots for
    diagnostic studies including forest and SROC plots.

*Investigating heterogeneity*

-   Confidence intervals for the heterogeneity parameter are provided in
    `r pkg("metafor")` and `r pkg("psychmeta")`.
-   `r pkg("altmeta")` presents a variety of alternative
    methods for measuring and testing heterogeneity with a focus on
    robustness to outlying studies.
-   `r pkg("metaforest")` investigates heterogeneity using
    random forests. Note that it has nothing to do with forest plots.
-   `r pkg("mc.heterogeneity")` implements a Monte Carlo
    based test for heterogeneity.
-   `r pkg("boot.heterogeneity")` provides a bootstrp test
    for heterogeneity for mean differences, correlations, and odds
    ratios.

*Model criticism*

-   An extensive series of plots of diagnostic statistics is provided in
    `r pkg("metafor")`.
-   `r pkg("metaplus")` provides outlier diagnostics.
-   `r pkg("psychmeta")` provides leave-one-out methods.
-   `r pkg("EValue")` provides sensitivity analysis of the
    effect of unmeasured confounders
-   `r pkg("boutliers")` provides bootstrap distributions
    for outlier detection and influence diagnostics

*Investigating small study bias*

The issue of whether small studies give different results from large
studies has been addressed by visual examination of the funnel plots
mentioned above. In addition:

-   `r pkg("meta")` and `r pkg("metafor")`
    provide both the non-parametric method suggested by Begg and
    Mazumdar and a range of regression tests modelled after the approach
    of Egger.
-   `r pkg("xmeta")` provides a method in the context of
    multivariate meta-analysis.
-   `r pkg("metamisc")` provides funnel plots and tests for
    asymmetry.
-   An exploratory technique for detecting an excess of statistically
    significant studies is provided by `r pkg("PubBias")`.
-   `r pkg("puniform")` provides methods using only the
    statistically significant studies, methods for the special case of
    replication studies and sample size determinations.
-   `r pkg("PublicationBias")` performs sensitivity analysis
    of the number of unpublished studies needed to have a specified
    influence.
-   `r pkg("metafor")` provides a variety of selection
    models.

*Unobserved studies*

A recurrent issue in meta-analysis has been the problem of unobserved
studies.

-   Rosenthal's fail safe n is provided by `r pkg("MAc")`
    and `r pkg("MAd")`. `r pkg("metafor")`
    provides it as well as two more recent methods by Orwin and
    Rosenberg.
-   Duval's trim and fill method is provided by
    `r pkg("meta")` and `r pkg("metafor")`.
-   `r pkg("metasens")` provides Copas's selection model
    and also the method of limit meta-analysis (a regression based
    approach for dealing with small study effects) due to Rücker et al.
-   `r pkg("selectMeta")` provides various selection models:
    the parametric model of Iyengar and Greenhouse, the non-parametric
    model of Dear and Begg, and proposes a new non-parametric method
    imposing a monotonicity constraint.
-   `r pkg("SAMURAI")` performs a sensitivity analysis
    assuming the number of unobserved studies is known, perhaps from a
    trial registry, but not their outcome.
-   The `r pkg("metansue")` package allows the inclusion by
    multiple imputation of studies known only to have a non-significant
    result.
-   `r pkg("weightr")` provides facilities for using the
    weight function model of Vevea and Hedges.
-   `r pkg("publipha")` estimates models accounting for
    publication bias or p-hacking using a Bayesian framework
-   `r pkg("fsn")` computes the fail-safe number with
    confidence interval.
-   `r pkg("RobustBayesianCopas")` fits a robust version of
    the Copas selection model.
-   `r pkg("metafor")` provides the test of excess
    signifcance.

*Other study designs*

-   `r pkg("SCMA")` provides single case meta-analysis. It
    is part of a suite of packages dedicated to single-case designs.
-   `r pkg("joint.Cox")` provides facilities for the
    meta-analysis of studies of joint time-to-event and disease
    progression.
-   `r pkg("dfmeta")` provides meta-analysis of Phase I
    dose-finding clinical trials
-   `r pkg("metaRMST")` implements meta-analysis of trials
    with difference in restricted mean survival times
-   `r pkg("metamisc")` provides for meta-analysis of
    prognostic studies
-   `r pkg("metamicrobiomeR")` provides meta-analysis of
    zero-inflated beta microbiome data fitted with GAMLSS models.
-   `r pkg("metaSurvival")` estimates the survival curves
    from data extracted from primary study survival curves.

*Meta-analysis of significance values*

-   Fisher's method and Lancaster's are available in
    `r pkg("aggregation")`, `r pkg("metap")`,
    and `r pkg("poolr")`.
-   Stouffer's method, Tippett's and Wilkinson's are available in
    `r pkg("metap")` and `r pkg("poolr")`.
-   Edgington's method, inverse-t, logit, mean of p, and mean of z are
    all available in `r pkg("metap")`.

In all cases `r pkg("poolr")` considers correlated p-values
in addition to independent. The others above do not.

-   `r pkg("TFisher")` provides Fisher's method using both
    hard and soft thresholding for the p-values. There is a wrapper in
    `r pkg("metap")` for the hard threshold case.
-   `r pkg("harmonicmeanp")` uses the method of harmonic
    mean of p-values which is robust to correlation between the
    p-values.
-   `r pkg("amanida")` provides meta-analysis of metabolite
    data using p-values and fold change.
-   `r pkg("metap")` provides simple graphics.

Some methods are also provided in some of the genetics packages
mentioned below.

#### Multivariate meta-analysis

Standard methods outlined above assume that the effect sizes are
independent. This assumption may be violated in a number of ways: within
each primary study multiple treatments may be compared to the same
control, each primary study may report multiple endpoints, or primary
studies may be clustered for instance because they come from the same
country or the same research team. In these situations where the outcome
is multivariate:

-   `r pkg("mvmeta")` assumes the within study covariances
    are known and provides a variety of options for fitting random
    effects. `r pkg("metafor")` provides fixed effects and
    likelihood based random effects model fitting procedures. Both these
    packages include meta-regression, `r pkg("metafor")`
    also provides for clustered and hierarchical models.
-   `r pkg("mvtmeta")` provides multivariate meta-analysis
    using the method of moments for random effects although not
    meta-regression,
-   `r pkg("metaSEM")` provides multivariate (and
    univariate) meta-analysis and meta-regression by embedding it in the
    structural equation framework and using OpenMx for the structural
    equation modelling. It can provide a three-level meta-analysis
    taking account of clustering and allowing for level 2 and level 3
    heterogeneity. It also provides via a two-stage approach
    meta-analysis of correlation or covariance matrices.
-   `r pkg("xmeta")` provides various functions for
    multivariate meta-analysis and also for detecting publication bias.
-   `r pkg("dosresmeta")` concentrates on the situation
    where individual studies have information on the dose-response
    relationship. `r pkg("MBNMAdose")` provides a Bayesian
    analysis using network meta-analysis of dose response studies.
-   `r pkg("robumeta")` provides robust variance estimation
    for clustered and hierarchical estimates.
-   `r pkg("CIAAWconsensus")` has a function for
    multivariate m-a in the context of atomic weights and estimating
    isotope ratios.

#### Meta-analysis of studies of diagnostic tests

A special case of multivariate meta-analysis is the case of summarising
studies of diagnostic tests. This gives rise to a bivariate, binary
meta-analysis with the within-study correlation assumed zero although
the between-study correlation is estimated. This is an active area of
research and a variety of methods are available including what is
referred to here as Reitsma's method, and the hierarchical summary
receiver operating characteristic (HSROC) method. In many situations
these are equivalent.

-   `r pkg("mada")` provides various descriptive statistics
    and univariate methods (diagnostic odds ratio and Lehman model) as
    well as the bivariate method due to Reitsma. Meta-regression is
    provided. Graphical facilities are also available.
-   `r pkg("Metatron")` provides a method for the Reitsma
    model incuding the case of an imperfect reference standard.
-   `r pkg("bamdit")` provides Bayesian meta-analysis with a
    bivariate random effects model (using JAGS to implement the MCMC
    method). Graphical methods are provided.
-   `r pkg("meta4diag")` provides Bayesian inference
    analysis for bivariate meta-analysis of diagnostic test studies and
    an extensive range of graphical methods.
-   `r pkg("CopulaREMADA")` uses a copula based mixed model
-   `r pkg("diagmeta")` considers the case where the primary
    studies provide analysis using multiple cut-offs. Graphical methods
    are also provided.
-   `r pkg("CopulaDTA")` uses the beta-binomial model to
    yield marginal mean sensitivity and specificity. Graphical
    facilities are available.
-   `r pkg("NMADiagT")` provides network meta-analysis of
    diagnostic tests in a Bayesian framework using Stan as the engine,
    graphical output is provided.
-   `r pkg("DTAplots")` prouces various plots for diagnostic
    studies including forest and SROC plots.

#### Meta-regression

Where suitable moderator variables are available they may be included
using meta-regression. All these packages are mentioned above, this just
draws that information together.

-   `r pkg("metafor")` provides meta-regression (multiple
    moderators are catered for). Various packages rely on
    `r pkg("metafor")` to provide meta-regression
    (`r pkg("meta")`, `r pkg("MAc")`, and
    `r pkg("MAd")`) and all three of these provide bubble
    plots. `r pkg("psychmeta")` also uses
    `r pkg("metafor")`.
-   `r pkg("metaLik")`, `r pkg("metansue")`,
    `r pkg("metaSEM")`, and `r pkg("metatest")`
    also provide meta-regression.
-   `r pkg("mvmeta")` provides meta-regression for
    multivariate meta-analysis as do `r pkg("metafor")` and
    `r pkg("metaSEM")`.
-   `r pkg("mada")` provides for the meta-regression of
    diagnostic test studies.
-   `r pkg("GENMETA")` uses generalised meta-analysis to
    handle the situation where the studies do not all use the same
    regressors
-   `r pkg("jarbes")` uses the Bayesian approach of
    hierarchical meta-regression
-   `r pkg("metacart")` uses classification and regression
    trees to identify interactions between moderators

#### Individual participant data (IPD)

Where all studies can provide individual participant data then software
for analysis of multi-centre trials or multi-centre cohort studies
should prove adequate and is outside the scope of this task view. Other
packages which provide facilities related to IPD are:

-   `r pkg("ecoreg")` which is designed for ecological
    studies enables estimation of an individual level logistic
    regression from aggregate data or individual data.
-   `r pkg("multinma")` provides network meta-analysis and
    network meta-regression models for aggregate data, individual
    patient data, and mixtures of both individual and aggregate data
-   `r pkg("MetaIntegration")` combines IPD data with
    external models.

#### Network meta-analysis

Also known as multiple treatment comparison. This is a very active area
of research and development. Note that some of the packages mentioned
above under multivariate meta-analysis can also be used for network
meta-analysis with appropriate setup.

-   `r pkg("netmeta")` works in a frequentist framework. It
    provides an extensive range of graphical and other displays
    including network graphs and a heatmap for displaying inconsistency
    and heterogeneity. A frequentist analogue of SUCRA is also
    available.
-   A Bayesian framework by `r pkg("pcnetmeta")`, which uses
    JAGS. It provides a number of data-sets.
    `r pkg("nmaINLA")` uses integrated nested Laplace
    approximations as an alternative to MCMC. It provides a number of
    data-sets. `r pkg("NMADiagT")` provides network
    meta-analysis of diagnostic tests in a Bayesian framework using Stan
    as the engine, graphical output is provided.
    `r pkg("gemtc")`, which acts as a front-end to BUGS or
    JAGS, `r pkg("bnma")` provides arm-based methods using
    JAGS as the engine, `r pkg("metapack")` provides methods
    using built-in MCMC code
-   `r pkg("multinma")` provides network meta-analysis and
    network meta-regression models for aggregate data, individual
    patient data, and mixtures of both individual and aggregate data
-   `r pkg("nmathresh")` provides decision-invariant bias
    adjustment thresholds and intervals the smallest changes to the data
    that would result in a change of decision.
    `r pkg("NMAoutlier")` detects outliers in NMA using
    forward search,
-   `r pkg("pcnetmeta")` provides network graphs.
    `r pkg("nmaplateplot")` displays the results from an NMA
    using a heatplot style and also displays SUCRA.
-   `r pkg("nmarank")` evaluates hierarchies of evidence in
    network meta-analysis

#### Genetics

There are a number of packages specialising in genetic data:
`r pkg("catmap")` combines case-control and family study
data, graphical facilities are provided, `r pkg("CPBayes")`
uses a Bayesian approach to study cross-phenotype genetic associations,
`r pkg("etma")` proposes a new statistical method to detect
epistasis, `r pkg("gap")` combines p-values,
`r pkg("getmstatistic")` quantifies systematic
heterogeneity, `r pkg("getspres")` uses standardised
predictive random effects to explore heterogeneity in genetic
association meta-analyses, `r pkg("GMCM")` uses a Gaussian
mixture copula model for high-throughput experiments,
`r pkg("MBNMAtime")` provides methods for analysis of
repeated measures network meta-analysis,
`r pkg("MendelianRandomization")` provides several methods
for performing Mendelian randomisation analyses with summarised data,
`r pkg("MetABEL")` provides meta-analysis of genome wide SNP
association results, `r pkg("MetaIntegrator")` provides
meta-analysis of gene expression data, `r pkg("metaMA")`
provides meta-analysis of p-values or moderated effect sizes to find
differentially expressed genes, `r pkg("MetaPath")` performs
meta-analysis for pathway enrichment, `r pkg("metaRNASeq")`
meta-analysis from multiple RNA sequencing experiments,
`r pkg("MetaSubtract")` uses leave-one-out methods to
validate meta-GWAS results, `r pkg("ofGEM")` provides a
method for identifying gene-environment interactions using
meta-filtering, `r pkg("RobustRankAggreg")` provides methods
for aggregating lists of genes, `r pkg("SPAtest")` combines
association results.

#### Data-sets

-   `r pkg("metadat")` provides a large number of data-sets
    used in meta-analysis
-   `r pkg("nmadb")` provides access to a database of
    network meta-analyses

#### Interfaces

-   Plug-ins for Rcmdr are provided by:
    `r pkg("RcmdrPlugin.EZR")` which uses
    `r pkg("meta")` and `r pkg("metatest")`,
    `r pkg("RcmdrPlugin.MA")` which uses
    `r pkg("MAd")` and `r pkg("metafor")`, and
    `r pkg("RcmdrPlugin.RMTCJags")` for network
    meta-analysis using BUGS code.

#### Simulation

`r pkg("psychmeta")` provides facilities for simulation of
psychometric data-sets.

#### Others

`r pkg("CRTSize")` provides meta-analysis as part of a
package primarily dedicated to the determination of sample size in
cluster randomised trials in particular by simulating adding a new study
to the meta-analysis.

`r pkg("CAMAN")` offers the possibility of using finite
semiparametric mixtures as an alternative to the random effects model
where there is heterogeneity. Covariates can be included to provide
meta-regression.

`r pkg("KenSyn")` provides data-sets to accompany a French
language book on meta-analysis in the agricultural sciences.

`r pkg("PRISMAstatement")` generates a flowchart conforming
to the PRISMA statement.

`r pkg("metabolic")` provides data and code to support a book.

