# Technical Overview

## Introduction

Geostatistics refers to the study of a spatially distributed variable of
interest, which in theory is defined at every point over a bounded study
region of interest. Statistical modelling and analysis for spatially
oriented point-referenced outcomes play a crucial role in diverse
scientific applications such as earth and environmental sciences,
ecology, epidemiology, and economics. With the advent of Markov chain
Monte Carlo (MCMC) algorithms, Bayesian hierarchical models have gained
massive popularity in analyzing such point-referenced or, geostatistical
data. These models involve latent spatial processes characterized by
spatial process parameters, which besides lacking substantive relevance
in scientific contexts, are also weakly identified and hence, impedes
convergence of MCMC algorithms. Thus, even for moderately large datasets
(~$10^{3}$ or higher), the computation for MCMC becomes too onerous for
practical use.

We introduce the R package `spStack` that implements Bayesian inference
for a class of geostatistical models, where we obviate the issues
mentioned by sampling from analytically available posterior
distributions conditional upon some candidate values of the spatial
process parameters and, subsequently assimilate inference from these
individual posterior distributions using Bayesian predictive stacking.
Besides delivering competitive predictive performance as compared to
fully Bayesian inference using MCMC, our proposed algorithm is
embarrassingly parallel, thus drastically improves runtime and elevating
the utility of the package for a diverse group of practitioners with
limited computational resources at their disposal. This package, to the
best of our knowledge, is the first to implement stacking for Bayesian
analysis of spatial data.

Technical details surrounding the methodology can be found in the
articles Zhang, Tang, and Banerjee ([2025](#ref-zhang2024stacking))
which discuss the case where the distribution of the point-referenced
outcomes are Gaussian, and, in Pan et al. ([2025](#ref-pan2024stacking))
where the case of non-Gaussian outcomes is explored. The code for this
package is written primarily in C/C++ with additional calls to FORTRAN
routines for optimized linear algebra operations. We leverage the
`F77_NAME` macro to interface with legacy FORTRAN functions in
conjunction with efficient matrix computation libraries such as
[BLAS](https://netlib.org/blas/) (Basic Linear Algebra Subprograms) and
[LAPACK](https://netlib.org/lapack/) (Linear Algebra Package) to
implement our stacking algorithm.

The remainder of the vignette evolves as follows - the next two sections
discuss Bayesian hierarchical spatial models for Gaussian and
non-Gaussian outcomes, which is followed by a section providing brief
details on predictive stacking and a section dedicated for illustration
of functions in the package.

## Bayesian Gaussian spatial regression models

Let $\chi = \{ s_{1},\ldots,s_{n}\} \in \mathcal{D}$ be a be a set of
$n$ spatial locations yielding measurements
$y = \left( y_{1},\ldots,y_{n} \right)^{\top}$ with known values of
predictors at these locations collected in the $n \times p$ full rank
matrix
$X = \left\lbrack x\left( s_{1} \right),\ldots,x\left( s_{n} \right) \right\rbrack^{\top}$.
A customary geostatistical model is
$$y_{i} = x\left( s_{i} \right)^{\top}\beta + z\left( s_{i} \right) + \epsilon_{i},\quad i = 1,\ldots,n,$$
where $\beta$ is the $p \times 1$ vector of slopes,
$z(s) \sim {\mathsf{G}\mathsf{P}}\left( 0,R\left( \cdot , \cdot ;\theta_{\text{sp}} \right) \right)$
is a zero-centered spatial Gaussian process on $\mathcal{D}$ with
spatial correlation function
$R\left( \cdot , \cdot ;\theta_{\text{sp}} \right)$ characterized by
process parameters $\theta_{\text{sp}}$, $\sigma^{2}$ is the spatial
variance parameter (“partial sill”) and
$\epsilon_{i} \sim \mathsf{N}\left( 0,\tau^{2} \right),i = 1,\ldots,n$
are i.i.d. with variance $\tau^{2}$ (“nugget”) capturing measurement
error. The spatial process $z( \cdot )$ is assumed to be independent of
the measurement errors $\{\epsilon_{i},i = 1,\ldots,n\}$. Let
$z = \left( z\left( s_{1} \right),\ldots,z\left( s_{n} \right) \right)^{\top}$
denotes the realization of the spatial process on $\chi$ and
$n \times n$ correlation matrix
$R\left( \chi;\theta_{\text{sp}} \right) = \left( R\left( s_{i},s_{j}\theta_{\text{sp}} \right) \right)_{1 \leq i,j \leq n}$.
We build a conjugate Bayesian hierarchical spatial model,
$$\begin{array}{r}
\begin{aligned}
{y \mid z,\beta,\sigma^{2}} & {\sim \mathsf{N}\left( X\beta + z,\delta^{2}\sigma^{2}I_{n} \right),} \\
{z \mid \sigma^{2}} & {\sim \mathsf{N}\left( 0,\sigma^{2}R\left( \chi;\theta_{\text{sp}} \right) \right),} \\
{\beta \mid \sigma^{2}} & {\sim \mathsf{N}\left( \mu_{\beta},\sigma^{2}V_{\beta} \right),\quad\sigma^{2} \sim {\mathsf{I}\mathsf{G}}\left( a_{\sigma},b_{\sigma} \right),}
\end{aligned}
\end{array}$$

where we fix the noise-to-spatial variance ratio
$\delta^{2} = \tau^{2}/\sigma^{2}$, the process parameters
$\theta_{\text{sp}}$ and the hyperparameters $\mu_{\beta}$, $V_{\beta}$,
$a_{\sigma}$ and $b_{\sigma}$. In this package, we use the Matern
covariogram specified by spatial decay parameter $\phi$ and smoothness
parameter $\nu$ i.e., $\theta_{\text{sp}} = \{\phi,\nu\}$, given by
$$R\left( s,s\prime;\theta_{\text{sp}} \right) = \frac{\left( \phi|s - s\prime| \right)^{\nu}}{2^{\nu - 1}\Gamma(\nu)}K_{\nu}\left( \phi|s - s\prime| \right)).$$
We utilize a composition sampling strategy to sample the model
parameters from their joint posterior distribution which can be written
as
$$p\left( \sigma^{2},\beta,z \mid y \right) = p\left( \sigma^{2} \mid y \right) \times p\left( \beta \mid \sigma^{2},y \right) \times p\left( z \mid \beta,\sigma^{2},y \right).$$
We proceed by first sampling $\sigma^{2}$ from its marginal posterior,
then given the samples of $\sigma^{2}$, we sample $\beta$ and
subsequently, we sample $z$ conditioned on the posterior samples of
$\beta$ and $\sigma^{2}$([Banerjee 2020](#ref-banerjee_massivespatial)).
More details can be found in Zhang, Tang, and Banerjee
([2025](#ref-zhang2024stacking)).

The function
[`spLMexact()`](https://span-18.github.io/spStack-dev/reference/spLMexact.md)
delivers samples from this posterior distribution conditional on fixed
hyperparameters. For predictive stacking, use the function
[`spLMstack()`](https://span-18.github.io/spStack-dev/reference/spLMstack.md).

## Bayesian non-Gaussian spatial regression models

Analyzing non-Gaussian spatial data typically requires introducing
spatial dependence in generalized linear models through the link
function of an exponential family distribution. Let $y(s)$ be the
outcome at location $s \in \mathcal{D}$ endowed with a probability law
from the natural exponential family, which we denote by
$$y(s) \sim {\mathsf{E}\mathsf{F}}\left( x(s)^{\top}\beta + z(s);b,\psi_{y} \right)$$
for some positive parameter $b > 0$ and unit log partition function
$\psi_{y}$. Fixed effects regression and spatial dependence, e.g.,
$x(s)^{\top}\beta + z(s)$, is introduced in the natural parameter, where
$x(s)$ is a $p \times 1$ vector of predictors referenced with respect to
$s$, $\beta$ is a $p \times 1$ vector of slopes measuring the trend,
$z(s)$ is a zero-centered spatial process on $\mathcal{D}$ specified by
a scale parameter $\sigma_{z}$ and a spatial correlation function
$R\left( \cdot , \cdot ;\theta_{\text{sp}} \right)$ with
$\theta_{\text{sp}}$ consisting of spatial-temporal decay and smoothness
parameters.

Unlike in Gaussian likelihoods, inference is considerably encumbered by
the inability to analytically integrate out the random effects and
reduce the dimension of the parameter space. Iterative algorithms such
as Markov Chain Monte Carlo (MCMC), thus attempt to sample from a very
high-dimensional posterior distribution, and convergence is often
hampered by high auto-correlations and weakly identified spatial process
parameters $\theta_{\text{sp}}$.

This model is implemented using the function
[`spGLMexact()`](https://span-18.github.io/spStack-dev/reference/spGLMexact.md)
when using fixed hyperparameters, and
[`spGLMstack()`](https://span-18.github.io/spStack-dev/reference/spGLMstack.md)
when using predictive stacking.

We consider the following three point-referenced data -

- **Poisson count data**: Here $b = 1$ and $\psi_{y}(t) = e^{t}$.
  $$\begin{aligned}
  {y\left( s_{i} \right)} & {\sim {\mathsf{P}\mathsf{o}\mathsf{i}\mathsf{s}\mathsf{s}\mathsf{o}\mathsf{n}}\left( \lambda\left( s_{i} \right) \right),\quad i = 1,\ldots,n.} \\
  {\lambda\left( s_{i} \right)} & {= \exp\left( x\left( s_{i} \right)^{\top}\beta + z\left( s_{i} \right) \right)}
  \end{aligned}$$ This is accessed by setting `family = "poisson"` in
  the above functions.

- **Binomial count data**: Here $b = m\left( s_{i} \right)$ for each $i$
  and $\psi_{y}(t) = \log\left( 1 + e^{t} \right)$. $$\begin{aligned}
  {y\left( s_{i} \right)} & {\sim {\mathsf{B}\mathsf{i}\mathsf{n}\mathsf{o}\mathsf{m}\mathsf{i}\mathsf{a}\mathsf{l}}\left( m\left( s_{i} \right),\pi\left( s_{i} \right) \right),\quad i = 1,\ldots,n.} \\
  {\pi\left( s_{i} \right)} & {= {ilogit}\left( x\left( s_{i} \right)^{\top}\beta + z\left( s_{i} \right) \right)}
  \end{aligned}$$ This is accessed by setting `family = "binomial"` in
  the above functions.

- **Binary data**: Here $b = 1$ and
  $\psi_{y}(t) = \log\left( 1 + e^{t} \right)$. $$\begin{aligned}
  {y\left( s_{i} \right)} & {\sim {\mathsf{B}\mathsf{e}\mathsf{r}\mathsf{n}\mathsf{o}\mathsf{u}\mathsf{l}\mathsf{l}\mathsf{i}}\left( \pi\left( s_{i} \right) \right),\quad i = 1,\ldots,n.} \\
  {\pi\left( s_{i} \right)} & {= {ilogit}\left( x\left( s_{i} \right)^{\top}\beta + z\left( s_{i} \right) \right)}
  \end{aligned}$$ This is accessed by setting `family = "binary"` in the
  above functions.

Following Bradley and Clinch ([2024](#ref-bradleyclinch2024)), we
introduce a Bayesian hierarchical spatial model as $$\begin{aligned}
{y\left( s_{i} \right) \mid \beta,z,\xi} & {\sim {\mathsf{E}\mathsf{F}}\left( x\left( s_{i} \right)^{\top}\beta + z\left( s_{i} \right) + \xi_{i} - \mu_{i};b_{i},\psi_{y} \right),i = 1,\ldots,n} \\
{\beta \mid \sigma_{\beta}^{2}} & {\sim \mathsf{N}\left( 0,\sigma_{\beta}^{2}V_{\beta} \right),\quad\sigma_{\beta}^{2} \sim {\mathsf{I}\mathsf{G}}\left( \nu_{\beta}/2,\nu_{\beta}/2 \right)} \\
{z \mid \sigma_{z}^{2}} & {\sim \mathsf{N}\left( 0,\sigma_{z}^{2}R\left( \chi;\theta_{\text{sp}} \right) \right),\quad\sigma_{z}^{2} \sim {\mathsf{I}\mathsf{G}}\left( \nu_{z}/2,\nu_{z}/2 \right),} \\
{\xi \mid \beta,z,\sigma_{\xi}^{2},\alpha_{\epsilon}} & {\sim {\mathsf{G}\mathsf{C}\mathsf{M}_{\mathsf{c}}}\left( {\widetilde{\mu}}_{\xi},H_{\xi},\epsilon,\kappa_{\xi};\psi_{\xi} \right),}
\end{aligned}$$ where
$\mu = \left( \mu_{1},\ldots,\mu_{n} \right)^{\top}$ denotes the
discrepancy parameter. We fix the spatial process parameters
$\theta_{\text{sp}}$, the boundary adjustment parameter $\epsilon$ and
the hyperparameters $V_{\beta}$, $\nu_{\beta}$, $\nu_{z}$ and
$\sigma_{\xi}^{2}$. The term $\xi$ is known as the fine-scale variation
term which is given a conditional generalized conjugate multivariate
distribution ($GCM_{c}$) as prior. For details, see Pan et al.
([2025](#ref-pan2024stacking)).

## Bayesian non-Gaussian spatial-temporal regression model

We consider a rich family of Bayesian spatial-temporal model with
spatially-temporally varying regression coefficients. Suppose
$\ell = (s,t)$, with location $s \in \mathcal{D}$ and time
$t \in \mathcal{T}$, denote a spatial-temporal coordinate in
$\mathcal{L} = \mathcal{D} \times \mathcal{T}$.

Let $\mathcal{L} = \{\ell_{1},\ldots,\ell_{n}\}$ be a fixed set of $n$
distinct space-time coordinates in $\mathcal{D}$, where
$y(\mathcal{L}) = \left( y\left( \ell_{1} \right),\ldots,y\left( \ell_{n} \right) \right)^{\top}$,
which we simply denote by $y$, is the vector of observed outcomes, each
distributed as a member of the natural exponential family with log
partition function $\psi_{y}$. Suppose, $x\left( \ell_{i} \right)$ is a
$p \times 1$ vector of predictors, $\beta$ is the corresponding
$p \times 1$ vector of slopes (fixed effects),
$\widetilde{x}\left( \ell_{i} \right)$ is $r \times 1$ ($r \leq p$)
consisting of predictors in $x\left( \ell_{i} \right)$ that are posited
to have spatially-temporally varying regression coefficients
$z\left( \ell_{i} \right) = \left( z_{1}\left( \ell_{i} \right),\ldots,z_{r}\left( \ell_{i} \right) \right)^{\top}$,
where each $z_{j}\left( \ell_{i} \right)$ is a spatially-temporally
varying coefficient for the predictor
${\widetilde{x}}_{j}\left( \ell_{i} \right)$, $\xi_{i}$ is a fine-scale
variation term and $\mu_{i}$ is the discrepancy parameter (see above).
We introduce spatially-temporally varying coefficients in
$\eta\left( \ell \right)$ as $$\begin{aligned}
{y\left( \ell_{i} \right)} & {\mid \beta,z\left( \ell_{i} \right),\xi_{i},\mu_{i}\overset{\text{ind}}{\sim}{\mathsf{E}\mathsf{F}}\left( \eta\left( \ell_{i} \right) + \xi_{i} - \mu_{i};b_{i},\psi_{y} \right),\ i = 1,\ldots,n\;,} \\
{\eta\left( \ell \right)} & {= x\left( \ell \right)^{\top}\beta + \widetilde{x}\left( \ell \right)^{\top}z\left( \ell \right),\quad\beta \mid \sigma_{\beta}^{2},\mu_{\beta},V_{\beta} \sim \mathsf{N}\left( \mu_{\beta},\sigma_{\beta}^{2}V_{\beta} \right),\quad\sigma_{\beta}^{2} \sim \pi_{\beta}\left( \sigma_{\beta}^{2} \right)\;,} \\
{z\left( \ell \right)} & {\mid \theta_{z},\theta_{\text{sp}} \sim {\mathsf{G}\mathsf{P}}\left( 0,C_{z}\left( \cdot , \cdot ;\theta_{\text{sp}},\theta_{z} \right) \right)\;,\quad\theta_{z} \sim \pi_{z}\left( \theta_{z} \right)\;,} \\
\xi & {\mid \beta,z,\mu,\alpha_{\epsilon},\kappa_{\epsilon},\sigma_{\xi}^{2} \sim {\mathsf{G}\mathsf{C}\mathsf{M}_{\mathsf{c}}}\left( {\widetilde{\mu}}_{\xi},H_{\xi},\alpha_{\xi},\kappa_{\xi},D_{\xi},\pi_{\xi};\psi_{\xi} \right),\ \sigma_{\xi}^{2} \sim \pi_{\xi}\left( \sigma_{\xi}^{2} \right),\ p(\mu) \propto 1\;,}
\end{aligned}$$ where
$z\left( \ell \right) = \left( z_{1}\left( \ell \right),\ldots,z_{r}\left( \ell \right) \right)^{\top}$
is a multivariate Gaussian process with a separable cross-covariance
function
$C_{z}\left( \cdot , \cdot ;\theta_{\text{sp}},\theta_{z} \right)$,
characterized by process parameters $\theta_{\text{sp}}$ which controls
the within-process spatial-temporal correlation, and $\theta_{z}$ which
controls the between-process covariance matrix . Given $\theta_{z}$, the
$nr \times 1$ vector
$z = \left( z_{1}^{\top},\ldots,z_{r}^{\top} \right)^{\top}$, where
$z_{j} = \left( z_{j}\left( \ell_{1} \right),\ldots,z_{j}\left( \ell_{n} \right) \right)^{\top}$
for $j = 1,\ldots,r$, follows a multivariate Gaussian distribution with
mean $0$ and $nr \times nr$ covariance matrix
$C_{z}\left( \mathcal{L};\theta_{\text{sp}},\theta_{z} \right)$.

This model is implemented by the function
[`stvcGLMexact()`](https://span-18.github.io/spStack-dev/reference/stvcGLMexact.md)
under fixed hyperparameters, and
[`stvcGLMstack()`](https://span-18.github.io/spStack-dev/reference/stvcGLMstack.md)
when using predictive stacking. We also implement different
specifications for $C_{z}$ in this package as follows.

1.  **Independent spatial-temporal process**: We consider $r$ Gaussian
    spatial-temporal processes $$\begin{aligned}
    {{\text{Independent processes:}\mspace{6mu}}z_{j}\left( \ell \right) \mid \sigma_{z_{j}}^{2},{\theta_{\text{sp}}}_{j}} & {\overset{\text{ind}}{\sim}{\mathsf{G}\mathsf{P}}\left( 0,\sigma_{z_{j}}^{2}R_{j}\left( \cdot , \cdot ;{\theta_{\text{sp}}}_{j} \right) \right),} \\
    \sigma_{z_{j}}^{2} & {\sim {\mathsf{I}\mathsf{G}}\left( \nu_{z_{j}}/2,\nu_{z_{j}}/2 \right),\quad j = 1,\ldots,r,}
    \end{aligned}$$ where $\sigma_{z_{j}}^{2}$ is the variance parameter
    corresponding to process $z_{j}\left( \ell \right)$. This
    corresponds to the covariance matrix
    $C_{z}\left( \mathcal{L};\theta_{\text{sp}},\theta_{z} \right) = \oplus_{j = 1}^{r}\sigma_{z_{j}}^{2}R_{j}\left( {\theta_{\text{sp}}}_{j} \right)$
    with
    $\theta_{\text{sp}} = \{{\theta_{\text{sp}}}_{j}:j = 1,\ldots,r\}$,
    where ${\theta_{\text{sp}}}_{j}$ denotes covariance kernel
    parameters for the $j$th process, and
    $\theta_{z} = \{\sigma_{z_{1}}^{2},\ldots,\sigma_{z_{r}}^{2}\}$.
    This is accessed by setting the option
    `process.type = "independent"` in the above functions.

2.  **Independent shared spatial-temporal process**: This corresponds to
    the above with ${\theta_{\text{sp}}}_{j} = \theta_{\text{sp}}$ and
    $\sigma_{zj}^{2} = \sigma_{z}^{2}$ for all $j = 1,\ldots,r$. This is
    accessed by setting the option `process.type = "independent.shared"`
    in the above functions.

3.  **Multivariate spatial-temporal process**: We can introduce
    dependence among the elements of the $r \times 1$ vector
    $z\left( \ell \right)$ using
    $${\text{Multivariate process:}\mspace{6mu}}z\left( \ell \right) \mid \Sigma \sim {\mathsf{G}\mathsf{P}}\left( 0,R\left( \cdot , \cdot ;\theta_{\text{sp}} \right)\Sigma \right),\quad\Sigma \sim {\mathsf{I}\mathsf{W}}\left( \nu_{z} + 2r,\Psi \right)\;,$$
    where
    ${\mathcal{G}\mathcal{P}}\left( 0,R\left( \cdot , \cdot ;\theta_{\text{sp}} \right)\Sigma \right)$
    is an $r \times 1$ multivariate Gaussian process with matrix-valued
    cross-covariance function
    $R\left( \cdot , \cdot ;\theta_{\text{sp}} \right)\Sigma$ and
    $\Sigma$ is an $r \times r$ positive definite random matrix. This
    corresponds to the spatial-temporal covariance matrix
    $C_{z}\left( \mathcal{L};\theta_{\text{sp}},\theta_{z} \right) = \Sigma \otimes R\left( \theta_{\text{sp}} \right)$
    with $\theta_{z} = \Sigma$. We place an inverse-Wishart prior on the
    scale parameter with shape $\nu_{z} + 2r$ and $r \times r$ positive
    definite scale matrix $\Psi$, given by
    $\pi\left( \theta_{z} \right) = {\mathsf{I}\mathsf{W}}\left( \Sigma \mid \nu_{z} + 2r,\Psi \right)$.
    This is accessed by setting the option
    `process.type = "multivariate"` in the above functions.

## Predictive stacking

Following Yao et al. ([2018](#ref-yao2018stacking)), we consider a set
of candidate models based on a grid of values of the parameters in
$\{\theta_{\text{sp}},\delta^{2}\}$ for the Gaussian case, and
$\{\theta_{\text{sp}},\epsilon\}$ for the non-Gaussian case, as will be
supplied by the user. We build a set of candidate models based on the
Cartesian product of the collection of values for each individual
parameter as $\mathcal{M} = \{ M_{1},\ldots,M_{G}\}$. Then, for each
$g = 1,\ldots,G$, we sample from the posterior distribution
$p\left( \sigma^{2},\beta,z \mid y,M_{g} \right)$ under the model
$M_{g}$ and find leave-one-out predictive densities
$p\left( y_{i} \mid y_{- i},M_{g} \right)$. Then we solve the
optimization problem $$\begin{aligned}
\max\limits_{w_{1},\ldots,w_{G}} & {\,\frac{1}{n}\sum\limits_{i = 1}^{n}\log\sum\limits_{g = 1}^{G}w_{g}p\left( y_{i} \mid y_{- i},M_{g} \right)} \\
\text{subject to} & {\quad w_{g} \geq 0,\sum\limits_{g = 1}^{G}w_{g} = 1}
\end{aligned}$$ to find the optimal stacking weights
${\widehat{w}}_{1},\ldots,{\widehat{w}}_{G}$. After obtaining the
optimal stacking weights, posterior inference of any quantity of
interest subsequently proceed from the *stacked* posterior,
$$\widetilde{p}( \cdot \mid y) = \sum\limits_{g = 1}^{G}{\widehat{w}}_{g}p\left( \cdot \mid y,M_{g} \right).$$

## References

Banerjee, Sudipto. 2020. “Modeling Massive Spatial Datasets Using a
Conjugate Bayesian Linear Modeling Framework.” *Spatial Statistics* 37:
100417. <https://doi.org/10.1016/j.spasta.2020.100417>.

Bradley, Jonathan R., and Madelyn Clinch. 2024. “Generating Independent
Replicates Directly from the Posterior Distribution for a Class of
Spatial Hierarchical Models.” *Journal of Computational and Graphical
Statistics* 0 (0): 1–17.
<https://doi.org/10.1080/10618600.2024.2365728>.

Pan, Soumyakanti, Lu Zhang, Jonathan R. Bradley, and Sudipto Banerjee.
2025. “Bayesian Inference for Spatial-Temporal Non-Gaussian Data Using
Predictive Stacking.” *Bayesian Analysis* (In press).
<https://doi.org/10.1214/25-BA1582>.

Yao, Yuling, Aki Vehtari, Daniel Simpson, and Andrew Gelman. 2018.
“Using Stacking to Average Bayesian Predictive Distributions (with
Discussion).” *Bayesian Analysis* 13 (3): 917–1007.
<https://doi.org/10.1214/17-BA1091>.

Zhang, Lu, Wenpin Tang, and Sudipto Banerjee. 2025. “Bayesian
Geostatistics Using Predictive Stacking.” *Journal of the American
Statistical Association* (In press) (January).
<https://doi.org/10.1080/01621459.2025.2566449>.
