---
title: "spStack"
---

## A brief introduction

The R package spStack delivers Bayesian inference for point-referenced spatial data by assimilating posterior inference over a collection of candidate models using stacking of predictive densities. Currently, it supports point-referenced Gaussian, Poisson, binomial and binary outcomes. Users can supply candidate values of spatial process parameters and certain auxiliary model parameters, based on which the collection of models will be created. spStack utilizes the Bayesian conjugate linear modelling framework for Gaussian data and the generalized conjugate multivariate distribution theory for non-Gaussian exponential family data. Technical details are available in [Zhang, Tang, and Banerjee 2024](doi:10.48550/arXiv.2304.12414) and [Pan, Zhang, Bradley, and, Banerjee 2024](doi:10.48550/arXiv.2406.04655).
