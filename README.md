# ParamHermite (A Maple library for Parametric Hermite matrices) 
### (in testing)

## About

ParamHermite.mla is a Maple library for real root classification using parametric Hermite matrices.  
The software is under development within the scope of my PhD at PolSys team - LIP6 - Sorbonne Universite.

Given a system of polynomial equations `F` in variables `vars` and parameters `params`,
a parametric Hermite matrix `H` is a matrix whose coefficients are rational functions in `params` with the following properties

* The rank of `H` when `params` is instantiated to generic complex values equals the number of complex solutions of `F`
* The signature of `H` when `params` is instantiated to generic real values equals the number of real solutions of `F`

## Installation

### Dependencies

It requires the following software to function correctly
* [FGb](https://www-polsys.lip6.fr/~jcf/FGb/) (v. 1.70) for Grobner basis computations.
* [Maple 2018+](https://www.maplesoft.com/)

### Setup for Maple

The library is provided as a file `ParamHermite.mla` which should be loaded in Maple.

1. Add the following line in your mapleinit file (usually found at $HOME/.mapleinit)

```
libname := "PATH/TO/ParamHermite.mla", libname:
```

2. In each Maple session, you can load the library by

```
with(ParamHermite):
```

## Features

This library allows to compute a similar output as the function DiscriminantVariety in RootFinding:-Parametric library in Maple in a different (faster) way. It gives a list of polynomials whose sets of solutions form a boundary that separates the space of `params` into regions over each of which the number of real solutions in `vars` of `F` is constant. Then, this library computes at least one point for each region and semi-algebraic defining formulas of those regions.

(An example will be added soon)

## Main functions

Currently, the library supports only a system of polynomial equations (denoted by `F`) in a list of variables (`vars`) and parameters (`params`).

Note: The function for handling one inequality exists but only for testing.

Here is a list of main functions:
* Compute a parametric Hermite matrix w.r.t. DRL ordering:
```
DRLMatrix(F,vars,params,char):
```
* Compute a parametric Hermite matrix w.r.t. DRL ordering with some additional optimizations to remove the denominators:
```
DRLMatrixNoDenom(F,vars,params,char):
```
* A polynomial in params that encodes the non-specialization locus of the DRL matrix:
```
NonSpecializationPolynomial(F,vars,params,char)
```
* Compute a monomial basis w.r.t. DRL ordering of the quotient ring:
```
QuotientBasis(F,vars,params,char):
```
* Compute at least one point per connected component of the semi-algebraic set defined by `w != 0` where `w` is a polynomial in `params`:
```
SamplePoints(w,params):
```
* Classify the real solutions of a given system of polynomial equations:
```
RealRootClassification(F,vars,params):
```

More details will be provided later.
