# ChaoticShearAlfvenWaves (CSAW)

Package for computing the spectrum of Shear Alfven Waves in fusion plasmas with chaotic magentic fields.
This packages solved the generalised eigenvalue problem of the reduced ideal MHD equations for shear Alfven waves.

A description of this package is available at https://arxiv.org/abs/2511.21976

## Installation 

To install, within julia, run
```julia
] add https://github.com/Syndrius/ChaoticShearAlfvenWaves.git
```

## Companion Packages

There are three additional packages to be used with this one.
 - CSAWViz : Visualisation and plotting of solutions.
 - CSAWParallel : Parallel implementation, required for large resolutions.
 - CSAWCantori : Computation of QFM surfaces.

 These can be installed via
```julia
] add https://github.com/Syndrius/CSAWViz.git
] add https://github.com/Syndrius/CSAWParallel.git
] add https://github.com/Syndrius/CSAWCantori.git
```