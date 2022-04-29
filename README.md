# BattPhase.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bradyplanden.github.io/LiMetalPhaseFields.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bradyplanden.github.io/LiMetalPhaseFields.jl/dev)
[![Build Status](https://github.com/bradyplanden/LiMetalPhaseFields.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bradyplanden/LiMetalPhaseFields.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/bradyplanden/LiMetalPhaseFields.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/bradyplanden/LiMetalPhaseFields.jl)

BattPhase.jl provides a Julia framework for solving lithium-metal 2D phase field problems. This work is presented to support pre-print: [10.1149/osf.io/k2vu6]

&nbsp;

Install (Julia 1.7 and later)
-----------------------------

```julia
(v1.7) pkg> add https://github.com/BradyPlanden/BattPhase.jl
```

(Type `]` to enter package mode.)

<!-- &nbsp;
## Examples 
Run the semi-circle example via,
```julia
include("examples/Semi-example.jl")
```
-->

&nbsp;
## Results
<p align="center">
<img src="examples/semicircle_fps15.gif" width="900" align="center"  />
</p>

[10.1149/osf.io/k2vu6]: https://ecsarxiv.org/k2vu6/


&nbsp;
## Code Timing - Intel 10980XE - Upwind Scheme

<div align="center">
  
|Number of Node Points|Runge-Kutta 3 |Runge-Kutta 3 Approximation|
|:-:|:-:|:-:|
| 10<sup>2</sup>  |  0.0037 | 0.0013  |
|  20<sup>2</sup> | 0.0279  | 0.0091  |
|  40<sup>2</sup>|  0.183 | 0.0614 |
|  80<sup>2</sup>| 1.42  | 0.487 |
|  160<sup>2</sup> |  15.0 | 5.19  |
|  320<sup>2</sup> |  170.6 | 57.5 |
  
</div>
