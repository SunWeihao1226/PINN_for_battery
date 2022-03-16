# PINN_for_battery

## Introduction

This repository contains the code and results for our REX research project. The code builds a Physics Informed Neural Network (PINN) for a Single Particle Model (SPM) for Lithium-Ion batteries. The goal of the model is to solve the partial differential equations proposed by Marquis et al. [[4]](#3) to improve the current battery models. 

The built-in PINN in the NeuralPDE package of Julia [[1]](#1) is used in this project. The `pinn_spm.jl` models our PINN to solve the positive and negative electrode concentrations (i.e., $C_{s,p}$ and $C_{s,n}$) with respect to scaled time `t` and particle radius `r`. The Neural Network contains three layers including the input and output, and the hidden layer consists of 15 neurons in this model. 10,000 iterations were performed to solve the equations. Both time and radius are discretized into 11 points in $[0,1]$.  The `ode_spm` provides the numerical discretization of the same equations using the MethodOfLines package. [[2]](#2) These results are useful to calculate errors. `analysis.jl` includes the pipeline of running the PINN model and MethodOfLines discretization, as well as the errors and plots for each result. The error is calculated as the absolute value of the difference of two results (i.e., $|C_{s,r\space MethodOfLine} - C_{s,r\space PINN}, r\in  \{n,p\}|$). The analytical plots are generated and saved in the `plots` folder. The gif plots are the animation of concentrations against $dt$ in changes of $dr$. The `c_sn_error` and `c_sp_error` are heat maps of errors. When $dr$ stays the same, the error is larger when $dt$ is larger.

The parameters in the models related to the network structure are changeable to make sure the code is robust, which makes the model is easy to maintain in future improvement and usage of solving similar equation problems. Detailed explanation is provided in the Usage part.


## Usage
1. Make sure that the Julia environment (at least 1.7.1) is installed. If not, the packages can be downloaded [here](https://julialang.org/downloads/).
2. Clone this repo using ```
git clone https://github.com/luckyberen/PINN_for_battery.git```
3. If you are using VS Code, install the [Julia Extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia)
4. You can change the maximum iterations, number of hidden layer neurons, and discretization levels in `pinn_spm.jl` and `ode_spm.jl` if you want to explore the model more. Make sure that the `dt` and `dr` in the two files are the same.
5. Start the Julia environment. Enter the repository root, use `cd("src")` in Julia terminal to enter the code path.
6. Run it by `include("analysis.jl")`. The analysis results will automatically generated.

## Contribute
The other contributors of this Rex project: Maricela, Chloe, Harry, Emerald, Roy.




## reference

<a id=1>[1]</a> The PINN code structures are based on the Julia tutorial:
https://neuralpde.sciml.ai/stable/pinn/system/

<a id=1>[2]</a>MethodOfLines discretization: 
https://github.com/SciML/MethodOfLines.jl

<a id=2>[3]</a>Cen, Z, Kubiak, P. Lithium-ion battery SOC/SOH adaptive estimation via simplified single particle model. Int J Energy Res. 2020; 44: 12444– 12459. https://doi.org/10.1002/er.5374

<a id=3>[4]</a>Marquis, S. G., Sulzer, V., Timms, R., Please, C. P., & Chapman, S. J. (2019). An asymptotic derivation of a single particle model with electrolyte. arXiv [physics.chem-ph]. Opgehaal van http://arxiv.org/abs/1905.12553

