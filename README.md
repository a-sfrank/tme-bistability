# tme-bistability


Copyright (c) 2025 by the authors. All rights reserved.

## Citation

If you use this code or parts of it, please cite the following references:

(Placeholder for now: 1. **Frank, A. S., Larripa, K., Ryu, H., and Röblitz S.** (2022). Macrophage phenotype transitions in a stochastic gene-regulatory network model. [bioRxiv preprint:10.1101/2022.10.21.513139](https://biorxiv.org/cgi/content/short/2022.10.21.513139v1).)

## See LICENSE


## Questions or support contact:

- **Anna-Simone Frank** (email: asfrank88@gmail.com)


## Download the code

You can download the code at: [https://github.com/a-sfrank/tme-bistability.git](https://github.com/a-sfrank/tme-bistability.git)

## General Information

The code studies bistable dynamics in macrophage–tumor interactions within the tumor microenvironment. It solves a system of ODEs, calculates steady states and their stability. In addition, it determines the basin of attractions  and produces bifurcation diagrams.

The code generates the results presented in the article:

> (Placeholder: **Frank, A. S., Larripa, K., Ryu, H., and Röblitz S.** (2022). Macrophage phenotype transitions in a stochastic gene-regulatory network model. [bioRxiv preprint:10.1101/2022.10.21.513139](https://biorxiv.org/cgi/content/short/2022.10.21.513139v1).)

The code implementation is based on the following articles:

- **Frank, A. S., Larripa, K., Ryu, H., and Röblitz S.** (2022). Macrophage phenotype transitions in a stochastic gene-regulatory network model. [bioRxiv preprint:10.1101/2022.10.21.513139](https://biorxiv.org/cgi/content/short/2022.10.21.513139v1).

- **Frank, A. S., Larripa, K., Ryu, H., Snodgrass, R. G., & Röblitz, S.** (2021). [Bifurcation and sensitivity analysis reveal key drivers]()

## Included Code Files and Descriptions

**Folders**:

1. `Basin_of_Attraction` folder contains following files:

- **parameters.m:** Parameter file with baseline parameter values and specific cases.
- **population_model_v2.m [Main run file]:** Files solves ODE system, determines steady states and calculates stability.It also plots time-series of dynamics.
- **Jacobian_Sym_population_model:** Files symbolically calculates the Jacobian matrix of the ODE system, which is needed for stability analysis.
- **Basin_of_attraction [Main run file]:** File calculates the Basin of attractions and takes as input the steady states determined from population_model_v2.m. You are prompted to specify the parameter case.

2. 'Bifurcation_stability_analysis` folder contains the following files:

- **parameters.m:** Parameter file with baseline parameter values and specific cases.
- **odefun:** Specifies the ODE equations with adaptations to bifrucation analysis.
- **Bifurcation_Analysis_ODE23s.m [Main run file]:** Files calculates the bifurcation plots for different bifurcation parameters, tumor inital conditiosn and parameter cases. This information has to be `specified within the file` before running it.


## How to Run the Code

- Open one of the [Main run files]_ 
  - For **population_model_v2.m** you are promoted to specify which parameter case to run.
  - For **Basin_of_attraction** you are prompted to specify which plot you want, e.g., `M0 vs T` or `M2 vs T`. By default it always runs parameter case: `Case = 2`. 
  - For **Bifurcation_Analysis_ODE23s.m** specify specific parameters wihtin the file before running manually.

For additional details, check the artcile or contact the authors.
