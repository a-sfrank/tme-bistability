[![DOI](https://zenodo.org/badge/988167516.svg)](https://doi.org/10.5281/zenodo.15772301)

# tme-bistability

Copyright (c) 2025 by the authors. All rights reserved.

### Citation

If you use this code or parts of it, please cite the following references:

**Ryu, H., Röblitz S., Larripa, K. and Frank, A. S.** (2025). Modeling Bistable Dynamics Arising from Macrophage–Tumor Interactions in the Tumor Microenvironment. bioRXiv preprint. doi: [https://doi.org/10.1101/2025.06.24.661275](https://doi.org/10.1101/2025.06.24.661275).

#### See LICENSE for more details

### Questions or support contact:

- **Anna-Simone Frank** (email: asfrank88@gmail.com)

### Download the code

You can download the code at: [https://github.com/a-sfrank/tme-bistability.git](https://github.com/a-sfrank/tme-bistability.git)

### General Information

The code studies bistable dynamics in macrophage–tumor interactions within the tumor microenvironment. It solves a system of ODEs, calculates steady states and their stability. In addition, it determines the basin of attractions and produces bifurcation diagrams.

The code generates the results presented in the article:

**Ryu, H., Röblitz S., Larripa, K. and Frank, A. S.** (2025). Modeling Bistable Dynamics Arising from Macrophage–Tumor Interactions in the Tumor Microenvironment. bioRXiv preprint. doi: [https://doi.org/10.1101/2025.06.24.661275](https://doi.org/10.1101/2025.06.24.661275).

### Included Code Files and Descriptions

**Folders**:

1. `Basin_of_Attraction` folder contains following files:

- **parameters.m:** Parameter file with baseline parameter values and specific cases.
- **population_model_v2.m [Main run file]:** Files solves ODE system, determines steady states and calculates stability. It also plots time-series of dynamics.
- **Jacobian_Sym_population_model:** Files symbolically calculates the Jacobian matrix of the ODE system, which is needed for stability analysis.
- **Basin_of_attraction [Main run file]:** File calculates the Basin of attractions and takes as input the steady states determined from population_model_v2.m. You are prompted to specify the parameter case.

2. `Bifurcation_stability_analysis` folder contains the following files:

- **parameters.m:** Parameter file with baseline parameter values and specific cases.
- **odefun:** Specifies the ODE equations with adaptations to bifrucation analysis.
- **Bifurcation_Analysis_ODE23s.m [Main run file]:** Files calculates the bifurcation plots for different bifurcation parameters, tumor inital conditions and parameter cases. This information has to be `specified within the file` before running it.

### Output:

- Results are mainly figures which are saved to Figures/-folders.

### How to Run the Code

- Open one of the [Main run files]:
  - For **population_model_v2.m** you are promoted to specify which parameter case to run.
  - For **Basin_of_attraction** you are prompted to specify which plot you want, e.g., `M0 vs T` or `M2 vs T`. By default it always runs parameter case: `Case = 2`. 
  - For **Bifurcation_Analysis_ODE23s.m** specify specific parameters within the file before running manually.

For additional details, check the artcile or contact the authors.
