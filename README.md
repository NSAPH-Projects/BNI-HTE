# Environmental Justice Implications of Power Plant Emissions Control Policies: Heterogeneous Causal Effect Estimation under Bipartite Network Interference

### Kevin Chen, Falco Bargagli-Stoffi, Raphael Kim, Rachel Nethery

Updated 8/9/2024

Repository for simulation and analysis code used in the paper entitled "Environmental Justice Implications of Power Plant Emissions Control Policies: Heterogeneous Causal Effect Estimation under Bipartite Network Interference", authored by Kevin L. Chen, Falco J. Bargagli Stoffi, Raphael C. Kim, and Rachel C. Nethery. The preprint of this paper can be found on arXiv at https://arxiv.org/abs/2304.12500.

**Overview:** Estimation of heterogeneous treatment effects in the setting of bipartite network interference using G-computation and augmented inverse propensity weighting methods. These methods are empirically evaluated using a proposed empirical Monte Carlo simulation approach and applied to measure the causal effects of flue gas desulfurization equipment ("scrubber") installation on coal-fired power plants on ischemic heart disease hospitalizations among the Medicare population in the United States.

***

[`sim_data_generation.R`](https://github.com/NSAPH-Projects/emissions-ihd-bipartite/blob/master/code/sim_data_generation.R) contains the code used to generate semi-synthetic datasets for simulations according to the simulation approach diagram:

![sim_approach_diagram(1) drawio](https://user-images.githubusercontent.com/42856787/229964368-760977d1-0624-4845-9b08-2dbb60e0150e.png)

Datasets involve real-world interference structures and covariates. Treatments are generated synthetically and propensity score estimates are computed according to a number of different misspecification scenarios.

[`simulations.R`](https://github.com/NSAPH-Projects/emissions-ihd-bipartite/blob/master/code/simulations.R) contains the code used to run simulations according to various heterogeneity, misspecification, effect size, and outcome model error variance scenarios, as well as to plot and evaluate results.

[`application_analysis.R`](https://github.com/NSAPH-Projects/emissions-ihd-bipartite/blob/master/code/application_analysis.R) contains the code used to process data and conduct the analysis of scrubber installation effects on ischemic heart disease hospitalization rates using proposed G-computation and augmented inverse propensity weighting methods. File also includes bootstrap analysis code for evaluation of population and subgroup-specific treatment effects.

***



