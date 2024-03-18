SASC code
================

This GitHub repository contains code to reproduce simulation results, figures, and real data analysis results from the paper "Immune profiling among colorectal cancer subtypes using dependent mixture models".

## File directories

### code folder

#### R functions

SASC_func.R contains the main functions of the proposed SASC method.

realdata_func.R contains functions for pre- and post-processing of real data analysis and visualization.

comp_func.R contains contains other functions for simulation visualization and comparison with benchmark methods.

#### Real data analysis

SASC_CRC.Rmd reproduces real data analysis results. (Figure 1, 5, S7, and S8) See [here]() for the corresponding html output file as a tutorial. 

#### Simulation studies

simulation_sc1.Rmd reproduces simulation study scenario 1 results. (Figure 2, 3, 4 and Table 2) See [here]() for the corresponding html file as a tutorial. 

Similarly, see simulation_sc2.Rmd [here](), simulation_sc3.Rmd [here]() for scenario 2 and 3. (Figures in Appendix D)

### data folder

The data folder contains real data and data generated in simulation studies.

simulate_data_sc1.Rmd [here]() reproduces the data generation process in simulation study scenario 1. Codes are also provided for scenario 2 and 3 in simulate_data_sc2.Rmd and simulate_data_sc3.Rmd.

### output folder

The output folder contains all output data, figures, and tables. 











