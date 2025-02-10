# Repository Overview

Welcome to this GitHub repository, which houses the code corresponding to the article titled "A Flux-Limited Model for Glioma Patterning with Hypoxia-Induced Angiogenesis", DOI: https://www.mdpi.com/2073-8994/12/11/1870. The results are used in Chapter 3 of the thesis: https://kluedo.ub.rptu.de/frontdoor/index/index/docId/6573.

## Key Features

The codebase encompasses a variety of functionalities, including:

- Full model simulations
- Comparison between different models
- Pattern formation analysis

## Project Structure

### Scripts

The main script `main_all.m` runs all the below-mentioned codes, whose results are used in Chapter 3 of the thesis ("Diss_Kumar_Pawan.pdf" present in the parent directory).

### Full_model_I_II (figure 3.2 and 3.3)

This directory contains the codes and results for experiments 1 and 2 mentioned in the 3rd Chapter of the thesis.

### Difference_minimal (Figure 3.4)

The codes & results in this directory are for the difference between the full model and the old model (e.g., glioma density full model - glioma density minimal/old model) in the framework of experiment 2.

### Difference_self_diffusion (figure 3.5)

In this directory, the codes and simulation represent the comparison (by difference) between the full model and the model without limited flux for glioma self-diffusion (gamma_2 = 0, the denominator of ph_taxis = 1) in the framework of experiment 2.

### Pattern (figure 3.6: first and second row)

This folder contains the 1D codes and videos for the full system in dimensional form, also the results for pattern formations.

### Pattern_difference (figure 3.7)

The pattern difference between the old and new model (experiment 3 in 1D).

### Pattern_diff_W_L_flux (figure 3.6: 3rd to 6th row)

The pattern formed by the model without self-diffusion and the difference between the full model and the model without self-diffusion.

Feel free to explore the respective directories to delve into specific aspects of the project. If you have any questions or need further clarification, please don't hesitate to reach out.


Each subdirectory of this repository contains a separate readme file.

### Author

Pawan Kumar: [@its-Pa1](https://github.com/its-Pa1)

© 2025, Pawan Kumar. All Rights Reserved.
