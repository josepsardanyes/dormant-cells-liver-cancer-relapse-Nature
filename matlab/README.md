# MATLAB ABC-SMC analysis

`abc_smc_cnio_bigger_range_muQS_muQR.m` contains the ODE model, ABC-SMC inference procedure, and diagnostic plotting code supplied by Tomas Alarcon.

Set `sample_id` near the beginning of the script to `1`, `2`, `3`, or `4`. The program reads the corresponding `volume_treated_X_NoZeros.dat` file and writes sample-specific outputs under `results/`.

The authors confirm successful execution with MATLAB R2024b. Twenty-five ABC-SMC populations take approximately 45–60 minutes. The principal displayed figure is `ABC-SMC Posterior Marginals (final population)`. The program supports Extended Data Fig. 4.

The script calls `randsample`, `ecdf`, `gscatter`, and `boxplot`, which MathWorks documents under the Statistics and Machine Learning Toolbox, and `heaviside`, documented under Symbolic Math Toolbox. These toolboxes are therefore listed as dependencies for the deposited version.
