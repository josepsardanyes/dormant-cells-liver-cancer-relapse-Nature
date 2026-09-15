# Data

`demo_regression_data.csv` is a simulated dataset created only to demonstrate that the regression functions run. It is not part of the study and must not be used to reproduce or interpret manuscript findings.

`tumor_growth_data.xlsx` contains the experimental data used for Extended Data Fig. 5a,b. Sheet `Sheet1` is arranged as follows:

- Column B: time after L/G treatment initiation, in days.
- Columns C–Z: eight untreated tumours, each represented by short diameter, long diameter, and volume (`mm^3`).
- Columns AA–BJ: twelve L/G-treated tumours, each represented by short diameter, long diameter, and volume (`mm^3`).

The Python analysis uses the volume column from each three-column tumour block. Blank entries and non-positive volumes are excluded.

`volume_treated_1_NoZeros.dat` through `volume_treated_4_NoZeros.dat` are the datasets supplied with the MATLAB ABC-SMC analysis. The deposited MATLAB script is configured for sample 3. In each file, the first column is used as observation time and the fourth column as tumour volume. The meanings and units of all five columns must be documented by the authors before submission.

Before public release, confirm the provenance and open-sharing status of every experimental dataset. Do not commit restricted or identifiable data to a public repository.
