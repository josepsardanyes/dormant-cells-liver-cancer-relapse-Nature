# Tracking dormant liver cancer cells reveals origins of local and metastatic relapses with distinct therapeutic avenues

This repository contains code supporting the computational analyses reported in the manuscript **“Tracking dormant liver cancer cells reveals origins of local and metastatic relapses with distinct therapeutic avenues.”**

## Code contributors

- **Josep Sardanyes** — Python growth-rate estimation code
- **Tomas Alarcon** — MATLAB ABC-SMC code

Both contributors are affiliated with the Centre de Recerca Matemàtica, Edifici C, Campus de Bellaterra, 08193 Cerdanyola del Valles, Barcelona, Spain.

## Repository status

The Python component estimates the untreated and relapse-phase growth rates, `gS` and `gR`. The MATLAB component implements the ABC-SMC parameter-estimation analysis. The manuscript mapping and author validation status are documented below.

## Contents

- `python/fit_gS_untreated.py`: estimates the exponential growth rate, `gS`, from untreated tumour-volume measurements by ordinary least-squares regression of `log(volume)` against time.
- `python/fit_gR_relapse.py`: estimates the relapse-phase exponential growth rate, `gR`. It identifies the minimum measured volume on or after a user-specified starting day and regresses `log(volume)` against time from that observation onward.
- `python/run_demo.py`: runs both regression functions on the included simulated demonstration data.
- `python/generate_extended_data_fig5a.py`: reads the experimental workbook, estimates `gS` for eight untreated tumours, and generates Extended Data Fig. 5a outputs.
- `python/generate_extended_data_fig5b.py`: identifies relapse phases for the nine tumours listed in Supplementary Table 8, estimates `gR`, and generates Extended Data Fig. 5b outputs.
- `data/tumor_growth_data.xlsx`: experimental tumour measurements used for the `gS` and `gR` fits.
- `data/demo_regression_data.csv`: small simulated dataset used solely to test the regression code. It is not manuscript data.
- `matlab/abc_smc_cnio_bigger_range_muQS_muQR.m`: ABC-SMC parameter inference for the five-state ODE model.
- `data/volume_treated_[1-4]_NoZeros.dat`: treated tumour-volume datasets supplied with the MATLAB analysis.
- `docs/NATURE_CHECKLIST.md`: current status against Nature's Code and Software Submission Checklist.

## Correspondence with the manuscript

| Manuscript item | Code |
|---|---|
| Extended Data Fig. 5a | `python/generate_extended_data_fig5a.py` using `fit_gS_untreated.py` |
| Extended Data Fig. 5b | `python/generate_extended_data_fig5b.py` using `fit_gR_relapse.py` |
| Extended Data Fig. 4 | `matlab/abc_smc_cnio_bigger_range_muQS_muQR.m` |

The Python programs reproduce the growth-rate regressions and generate panels for Extended Data Fig. 5a,b. ChatGPT assisted with development of the regression and visualization code; the results were reviewed and validated by the authors. The MATLAB program performs the ABC-SMC fitting underlying Extended Data Fig. 4.

## System requirements

### Software

- Python 3.10 or later
- NumPy 1.26 or later
- SciPy 1.11 or later
- Matplotlib 3.8 or later
- openpyxl 3.1 or later
- MATLAB R2024b
- Statistics and Machine Learning Toolbox (used functions include `randsample`, `ecdf`, `gscatter`, and `boxplot`)
- Symbolic Math Toolbox (the script calls `heaviside`)

The code does not require non-standard hardware and runs on a normal desktop or laptop computer. The authors confirmed successful execution of the MATLAB analysis in MATLAB R2024b. Although no specific toolbox requirement was initially reported, the current script calls functions documented by MathWorks under the Statistics and Machine Learning Toolbox and Symbolic Math Toolbox; these are therefore listed as dependencies.

The Python workflows were tested successfully on Linux x86_64 with Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, and openpyxl 3.1.5.

## Installation

From the repository root, create and activate a virtual environment and install the dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate the environment using `.venv\\Scripts\\activate`.

Typical installation time is expected to be less than five minutes on a normal desktop computer with an existing Python installation and internet connection.

## Demonstration

Run:

```bash
python python/run_demo.py
```

Expected output (allowing for insignificant formatting differences):

```text
Untreated fit: gS = 0.100000 day^-1, R^2 = 1.000000
Relapse fit:   gR = 0.179995 day^-1, R^2 = 1.000000
```

The demo should finish in less than one minute. It verifies that the regression functions execute correctly; it does not reproduce the manuscript results.

## Extended Data Figure 5

From the repository root, run:

```bash
python python/generate_extended_data_fig5a.py
python python/generate_extended_data_fig5b.py
```

The first command analyses all eight untreated tumours and prints the mean and sample standard deviation of `gS`. It creates `Extended_Data_Fig_5a.pdf`, `Extended_Data_Fig_5a.png`, and `gS_estimates.csv` under `results/`.

The second command analyses treated tumours 1–6, 9, 10, and 12, matching Supplementary Table 8. For each tumour, the relapse phase begins at the minimum positive tumour volume observed from day 26 onward and continues to its final positive observation. It creates `Extended_Data_Fig_5b.pdf`, `Extended_Data_Fig_5b.png`, and `gR_estimates.csv`.

Across all nine relapsing tumours, the code gives `gR = 0.180 ± 0.076 day^-1` (mean ± sample SD). The manuscript value `gR = 0.176 ± 0.078 day^-1` is calculated using the six tumours with more than three relapse-phase observations (tumours 1–4, 9, and 10), as specified in Supplementary Table 8.

## MATLAB ABC-SMC analysis

Set `sample_id` near the beginning of the MATLAB script to `1`, `2`, `3`, or `4` to select the corresponding experimental dataset. From MATLAB, change to the repository root and run:

```matlab
run('matlab/abc_smc_cnio_bigger_range_muQS_muQR.m')
```

The first column of the selected `.dat` file supplies observation times and the fourth column supplies tumour volumes. The script fixes the random-number seed at 42, evolves a five-state ODE model, performs ABC-SMC with 300 particles and 25 populations (generations), displays diagnostic figures, and writes numerical results to `results/`.

The authors confirm that the program completes without errors in MATLAB R2024b. A run of 25 populations takes approximately 45–60 minutes. The principal displayed figure is **“ABC-SMC Posterior Marginals (final population)”**. For sample `X`, where `X = 1, 2, 3, 4`, the principal numerical output files are:

- `quiescentcells_sampleX.dat`
- `distances_sampleX.dat`
- `RootMeanSquaredError_sampleX.dat`
- `final_particles_sampleX.dat`
- `BestParsSampleX.dat`
- `MinQuiescSampleX.dat`

The deposited version replaces author-specific absolute paths with repository-relative paths and includes a `sample_id` selector.

## Reproducing the manuscript results

The correspondence between programs and Extended Data Figs. 4 and 5a,b is documented above. The Python analyses read the included `data/tumor_growth_data.xlsx` workbook. The MATLAB analysis reads the included experimental sample files and produces the outputs listed above.

## Data availability and confidentiality

The included CSV file is simulated demonstration data. `tumor_growth_data.xlsx` contains the experimental measurements used for the `gS` and `gR` fits. Four treated tumour-volume `.dat` files are included for the MATLAB analysis. Before public release, the authors must confirm that all experimental datasets are permitted for open distribution and contain no identifiable, confidential, or access-controlled information.

## AI-assisted code development

Draft disclosure for author and journal review:

> OpenAI's ChatGPT was used to assist in writing Python code for standard regression analyses and data visualization. The analytical methods, generated code, and results were reviewed and validated by the authors, who take full responsibility for the analyses and reported results.

## Licence

The code is released under the MIT License. See `LICENSE.txt`.

## Citation

Citation information will be added after the manuscript bibliographic details and repository authors have been confirmed.
