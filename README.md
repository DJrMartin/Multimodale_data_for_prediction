# Multi-omics benchmark spectral datasets in a murine model of liver

This repository provides a comprehensive collection of multi-omics datasets---including metallomics, Raman spectroscopy, and mid-infrared (MIR) spectroscopy---acquired from a murine experiment, involving 100 mice divided into four experimental groups. This study aim to induce hepatic steatosis, fibrosis, and inflammation. These pathological states were then assessed through various omics-based analytical techniques.

## Installation

git clone <https://github.com/ton-nom-utilisateur/NomDuProjet.git>

Developed in **R**. Required packages:

``` r
install.packages(c("stringr" , "FactoMineR" , "prospectr", "fda", "randomForest", "caret"))
```

## Experimental design

A total of 100 mice were distributed across four experimental groups:

-   **Control (CTL)** -- standard diet
-   **High-Fat Diet (HF)** -- to induce metabolic stress
-   **Iron dextran (Iron)** -- to induce hepatic iron overload
-   **Iron + High-Fat diet (Iron-HF)** -- combined model of metabolic and iron-induced liver injury

The primary objective was to develop a robust preclinical model of liver dysfunction and to monitor dysfunction progression over time using complementary -omics approaches.

## Data modalities and time points

The dataset includes the following modalities:

-   **Mid-Infrared Spectroscopy (MIR)**: Molecular fingerprinting based on infrared absorption (collected at **2 months** and **5 months**).

-   **Raman spectroscopy**: Vibrational signatures of biochemical compounds (collected at **5 months** only).

-   **Metallomics**: Elemental profiling to assess iron and other trace metals (collected at **5 months** only).

This design allows for longitudinal assessment of MIR profiles and endpoint characterization through Raman and metallomic analyses. This dataset serves to: (1) Characterize liver disease progression using multi-omics techniques. (2) Enable integrative analysis of spectral and elemental data across experimental groups. (3) Support biomarker discovery and non-invasive diagnostic approaches for liver pathology.

## Repository structure

-   `data/`: Raw and/or preprocessed datasets organized by -omics type.
-   `scr/`: Preprocessing or analysis scripts (e.g., normalization, visualization).
    -   `synthetic_parameters_of_RF.R` sums up the liver dysfunction in mice into two quantitative variable: liver steatosis and liver inflammation. Synthetic variables are recorded in `data/data_preprocess/target_variables.rda`.
    -   `mean_of_represnetative_spectra.R` highlights mean spectral patterns for each samples.
    -   `tuning_parameters_of_RF.R` aims at defining the number of trees (`ntree`), the number of variables tried at each split (`mtry`), and the maximum of the nodes (`maxnodes`) in each decision treesof the randomForest.
    -   `performance_prediction_steatosis.R` computing the predictive performance of divers models.
-   `figures/`: Output plots from initial data exploration.
-   `tests/`: Tests to validate functions.
-   `functions`: custom functions for the project.
    -   `representative.mean()` calculates the mean of the spectra that are most representative of a given sample. It works by identifying a group of similar spectra that are most commonly observed within the sample, and then computing their average.
    -   `read.csv.transform()` is a function to import data and specify the transformation. It provides a computing of the second derivates, $\beta$-splines projection and Fourier transformaiton for spectral data, as well as Z-normalisation for other types of data.
    -   `boosting_rf_cv()` function implements a custom boosting model based on **Random Forests**, combined with **K-fold cross-validation**. The core idea is to sequentially model the residuals of a primary Random Forest trained on one dataset using additional datasets, mimicking a boosting-like behavior across different data sources.
-   `test_functions.R`: example usage and a validation test of the custom functions.

## License

This repository is made available under the [MIT License](LICENSE).

## Citation

If you use this resource in your research, please cite it as:

> [Martin D. et al.] *Multi-omics spectral datasets in a murine model of liver dysfunctions*. GitHub Repository, [2025].

## Contact

For questions, feedback, or collaboration opportunities, please contact:\
**David Martin (PhD)** -- *University of Rennes*\
[[david.martin.2\@univ-rennes.fr](mailto:david.martin.2@univ-rennes.fr){.email}]
