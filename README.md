**Huntington's Disease Gene Expression Analysis**

**MAT2007 – Introduction to Programming**  
Maastricht Science Programme

**Project**

For this project, I investigated how peripheral blood gene expression varies across different stages of Huntington's disease (HD) with Python (in Visual Studio Code).

The dataset contains expression for 54,675 probe sets across five disease-stage groups. Instead of choosing a specific gene beforehand, I compared the mean expression values of all probes and selected the one with the largest difference between groups. The selected probe was 213831_at, corresponding to *HLA-DQA1*. Then I compared its expression across the five disease stages using summary statistics and two plots.

Research question: How does peripheral blood *HLA-DQA1* expression vary across stages of Huntington's disease?

**Dataset**
The dataset comes from the NCBI Gene Expression Omnibus (GEO).

- **Dataset:** GSE135589
- **Link:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4020348
- **Original study:** Mitchell et al. (2020), *Longitudinal expression changes are weak correlates of disease progression in Huntington's disease*.
- **Publication:** https://doi.org/10.1093/braincomms/fcaa172

The dataset contains 178 peripheral whole-blood samples divided into five groups: Control, preHD A, preHD B, zHD stage 1, and zHD stage 2. Gene expression was measured using the Affymetrix Human Genome U133 Plus 2.0 Array. I used the processed GEO Series Matrix, which already contains log2-transformed and RMA-normalized expression values.

This dataset does not need to be downloaded manually, the script downloads it automatically.

**Requirements**

To run the code, you need Python, NumPy, pandas, and Matplotlib.

**How to run the code**

1. Download the Python script (`LMVL-Project.py`) from this repository.
2. Install the required packages (if you don't have them) by entering the following command in your terminal:

   `python3 -m pip install numpy pandas matplotlib`

3. Open the script in VS Code and select **Run Python File**, or run it from the terminal using:

   `python3 LMVL-Project.py`
   
**Data processing**

No manual cleaning of the dataset is required.

The script performs a few processing steps before the analysis:
- Reads the sample metadata and expression measurements from the compressed GEO file.
- Removes quotation marks from sample identifiers.
- Matches samples to the correct disease-stage groups using the metadata and checks whether any samples could not be classified.
- Calculates the mean expression for each probe within the five groups and selects the probe with the largest range.

For the selected probe, the code calculates means, medians, standard deviations, and differences relative to the Control group. Each disease-stage group is also divided into five subsamples. The standard deviation of their means is used to describe variability between subsamples.

**Outputs**

After running the script, the results are printed in the terminal, including the selected probe, sample counts, summary statistics, and expression differences relative to Control. Two plots are also generated and saved automatically in the `images` folder: `HD_expression_distribution.png` and `HD_expression_uncertainty.png`.

**Report**

The accompanying report explains the dataset, analysis, results, and interpretation in more detail.

This project is an exploratory analysis. The selected probe showed the largest observed difference between group means, but this does not necessarily mean it is the most biologically important gene or a reliable marker of Huntington's disease progression.gene or a reliable marker of Huntington's disease progression.
