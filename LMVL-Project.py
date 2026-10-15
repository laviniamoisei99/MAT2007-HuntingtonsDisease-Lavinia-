# ============================================================
# MAT2007 - PROJECT
# ============================================================

# Dataset: GSE135589

# ============================================================
# IMPORT LIBRARIES

import os
import gzip
import urllib.request

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# PROJECT FOLDERS

# this keeps the paths working even if someone runs the file from another folder
PROJECT_FOLDER = os.path.dirname(os.path.abspath(__file__))

DATA_FOLDER = os.path.join(PROJECT_FOLDER, "data")

IMAGES_FOLDER = os.path.join(PROJECT_FOLDER, "images")

# ============================================================
# FUNCTION 1: DOWNLOAD DATASET

def download_dataset():

    url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        "GSE135nnn/GSE135589/matrix/"
        "GSE135589_series_matrix.txt.gz"
    )

    if not os.path.exists(DATA_FOLDER):
        os.makedirs(DATA_FOLDER)
        print("\nData folder created.")

    file_path = os.path.join(
        DATA_FOLDER,
        "GSE135589_series_matrix.txt.gz"
    )

    # use the downloaded file if it is already there, so it does not download again
    if os.path.exists(file_path):

        print("\nDataset already available.")
        print("Using:", file_path)

    else:

        print("\nDataset not found on this computer.")
        print("Downloading GSE135589 from NCBI GEO...")
        print("Please wait until the download is complete.")

        urllib.request.urlretrieve(url, file_path)

        print("Download complete.")
        print("Dataset saved as:", file_path)

    return file_path


# ============================================================
# FUNCTION 2: READ SAMPLE METADATA

def read_metadata(file_path):

    metadata = {}

    # the geo file is compressed, so we read it without unzipping it manually
    with gzip.open(file_path, "rt") as file:

        for line in file:

            if line.startswith("!series_matrix_table_begin"):
                break

            if line.startswith("!Sample"):

                parts = line.strip().split("\t")
                category = parts[0]

                values = []

                for value in parts[1:]:
                    values.append(value.replace('"', ''))

                if category not in metadata:
                    metadata[category] = []

                metadata[category].append(values)

    return metadata


# ============================================================
# FUNCTION 3: READ GENE-EXPRESSION DATA

def read_expression_data(file_path):

    table_start = 0
    line_number = 0

    with gzip.open(file_path, "rt") as file:

        for line in file:

            line_number += 1

            if line.startswith("!series_matrix_table_begin"):
                table_start = line_number
                break

    # skip the metadata and load only the table with expression values
    expression_data = pd.read_csv(
        file_path,
        compression="gzip",
        sep="\t",
        skiprows=table_start,
        index_col=0,
        comment="!"
    )

    clean_columns = []

    for column in expression_data.columns:
        clean_columns.append(column.replace('"', ''))

    expression_data.columns = clean_columns

    return expression_data


# ============================================================
# FUNCTION 4: CREATE SAMPLE INFORMATION TABLE

def create_sample_information(metadata):

    # match each sample id with the disease stage written in its metadata
    sample_ids = metadata["!Sample_geo_accession"][0]

    characteristics = metadata["!Sample_characteristics_ch1"]

    groups = []

    for i in range(len(sample_ids)):

        sample_text = ""

        for row in characteristics:

            if i < len(row):
                sample_text += " " + row[i]

        # lowercase helps match the stage names even if their capitalization changes
        sample_text = sample_text.lower()

        group = "Unknown"

        if "disease stage: control" in sample_text:

            group = "Control"

        elif "disease stage: prehd a" in sample_text:

            group = "preHD A"

        elif "disease stage: prehd b" in sample_text:

            group = "preHD B"

        elif (
            "disease stage: zhd stage1" in sample_text
            or "diseasestage: zhd stage1" in sample_text
        ):

            group = "zHD stage1"

        elif (
            "disease stage: zhd stage2" in sample_text
            or "diseasestage: zhd stage2" in sample_text
        ):

            group = "zHD stage2"

        groups.append(group)

    sample_information = pd.DataFrame({
        "Sample": sample_ids,
        "Group": groups
    })

    return sample_information


# ============================================================
# FUNCTION 5: CALCULATE SUMMARY STATISTICS

def summary_statistics(data):

    mean = np.mean(data)

    median = np.median(data)

    standard_deviation = np.std(data)

    return mean, median, standard_deviation


# ============================================================
# FUNCTION 6: CALCULATE UNCERTAINTY

def calculate_uncertainty(data, n_subsamples):

    # Divide observations into five subsamples.
    # split the values into groups to check how much their averages vary

    subsamples = np.array_split(data, n_subsamples)

    subsample_means = []

    for subsample in subsamples:

        subsample_mean = np.mean(subsample)

        subsample_means.append(subsample_mean)

    # Unweighted mean of subsample means.

    mean = np.mean(subsample_means)

    # Standard deviation of subsample means.

    uncertainty = np.std(subsample_means)

    return mean, uncertainty


# ============================================================
# FUNCTION 7: GET SAMPLE IDS FOR EACH GROUP

def get_group_samples(sample_information):

    group_order = [
        "Control",
        "preHD A",
        "preHD B",
        "zHD stage1",
        "zHD stage2"
    ]

    group_samples = {}

    for group in group_order:

        samples = sample_information[
            sample_information["Group"] == group
        ]["Sample"].tolist()

        group_samples[group] = samples

    return group_samples


# ============================================================
# FUNCTION 8: FIND PROBE WITH LARGEST GROUP DIFFERENCE

def find_probe(expression_data, group_samples):

    group_order = [
        "Control",
        "preHD A",
        "preHD B",
        "zHD stage1",
        "zHD stage2"
    ]

    selected_probe = ""

    largest_range = -1

    # check every probe rather than picking a gene in advance
    for probe in expression_data.index:

        group_means = []

        for group in group_order:

            samples = group_samples[group]

            available_samples = []

            for sample in samples:

                if sample in expression_data.columns:
                    available_samples.append(sample)

            values = expression_data.loc[
                probe,
                available_samples
            ].to_numpy(dtype=float)

            group_mean = np.mean(values)

            group_means.append(group_mean)

        # compare the highest and lowest group averages for this probe
        expression_range = (
            max(group_means) - min(group_means)
        )

        if expression_range > largest_range:

            largest_range = expression_range

            selected_probe = probe

    return selected_probe, largest_range


# ============================================================
# FUNCTION 9: GET VALUES FOR SELECTED PROBE

def get_group_values(
    expression_data,
    group_samples,
    selected_probe
):

    group_order = [
        "Control",
        "preHD A",
        "preHD B",
        "zHD stage1",
        "zHD stage2"
    ]

    group_values = {}

    for group in group_order:

        samples = group_samples[group]

        available_samples = []

        for sample in samples:

            if sample in expression_data.columns:
                available_samples.append(sample)

        values = expression_data.loc[
            selected_probe,
            available_samples
        ].to_numpy(dtype=float)

        group_values[group] = values

    return group_values


# ============================================================
# FUNCTION 10: GET GENE NAME

def get_gene_name(selected_probe):

    # Probe 213831_at corresponds to HLA-DQA1.
    # The gene name is added after automatic probe selection.

    # only use the gene name if this exact probe was selected

    if selected_probe == "213831_at":

        gene_name = "HLA-DQA1"

    else:

        gene_name = selected_probe

    return gene_name


# ============================================================
# FUNCTION 11: PRINT SUMMARY STATISTICS

def print_results(group_values):

    group_order = [
        "Control",
        "preHD A",
        "preHD B",
        "zHD stage1",
        "zHD stage2"
    ]

    print("\nSUMMARY STATISTICS")
    print("------------------")

    for group in group_order:

        data = group_values[group]

        mean, median, standard_deviation = (
            summary_statistics(data)
        )

        mean_uncertainty, uncertainty = (
            calculate_uncertainty(data, 5)
        )

        print("\nGroup:", group)

        print("Number of samples:", len(data))

        print("Mean:", round(mean, 3))

        print("Median:", round(median, 3))

        print(
            "Standard deviation:",
            round(standard_deviation, 3)
        )

        print("Uncertainty:", round(uncertainty, 3))


# ============================================================
# FUNCTION 12: PRINT QUANTITATIVE GROUP DIFFERENCES

def print_group_differences(group_values):

    group_order = [
        "preHD A",
        "preHD B",
        "zHD stage1",
        "zHD stage2"
    ]

    # control is the reference for the differences printed below
    control_mean = np.mean(group_values["Control"])

    print("\nDIFFERENCES RELATIVE TO CONTROL")
    print("-------------------------------")

    print("Control mean:", round(control_mean, 3))

    largest_difference = 0

    largest_group = ""

    for group in group_order:

        group_mean = np.mean(group_values[group])

        difference = group_mean - control_mean

        print("\nControl ->", group)

        print("Group mean:", round(group_mean, 3))

        print(
            "Mean difference:",
            round(difference, 3),
            "log2 expression units"
        )

        if abs(difference) > abs(largest_difference):

            largest_difference = difference

            largest_group = group

    print("\nLargest difference from Control:")

    print(
        largest_group,
        "=",
        round(largest_difference, 3),
        "log2 expression units"
    )


# ============================================================
# FUNCTION 13: CREATE PLOT 1 - EXPRESSION DISTRIBUTION

def create_distribution_plot(
    group_values,
    selected_probe,
    gene_name
):

    group_order = [
        "Control",
        "preHD A",
        "preHD B",
        "zHD stage1",
        "zHD stage2"
    ]

    plot_data = []

    for group in group_order:

        plot_data.append(group_values[group])

    plt.figure(figsize=(9, 6))

    # this shows the medians, spread, and outliers for each stage
    plt.boxplot(
        plot_data,
        tick_labels=group_order,
        patch_artist=True,

        boxprops=dict(
            facecolor="white",
            edgecolor="#6F6F7A",
            linewidth=1.3
        ),

        medianprops=dict(
            color="#7563A8",
            linewidth=2.2
        ),

        whiskerprops=dict(
            color="#6F6F7A",
            linewidth=1.2
        ),

        capprops=dict(
            color="#6F6F7A",
            linewidth=1.2
        ),

        flierprops=dict(
            marker="o",
            markerfacecolor="#B7A9D6",
            markeredgecolor="#7563A8",
            markersize=5,
            alpha=0.6
        )
    )

    plt.xlabel("Huntington's disease stage")

    plt.ylabel("Log2 RMA expression intensity")

    plt.title(
        gene_name
        + " expression across Huntington's disease stages\n"
        + "Probe "
        + str(selected_probe)
    )

    plt.grid(axis="y", alpha=0.2)

    plt.tight_layout()

    plot_path = os.path.join(
        IMAGES_FOLDER,
        "HD_expression_distribution.png"
    )

    plt.savefig(plot_path, dpi=300)

    plt.show()

    plt.close()


# ============================================================
# FUNCTION 14: CREATE PLOT 2
# MEAN EXPRESSION WITH UNCERTAINTY
# ============================================================

def create_uncertainty_plot(
    group_values,
    selected_probe,
    gene_name
):

    group_order = [
        "Control",
        "preHD A",
        "preHD B",
        "zHD stage1",
        "zHD stage2"
    ]

    means = []

    uncertainties = []

    for group in group_order:

        data = group_values[group]

        mean, uncertainty = calculate_uncertainty(data, 5)

        means.append(mean)

        uncertainties.append(uncertainty)

    x = [0, 1, 2, 3, 4]

    plt.figure(figsize=(9, 6))

    # the bars show variation between subsample means, not confidence intervals
    plt.errorbar(
        x,
        means,
        yerr=uncertainties,
        fmt="o-",
        color="#7563A8",
        markerfacecolor="#A9C5E8",
        markeredgecolor="#7563A8",
        markersize=8,
        linewidth=1.8,
        capsize=5
    )

    for i in range(len(x)):

        label = (
            str(round(means[i], 2))
            + " ± "
            + str(round(uncertainties[i], 2))
        )

        plt.text(
            x[i],
            means[i] + uncertainties[i] + 0.08,
            label,
            ha="center",
            fontsize=9
        )

    plt.xticks(x, group_order)

    plt.xlim(-0.35, 4.35)

    plt.xlabel("Huntington's disease stage")

    plt.ylabel("Mean log2 RMA expression intensity")

    plt.title(
        "Mean "
        + gene_name
        + " expression with uncertainty\n"
        + "Probe "
        + str(selected_probe)
    )

    plt.grid(axis="y", alpha=0.2)

    plt.tight_layout()

    plot_path = os.path.join(
        IMAGES_FOLDER,
        "HD_expression_uncertainty.png"
    )

    plt.savefig(plot_path, dpi=300)

    plt.show()

    plt.close()


# ============================================================
# MAIN FUNCTION

def main():

    # STEP 1: START PROGRAM

    print("\nHUNTINGTON'S DISEASE GENE EXPRESSION ANALYSIS")
    print("=============================================")

    print("\nProject folder:")
    print(PROJECT_FOLDER)

    # STEP 2: CREATE IMAGES FOLDER

    if not os.path.exists(IMAGES_FOLDER):

        os.makedirs(IMAGES_FOLDER)

        print("\nImages folder created.")

    else:

        print("\nImages folder already available.")

    print("Images will be saved in:")
    print(IMAGES_FOLDER)

    # STEP 3: DOWNLOAD DATASET

    file_path = download_dataset()

    # STEP 4: READ METADATA

    print("\nReading sample information...")

    metadata = read_metadata(file_path)

    sample_information = create_sample_information(metadata)

    # STEP 5: READ EXPRESSION DATA

    print("Reading gene-expression data...")

    expression_data = read_expression_data(file_path)

    print("Data loaded successfully.")

    # STEP 6: DISPLAY DATASET INFORMATION

    print("\nDATASET INFORMATION")
    print("-------------------")

    print("Number of samples:", len(expression_data.columns))

    print("Number of probe sets:", len(expression_data.index))

    print("\nFirst five rows of the expression data:")

    print(expression_data.head())

    # STEP 7: DISPLAY AND CHECK SAMPLE GROUPS

    print("\nSAMPLE GROUPS")
    print("-------------")

    group_counts = sample_information["Group"].value_counts()

    print(group_counts)

    # stop if any sample could not be assigned to a stage
    unknown_samples = sample_information[
        sample_information["Group"] == "Unknown"
    ]

    if len(unknown_samples) > 0:

        print("\nERROR:")
        print("Some samples could not be assigned to a disease stage.")
        print(unknown_samples)

        return

    # STEP 8: CHECK TOTAL SAMPLE NUMBER

    if len(sample_information) != 178:

        print("\nWARNING:")
        print(
            "Expected 178 samples but found:",
            len(sample_information)
        )

    else:

        print("\nAll 178 samples were successfully classified.")

    # STEP 9: ORGANISE SAMPLE IDS BY GROUP

    group_samples = get_group_samples(sample_information)

    # STEP 10: SELECT PROBE

    print(
        "\nSearching the probe sets for the largest "
        "difference across disease stages..."
    )

    print("This may take a moment.")

    selected_probe, expression_range = find_probe(
        expression_data,
        group_samples
    )

    # STEP 11: IDENTIFY SELECTED GENE

    gene_name = get_gene_name(selected_probe)

    print("\nSELECTED PROBE")
    print("--------------")

    print("Probe:", selected_probe)

    print("Gene:", gene_name)

    print(
        "Range between highest and lowest group mean:",
        round(expression_range, 3),
        "log2 expression units"
    )

    print("\nThis is an exploratory selection.")

    print(
        "The probe was selected because it had the largest "
        "range between group means in this dataset."
    )

    # STEP 12: EXTRACT VALUES FOR SELECTED PROBE

    group_values = get_group_values(
        expression_data,
        group_samples,
        selected_probe
    )

    # STEP 13: PRINT SUMMARY STATISTICS

    print_results(group_values)

    # STEP 14: PRINT QUANTITATIVE DIFFERENCES

    print_group_differences(group_values)

    # STEP 15: CREATE PLOT 1

    print("\nCreating Plot 1:")
    print("Expression distribution across disease stages.")

    create_distribution_plot(
        group_values,
        selected_probe,
        gene_name
    )

    print(
        "Plot 1 saved as:",
        os.path.join(
            IMAGES_FOLDER,
            "HD_expression_distribution.png"
        )
    )

    # STEP 16: CREATE PLOT 2

    print("\nCreating Plot 2:")
    print("Mean expression with uncertainty.")

    create_uncertainty_plot(
        group_values,
        selected_probe,
        gene_name
    )

    print(
        "Plot 2 saved as:",
        os.path.join(
            IMAGES_FOLDER,
            "HD_expression_uncertainty.png"
        )
    )

    # STEP 17: FINISH

    print("\nANALYSIS COMPLETE")
    print("=================")

    print("\nExactly two plots were created.")

    print("Both plots are located in:")
    print(IMAGES_FOLDER)


# ============================================================
# RUN THE PROGRAM

main()