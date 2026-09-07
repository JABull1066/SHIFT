# SHIFT

SHIFT (Spatial Hierarchy of Interactions and Feature Transitions) is a Python package accompanying:

> Bull et al. (2026), *Integrating spatial descriptors quantifies changing cell interactions as colorectal neoplasia progresses*.

The package provides tools to:

1. Generate a consistent set of spatial summary statistics from annotated spatial data.
2. Calculate SHIFT scores describing changes in cell-cell interactions between biological states.

While SHIFT is compatible with statistics calculated from any quantitative spatial biology package, the implementation presented here is designed to facilitate analyses using our [MuSpAn package](www.muspan.co.uk) - please see the corresponding paper here: [https://doi.org/10.1038/s41467-026-75649-7](https://doi.org/10.1038/s41467-026-75649-7). This implementation (i.e., ``shift.compute_shift(...)``) can be used alongside dataframes generated using any generic pipeline, but care must be taken that the naming conventions match those produced by ``shift.generate_feature_dataframe(domain, ...)``. The aim of this repository is to provide a simple, installable, workflow to take as input a series of MuSpAn domains, and return SHIFT scores quantifying changes in cell-cell interactions between biological states.

---

## Installation

Clone the repository and install locally:

```bash
git clone https://github.com/<username>/SHIFT.git
cd SHIFT
pip install -e .
```

---

## Quick Start

### Step 1: Generate a feature dataframe

Create a MuSpAn `Domain` containing cell coordinates and cell type labels, then calculate the spatial summary statistics used by SHIFT:

```python
import shift

df = shift.generate_feature_dataframe(
    domain,
    cell_label_name="Celltype"
)
```

This returns a single-row dataframe containing spatial descriptors derived from the supplied domain.

### Step 2: Combine multiple ROIs

SHIFT operates on datasets containing multiple regions of interest (ROIs), where each row corresponds to one ROI and one column identifies the biological state associated with that ROI.

For example, if you have MuSpAn domains named ```"Adenoma ROI 1", "Adenoma ROI 2",... "Carcinoma ROI 1", "Carcinoma ROI 2"...``` etc, we can construct a DataFrame with the relevant MuSpAn metrics, labelled by a unique domain name and two (or more) shared "State" values, here ("Adenoma" or "Carcinoma"):

```python
rows = []

for domain in domains:
    df_row = shift.generate_feature_dataframe(
        domain,
        cell_label_name="Celltype"
    )

    df_row["name"] = domain.name
    df_row["State"] = domain.name.split("_")[0]

    rows.append(df_row)

df = pd.concat(rows, ignore_index=True)
```

### Step 3: Calculate SHIFT

Once ROI-level statistics have been assembled into a dataframe, SHIFT scores can be calculated directly:

```python
SHIFT_scores = shift.compute_shift(
    df,
    comparator_state_column="State",
    metadata_columns=["name"]
)
```

SHIFT compares descriptor distributions between two biological states and returns SHIFT scores for each cell-cell interaction pair.

---

## Input Requirements

The minimum information required for each ROI is:

- Cell coordinates (`x`, `y`)
- Cell type annotations
- A biological state assignment (e.g. normal, adenoma, carcinoma)

If using ```shift.generate_feature_dataframe(...)``` to calculate summary statistics, then data should first be loaded into a MuSpAn `Domain` for each ROI, after which this preprocessing function will automatically generate the required summary statistics.

---

## Worked Example

A complete worked example is provided in:

```text
examples/example_workflow.py
```

This example demonstrates the full workflow:

1. Generate synthetic spatial datasets.
2. Create MuSpAn domains.
3. Calculate spatial summary statistics using `generate_feature_dataframe`.
4. Combine statistics from multiple ROIs into a single dataframe.
5. Calculate SHIFT scores using `compute_shift`. 【3-ca6840】

The example generates two synthetic biological states:

- **State 1:** cell types A and B exhibit local co-localisation.
- **State 2:** cell type B is spatially excluded from regions occupied by cell type A.

SHIFT identifies these differences in spatial organisation through the resulting interaction scores.

---

## Main Functions

### `generate_feature_dataframe()`

Generates a consistent set of ROI-level spatial descriptors from a MuSpAn domain. Currently supported descriptors include:

- Cell type counts
- Quadrat correlation matrices
- Wasserstein distances
- Cross pair correlation functions
- Persistent homology descriptors
- Topographical correlation map filtrations
- Adjacency permutation tests
- Morisita-Horn indices

Individual descriptor classes can be selected or configured using the `statistics` and `args_for_statistics` arguments.

### `compute_shift()`

Calculates SHIFT scores from ROI-level descriptor data.

The function:

1. Compares descriptor distributions between two biological states using Mann-Whitney tests.
2. Combines evidence across descriptors using Fisher's method.
3. Assigns directional scores describing whether interactions become more or less spatially ordered between states.
---

## Citation

If you use SHIFT, please cite:

```text
Bull J.A. et al.

Integrating spatial descriptors quantifies changing cell
interactions as colorectal neoplasia progresses.

Cancer Research (2026).
```

If you use the function ```shift.generate_feature_dataframe``` to calculate statistics from your data as MuSpAn domains, or otherwise use MuSpAn within your pipeline, please also cite the MuSpAn paper:

```text
Bull J.A., Moore J.W., Corry S.M., Lin M., Belnoue-Davis H.L., Mulholland-Illingworth E.J., Leedham S.J., and Byrne H.M.

MuSpAn: a toolbox for multiscale spatial analysis. 

Nature Communications (2026). 
https://doi.org/10.1038/s41467-026-75649-7
```

---

## Support

This repository is intended as a lightweight implementation accompanying the manuscript. The primary resource for new users is the worked example in `examples/example_workflow.py`, which demonstrates the complete workflow from annotated spatial data through to calculation of SHIFT scores.