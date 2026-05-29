# Project Information

-   Analysis done on behalf of the MGH KFCCR Tumor Cartography Core for the Sgroi Lab (MGH).
-   Only code and results will be placed here. No data.
-   RNA: 380-plex preset panel + 100 custom genes (480 genes total).
-   Protein: multiplex IF staining for CD4, CD20, CD8A, CD3E, Vimentin, alphaSMA, CD45, and E-Cadherin.

# Analysis Summary

## 00 - 03

Just setting up the environment and the data for the project.

## 04 - 08

Iterative processing of the data.

-   First, I processed the data without any QC to identify possible low-expressing cell populations that might be filtered out by mistake.
-   This gives me the chance to "whitelist" these cells.
-   In this project, we did not find any such "whitelist" populations.
-   I performed QC with `scuttle`.
-   Integration (batch was set as the sample ID) and dimensionality reduction were achieved with `scvi-tools`.
-   Clustering performed with `BPCells`'s implementation of the leiden algorithm, run on an 50-NN adjacency matrix, weighted with SNN.
-   A resolution of 0.35 was chosen, after running the algorithm at many different resolutions and inspecting how well the results captured canonical marker gene expression patterns on the UMAP embedding.
-   The leiden clusters were manually annotated, after finding marker genes for each with the `glmm` branch of the `presto` R package.
-   The object produced was `sgroi-tnbc_filtered.h5ad`. Uncompressed size: 4.2 GB. Shape: `n_obs × n_vars = 1051558 × 480`.

## 09 - 10

Sub-clustering and re-annotating of the immune cells.

-   The data was subset to only include the immune cells. Done with the following line of code:

```         
adata = adata[adata.obs.celltype_level1.isin(["Lymphoid", "Myeloid", "pDC", "mregDC", "B", "Plasma"]),].copy()
```

-   Processing, clustering, and annotating were all done again (resolution = 0.45).
-   The object produced was `sgroi-tnbc_filtered_immune_subset.h5ad`. Uncompressed size: 1.3 GB. Shape: `n_obs × n_vars = 359467 × 480`.

## 11 - 12

Attempt to further annotate the T cells.

-   Similar methods to those in steps 09 - 10.
-   Results did not make much sense. This can happen in ST data due to contaminating signal.

## 13 - 14

First two rounds of DE (differential expression) analysis.

-   Done before we had totally sorted out the granular T cell typing. Updates will follow.
-   Mixed modeling done with `nebula`, setting patient ID as a random effect.
-   This is where we defined the HOXB13+/- patient groups, using qPCR data.

## 15

First round of DA (differential abundance) analysis.

-   Done before we had totally sorted out the granular T cell typing.

## 16 - 17

Attempt to further annotate the macrophages.

-   Similar methods to those in steps 09 - 10.
-   Results did not make much sense. This can happen in ST data due to contaminating signal.

## 18

Using `HieraType` to assign granular cell type annotations to the T cells and to the myeloid cells.

-   This method seemed to perform well.

## 19

Second round of DA analysis.

-   Done after figuring out the granular T cell typing.

## 20

Identifying spatial tissue niches with `Banksy`.

-   Note that I used the `Banksy` package to produce the necessary matrix and to do PCA.
-   I ran `Harmony` on this PCA matrix to integrate across samples (batch = sample ID).
-   I used `BPCells` to produce a 50-SNN graph and then clustered with the leiden algorithm (resolution = 0.35).
-   There were 4 singletons, which is fine. We will ignore those 4 cells during downstream analyses.

## 21

Organized all the data for the project.

-   I added globally-defined spatial coordinates.
-   I touched up the cell metadata so that the cell-typing makes more sense. `celltype_final_amended_1` holds the up-to-date cell-typing.
-   After using `HieraType` on the Myeloid cells, I added to the cell metadata again. `celltype_final_amended_2` holds the up-to-date cell-typing.

## 22

Characterizing the spatial tissue niches.

## 23

Finding spatially-correlated gene modules.

## 24

Performing NE (neighborhood enrichment) analysis with `cellcharter`.

## 25

Third round of DA analysis.

-   Done after figuring out the granular myeloid typing.
