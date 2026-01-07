rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

# The goal here is to make an AnnData (.h5ad file) object for each TMA with 
# the data that we want to consider for downstream analysis.

# We will not do QC, but we will subset down to the cells that we could be sure 
# were actually part of certain cores

## Read the metadata for each TMA
meta_list <- qs2::qs_read("01_core_mapping_outs/core-mapped_metadata.qs2")

## Paths to the counts matrices
h5mats <- dir(
  path = "../../sgroi-tnbc-data", 
  pattern = "cell_feature_matrix.h5", 
  recursive = T, 
  full.names = T
)

# Function that will read the data that we want (counts and feature metadata),
# keep cells that exist in all parts of the data, and then create an AnnData 
# object. I have also added the protein intensity data to the cell metadata.

# Per the 10X website: 
# The cell_feature_matrix includes scaled mean protein intensity values 
# (see calculation), 10x Genomics Xenium Protein marker IDs (TXP), and protein 
# marker abbreviations.
# Link to calculation: https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/algorithms-overview/protein-algorithms#mean-intensity

create_anndata_object <- function(tma_name) {
  
  # Read the data from h5
  h5mat_oi <- h5mats[grepl(pattern = tma_name, x = h5mats)]
  cts_tmp <- rhdf5::h5read(file = h5mat_oi, name = "matrix")
  
  # Get the mask for the RNA data
  gene_mask <- cts_tmp$features$feature_type == "Gene Expression"
  
  # Create the sparse counts matrix
  m <- Matrix::sparseMatrix(x = cts_tmp$data, 
                            dims = cts_tmp$shape, 
                            i = (cts_tmp$indices+1), 
                            p = cts_tmp$indptr, 
                            repr = "C", 
                            dimnames = list(cts_tmp$features$name, cts_tmp$barcodes)
  )
  
  # Subset to RNA and make sure the cells are the same across metadata and count
  m_c <- m[which(gene_mask), rownames(meta_list[[tma_name]])]
  
  # Get the mask for the protein data
  pr_mask <- cts_tmp$features$feature_type == "Protein Expression"
  
  # Subset to protein and make sure the cells are the same across metadata and intensity
  m_i <- m[which(pr_mask), rownames(meta_list[[tma_name]])]
  m_i <- as.matrix(m_i) |> as.data.frame() |> t()
  
  # Create the gene metadata & subset down to RNA
  genemeta <- lapply(cts_tmp$features, function(x) x[gene_mask])
  genemeta <- dplyr::bind_cols(genemeta) |> dplyr::select(name, id, feature_type) |> as.data.frame()
  rownames(genemeta) <- genemeta$name
  
  # Create the AnnData object
  adata <- anndataR::AnnData(
    X = Matrix::t(m_c), 
    obs = cbind(meta_list[[tma_name]], m_i),
    var =  genemeta
  )
  
  # Return the AnnData object
  return(adata)
  
}

## Make an object for each TMA
adata_list <- sapply(X = names(meta_list), FUN = create_anndata_object)

## Looks good... needed to check that the protein data was there.
adata_list$tma1$obs |> dplyr::glimpse()

## Save each AnnData object
filenames <- paste(names(adata_list), "h5ad", sep = ".")
purrr::map2(.x = adata_list, 
            .y = filenames, 
            .f = anndataR::write_h5ad, 
            compression = "none", 
            mode = "w")

