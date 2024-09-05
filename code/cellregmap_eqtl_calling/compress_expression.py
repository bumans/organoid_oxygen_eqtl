import scanpy as sc
import numpy as np
import pandas as pd
import xarray as xr

pseudocell_expression_loc = snakemake.input['pseudocell_data']
compressed_expression_loc = snakemake.output['exp']

## import the normalized pseudocell counts
phenotype = pd.read_csv(pseudocell_expression_loc, sep='\t')

#Save phenotype as binary file for fast read-in
phenotype_export = xr.DataArray(
    data=phenotype.values,
    dims=["pseudo_cell", "gene"],
    coords={"pseudo_cell": phenotype.index.values, "gene": phenotype.columns.values},
    name="expression"
    )
# assert all(phenotype.pseudo_cell.values == sample_mapping.index.values)

phenotype.to_netcdf(compressed_expression_loc)
