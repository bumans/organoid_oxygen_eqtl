library(tidyverse)
library(Matrix)
library(RcppParallel)
library(mashr)
library(udr)
library(matrixStats)

args = commandArgs(trailingOnly=TRUE)
celltype <- as.character(args[1])
filepath <- as.character(args[2])

print(celltype)

set.seed(1)

input <- readRDS(file = paste0(filepath, "MatrixEQTLSumStats_", celltype, "only.mash.rds"))

# estimate correlation for measurement overlap
data.temp = mash_set_data(input$random.b, input$random.s)
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

# https://stephenslab.github.io/udr/reference/ud_init.html
# Use exchangeable zscore model by setting alpha = 1.
data.random = mash_set_data(input$random.b, input$random.s, V=Vhat, alpha = 0)
data.strong = mash_set_data(input$strong.b, input$strong.s, V=Vhat, alpha = 0)
# why do we use strong.b and not strong.z if alpha=1?
U.c = cov_canonical(data.random)

V.em = mash_estimate_corr_em(data.random, Ulist = c(U.c), details = TRUE)
data.random <- mash_update_data(mashdata = data.random, V = V.em$V, ref = NULL)
data.strong <- mash_update_data(data.strong, V = V.em$V, ref = NULL)

fit0 = ud_init(data.strong, n_rank1 = 0, n_unconstrained = 20)
fit1 = ud_fit(fit0, control = list(unconstrained.update = "ed", maxiter  = 1e3), verbose=TRUE)
U.ed <- lapply(fit1$U,function (e) "[["(e,"mat"))

U.c = cov_canonical(data.random)

m = mash(data.random, Ulist = c(U.ed, U.c), outputlevel = 1)

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

saveRDS(m2, file = paste0(filepath, "MatrixEQTLSumStats_", celltype, "only_udr_yunqi_vhatem_EE.rds"))
