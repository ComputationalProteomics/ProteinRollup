# zrollup = function(peptides, combine_func=median) {
#     warning("Relatively untested, and no Grubbs outlier test performed")
#     
#     num_peps <- nrow(peptides)
#     #res <- matrix(NA, nrow=1, ncol=ncol(peptides))
#     
#     # Compute mean and sd of peptides
#     mds <- matrixStats::rowMedians(peptides, na.rm=TRUE)
#     sds <- matrixStats::rowSds(peptides, na.rm=TRUE)
#     
#     # Scale peptide data as pep_scaled = (pep - median) / sd
#     medians_mat <- matrix(mds, nrow=num_peps, ncol=ncol(peptides), byrow=FALSE)
#     standiv_mat <- matrix(sds, nrow=num_peps, ncol=ncol(peptides), byrow=FALSE)
#     
#     proteins_scaled <- apply((peptides - medians_mat) / standiv_mat, 2, combine_func, na.rm=TRUE)
#     proteins_scaled
# }

# # qrollup_thres: 0 - 1 value, peptides above threshold used for rollup
# qrollup = function(peptides, qrollup_thres, combine_func=median) {
#     
#     warning("Relatively untested, and no Grubbs outlier test performed")
#     num_peps <- nrow(peptides)
#     if (num_peps == 1) {
#         protein_val <- unlist(peptides)
#         peps_used <- 1
#     }
#     else {
#         # Step 1: Subset peptides with abundance >= qrollup_threshold
#         means <- Matrix::rowMeans(peptides, na.rm=TRUE)
#         quantil <- quantile(means, probs=qrollup_thres, na.rm=TRUE)
#         
#         quality_peps <- peptides[means >= quantil, ]
#         peps_used <- nrow(quality_peps)
#         
#         # Step 1b: If only 1 peptide, set protein value to that
#         if (nrow(quality_peps) == 1) {
#             protein_val <- unlist(quality_peps)
#         }
#         # Step 2: Set protein abundance to mean / median of peptide abundances
#         else {
#             protein_val <- apply(quality_peps, 2, combine_func, na.rm=TRUE)
#         }
#         
#         protein_val
#     }
# }