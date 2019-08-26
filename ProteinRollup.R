library(R6)
library(outliers)
library(tidyverse)

# Inspired by DANTE and pmartR
protein_rollup = function(protein_ids, pep_mat, rollup_func=rrollup, protein_col_name="Protein", get_debug_info=FALSE) {
    # protein_rollup = function(pep_ids, protein_ids, pep_mat, rollup_func=self$rrollup, protein_col_name="Protein") {
        
    if (typeof(pep_mat) != "double") {
        warning("pep_mat is expected to consist of a numeric matrix. You can obtain this by omitting your annotation columns and executing: as.matrix(my_df)")
    }
    
    unique_proteins <- unique(protein_ids)
    pep_counts <- list()
    prot_results <- list()
    scaled_pep_results <- list()
    scaled_pep_results_noout <- list()
    rollup_scores <- list()
    
    for (protein in unique_proteins) {
        
        row_inds <- which(protein_ids == protein)
        current_peps_only <- pep_mat[row_inds, , drop=FALSE]
        
        rollup_list <- rollup_func(current_peps_only, get_debug_info=TRUE)
        prot_results[[protein]] <- rollup_list$protein
        scaled_pep_results[[protein]] <- rollup_list$scaled_peptide
        scaled_pep_results_noout[[protein]] <- rollup_list$scaled_peptide_nooutlier
        rollup_scores[[protein]] <- rollup_list$rollup_score
        
        pep_counts[[protein]] <- nrow(current_peps_only)
    }

    annot_df <- data.frame(Protein=unique_proteins)
    colnames(annot_df) <- c(protein_col_name)
    
    # Generate the final annotated protein matrix
    out_df <- cbind(
        Protein=annot_df, 
        pep_count=pep_counts %>% unlist(), 
        rollup_score=rollup_scores %>% unlist(), 
        data.frame(do.call("rbind", prot_results))
    ) %>% arrange(Protein)
    rownames(out_df) <- NULL
    
    if (!get_debug_info) {
        out_df
    }
    else {
        list(
            protein=out_df,
            pep_scale=do.call("rbind", scaled_pep_results),
            pep_scale_noout=do.call("rbind", scaled_pep_results_noout)
        )
    }
}
    
rrollup = function(peptides, combine_func=median, min_overlap=3, get_debug_info=FALSE, minPs=5, outlier_pval=0.05) {
    
    num_peps <- nrow(peptides)
    if (num_peps == 1) {
        x_scaled <- unlist(peptides)
        x_scaled_nooutlier <- unlist(peptides)
        protein_val <- unlist(peptides)
        rollup_score <- NA
    }
    else {
        # 1. Select peptide with least number missing values as reference
        na_counts <- apply(is.na(peptides), 1, sum)
        least_na_pepindex <- which(na_counts == min(na_counts))
        
        # If multiple peptides have same number non-missing, select with highest median abundance
        if (length(least_na_pepindex) > 1) {
            mds <- rowSums(peptides, na.rm=TRUE)[least_na_pepindex]
            # mds <- matrixStats::rowMedians(peptides, na.rm=TRUE)[least_na_pepindex]
            least_na_pepindex <- least_na_pepindex[which(mds == max(mds))]
        }
        
        # Extract the reference peptide values
        reference_pep_vals <- unlist(peptides[least_na_pepindex, ])
        
        # 2. Ratio all peptides to the reference (in log scale, so simple substraction)
        ref_ratios <- matrix(reference_pep_vals, nrow=num_peps, ncol=ncol(peptides), byrow=TRUE) - peptides
        
        # Calculate how far off the median for each rows non-NA values are from the reference
        overlap_medians <- matrixStats::rowMedians(ref_ratios, na.rm=TRUE)
        
        # TODO Get rid of medians for sparse peptides, meaning these won't be weighted in the scaling calculations
        overlap_count <- rowSums(!is.na(ref_ratios))
        overlap_medians[which(overlap_count < min_overlap)] <- 0 # Should this not be assigned NA?
        
        # 3. Scale each peptide using the median differences in distance
        # Internal variations of peptides stay consistent, but focused around reference peptide
        x_scaled <- peptides + matrix(overlap_medians, nrow=num_peps, ncol=ncol(peptides))
        
        x_scaled_nooutlier <- remove_outliers(x_scaled, minPs=minPs, pvalue=outlier_pval)
        
        # 4. Calculate sample-wise medians across reference and scaled peptides
        protein_val <- apply(x_scaled_nooutlier, 2, combine_func, na.rm=TRUE)
        
        # TODO calculate rollup score
        rollup_score <- rollup.score(x_scaled, protein_val, "pearson")
    }
    
    if (!get_debug_info) {
        protein_val
    }
    else {
        list(
            protein=protein_val,
            scaled_peptide=x_scaled,
            scaled_peptide_nooutlier=x_scaled_nooutlier,
            rollup_score=rollup_score
        )
    }
}
    
rollup.score = function(currPepSel, currProtSel, method) {
    N1 <- dim(currPepSel)[1]
    N2 <- dim(currPepSel)[2]
    pepCorr <- rep(numeric(0),N1)
    ws <- rep(numeric(0),N1)
    
    for(i in 1:N1)#finds correlation values between each peptide profile and calculated protein profile
    {
        ws[i] <- sum(!is.na(currPepSel[i,]))/N2
        if(method=="pearson")
            pepCorr[i] <- cor(as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs")
        else if(method=="kendall")
            pepCorr[i] <- cor(as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs",
                              method="kendall")
        else
            pepCorr[i] <- cor(as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs",
                              method="spearman")
    }
    #meanCorr <- mean(pepCorr, na.rm=TRUE) #mean correlation value for each protein
    meanCorr <- weighted.mean(pepCorr, ws, na.rm=TRUE) #mean correlation value for each protein
    Penalty1 <- 1-1/N1
    #Penalty2 <- sum(!is.na(currPepSel))/(N1*N2)
    #Score <- meanCorr * Penalty1 * Penalty2
    Score <- meanCorr * Penalty1
    
    out <- Score
    return(out)
}
    
    
# Based on InfernoRDN
# Unable to target the reference peptide, is that correct?
# https://github.com/PNNL-Comp-Mass-Spec/InfernoRDN/blob/master/Rscripts/Rollup/RRollup.R
remove_outliers = function(pep_mat, minPs=5, pvalue_thres=0.05) {
    
    xPeptideCount <- colSums(!is.na(pep_mat))
    for (sample_i in 1:ncol(pep_mat)) {

        if (xPeptideCount[sample_i] >= minPs) {
            
            no_outliers_found <- FALSE
            iterations <- 0
            
            while (!no_outliers_found && iterations < 1000) {
                iterations <- iterations + 1
                
                grubbs <- outliers::grubbs.test(pep_mat[, sample_i])
                # browser()
                if ((grubbs$p.value < pvalue_thres) && (!is.nan(grubbs$statistic[2])) && grubbs$statistic[2] != 0) {
                    pep_mat[, sample_i] <- rm.outlier.1(pep_mat[, sample_i], fill=TRUE, select_func=median)
                }
                else {
                    no_outliers_found <- TRUE
                }
            }
        }
    }
    pep_mat
}
    
# modified R outlier functions from "outlier" package to handle missing
rm.outlier.1 = function (x, fill = FALSE, select_func = median, opposite = FALSE, na.rm = TRUE) {
    if (is.matrix(x)) {
        apply(x, 2, rm.outlier.1, fill = fill, median = median, opposite = opposite, na.rm = na.rm)
    }
    else if (is.data.frame(x)) {
        as.data.frame(
            sapply(x, rm.outlier.1, fill = fill, median = median, opposite = opposite, na.rm = na.rm)
        ) 
    }
    else {
        res <- x
        if (!fill) {
            res[-which(x == outlier.1(x, opposite))]
        }
        else {
            res[which(x == outlier.1(x, opposite))] <- select_func(x[-which(x == outlier.1(x, opposite))], na.rm = na.rm)
            res
        }
    }
}
    
# modified R outlier functions to handle missing
outlier.1 = function (x, opposite = FALSE, logical = FALSE, na.rm = TRUE)
{
    if (is.matrix(x)) {
        apply(x, 2, outlier.1, opposite = opposite, logical = logical)
    }
    else if (is.data.frame(x)) {
        sapply(x, outlier.1, opposite = opposite, logical = logical)
    }
    else {
        if (xor(((max(x, na.rm = na.rm) - mean(x, na.rm = na.rm)) <
                 (mean(x, na.rm = na.rm) - min(x, na.rm = na.rm))), opposite)) {
            if (!logical) {
                min(x, na.rm = na.rm)
            }
            else {
                x == min(x, na.rm = na.rm)
            }
        }
        else {
            if (!logical) {
                max(x, na.rm = na.rm)
            }
            else {
                x == max(x, na.rm = na.rm)
            } 
        }
    }
}

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

rollup.score = function(currPepSel, currProtSel, method) {
    
    N1 <- dim(currPepSel)[1]
    N2 <- dim(currPepSel)[2]
    pepCorr <- rep(numeric(0),N1)
    ws <- rep(numeric(0),N1)
    
    # Finds correlation values between each peptide profile and calculated protein profile
    for(i in 1:N1) {
        ws[i] <- sum(!is.na(currPepSel[i,])) / N2
        
        if (method=="pearson") {
            pepCorr[i] <- cor(
                as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs"
            )
        }
        else if (method=="kendall") {
            pepCorr[i] <- cor(
                as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs", method="kendall"
            )
        }
        else {
            pepCorr[i] <- cor(
                as.vector(currPepSel[i,]), as.vector(currProtSel), use="pairwise.complete.obs", method="spearman"
            )
        }
    }
    
    meanCorr <- weighted.mean(pepCorr, ws, na.rm=TRUE) # Mean correlation value for each protein
    Penalty1 <- 1 - 1 / N1
    Score <- meanCorr * Penalty1
    
    out <- Score
    return(out)
}


# Main run, would like to separate out this to other file
main <- function() {
    
    argv <- parse_input_params()
    # source(argv$protein_rollup_path)
    raw_rdf <- read_tsv(argv$rdf_fp, col_types=cols())
    
    if (!argv$one_column_mode) {
        ddf <- read_tsv(argv$ddf_fp, col_types=cols())
        rdf <- raw_rdf %>% filter(!is.na(UQ(as.name(argv$protein_col))))
        sdf <- rdf %>% dplyr::select(dplyr::one_of(ddf[[argv$sample_col]])) %>% as.matrix()
        protein_data <- rdf %>% dplyr::select(argv$protein_col) %>% unlist()
    }
    else {
        sdf <- raw_rdf[, -1] %>% as.matrix()
        protein_data <- raw_rdf[, 1] %>% unlist()
    }
    
    message("Performing rollup for ", nrow(raw_rdf), " peptides for ", length(unique(protein_data)), " unique protein IDs")
    
    prot_sdf <- protein_rollup(protein_data, as.matrix(sdf))
    write_tsv(prot_sdf, path=argv$out_fp)
    message("Resulting file written to ", argv$out_fp)
}

parse_input_params <- function() {
    parser <- arg_parser("Protein rollup, standalone wrapper")
    parser <- add_argument(parser, "--ddf_fp", help="Design matrix path", type="character")
    parser <- add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character")
    parser <- add_argument(parser, "--out_fp", help="Output matrix path", type="character")
    
    parser <- add_argument(parser, "--one_column_mode", help="If first column is the only annotation column, the --ddf_fp, --sample_col and --protein_col argument can be omitted", type="boolean", default=FALSE)
    
    parser <- add_argument(parser, "--sample_col", help="Design matrix sample column", type="character", default="sample")
    parser <- add_argument(parser, "--protein_col", help="Protein column in main data frame", type="character", default="Protein")
    
    parser <- add_argument(parser, "--out_protein_name", help="Name of protein column in output", type="character")
    # parser <- add_argument(parser, "--protein_rollup_path", help="CraftOmics protein tools path", type="character", default="ProteinRollup.R")
    
    argv <- parse_args(parser)
}

if (!interactive()) {
    library(argparser)
    main()
}
