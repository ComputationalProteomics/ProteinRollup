# Slightly reworked DANTE script
# https://github.com/PNNL-Comp-Mass-Spec/InfernoRDN/blob/master/Rscripts/Rollup/RRollup.R

#' Executes the protein rollup
#' 
#' @export
#' @return Produced protein matrix
#' @param protein_ids Vector containing protein IDs corresponding to each peptide
#' @param pep_mat Matrix containing peptide expression values
#' @param rollup_func What rollup function to apply
#' @param protein_col_name Name of protein column in input
#' @param get_debug_info Print debug information on the go
#' @param one_hit_wonders Whether proteins based in a single peptide should be retained
#' @param min_presence Value between 0 and 1. What fraction of values must be present in the most present peptide for the protein to be valid
#' @param min_overlap The number of non-missing values required from each peptide
#' @param debug_protein Send in a protein ID to print more extensive debug information
#' @import dplyr readr
#' @importFrom stats cor median weighted.mean
protein_rollup = function(protein_ids, pep_mat, rollup_func=rrollup, protein_col_name="Protein", get_debug_info=FALSE, one_hit_wonders=FALSE, min_presence=0, min_overlap=3, debug_protein=NA) {

    if (typeof(pep_mat) != "double") {
        warning("Input expected to be a 'double' matrix, found: ", typeof(pep_mat))
    }
    
    unique_proteins <- unique(protein_ids)
    pep_counts <- list()
    prot_results <- list()
    scaled_pep_results <- list()
    scaled_pep_results_noout <- list()
    rollup_scores <- list()
    low_presence_omitted <- 0
    passing_proteins <- c()
    
    debug_print <- FALSE
    
    for (protein in unique_proteins) {
        
        if (!is.na(debug_protein) && protein == debug_protein) {
            message("Found target protein: ", debug_protein)
            debug_print <- TRUE
        }
        
        row_inds <- which(protein_ids == protein)
        current_peps_only <- pep_mat[row_inds, , drop=FALSE]
        
        if (debug_print) {
            print(current_peps_only)
            print("Passing threshold")
            print(passes_presence_threshold(current_peps_only, min_presence, min_overlap, debug_print=TRUE))
        }
        
        if (passes_presence_threshold(current_peps_only, min_presence, min_overlap)) {
            rollup_list <- rollup_func(current_peps_only, get_debug_info=TRUE, min_overlap=min_overlap)
            prot_results[[protein]] <- rollup_list$protein
            scaled_pep_results[[protein]] <- rollup_list$scaled_peptide
            scaled_pep_results_noout[[protein]] <- rollup_list$scaled_peptide_nooutlier
            rollup_scores[[protein]] <- rollup_list$rollup_score

            pep_counts[[protein]] <- nrow(current_peps_only)
            passing_proteins <- c(passing_proteins, protein)
        }
        else {
            #message("Skipping protein ", protein, " due to low presence")
            low_presence_omitted <- low_presence_omitted + 1
        }
        
        if (debug_print) {
            stop("Debug print done")
        }
    }
    
    message(low_presence_omitted, " proteins omitted due to low protein coverage in peptides")

    annot_df <- data.frame(Protein=unique_proteins)
    colnames(annot_df) <- c(protein_col_name)
    
    # Generate the final annotated protein matrix
    out_df <- cbind(
        Protein=passing_proteins, 
        pep_count=pep_counts %>% unlist(), 
        rollup_score=rollup_scores  %>% unlist(),
        data.frame(do.call("rbind", prot_results), row.names = NULL)
    ) %>% arrange(.data$Protein)
    
    if (!one_hit_wonders) {
        out_df <- out_df %>% filter(.data$pep_count > 1)
    }
    
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
   
passes_presence_threshold <- function(peptides, presence_fraction, min_overlap, debug_print=FALSE, use_any_peptide_for_coverage=FALSE) {
    
    peptides_filtered <- peptides[rowSums(!is.na(peptides)) >= min_overlap, , drop=FALSE]
    
    if (debug_print) {
        print("Raw data")
        print(peptides)
        print("Filtered data")
        print(peptides_filtered)
        print("Numbers")
        print(rowSums(!is.na(peptides)))
        print(rowSums(!is.na(peptides)) >= 3)
        message("Total numbers of columns: ", ncol(peptides))
        print("use_any_peptide_for_coverage=FALSE")
        print(length(which(colSums(!is.na(peptides_filtered)) > 0)) / ncol(peptides_filtered))
        print("use_any_peptide_for_coverage=TRUE")
        print(max(rowSums(!is.na(peptides_filtered))) / ncol(peptides_filtered))
        print("Rowsums:")
        print(rowSums(!is.na(peptides_filtered)))
        print("New setup")
        print(length(which(colSums(!is.na(peptides_filtered)) > 0)) / ncol(peptides_filtered))
    }
    
    # Any peptide should match the coverage
    if (nrow(peptides_filtered) == 0){
        FALSE
    }
    else if (use_any_peptide_for_coverage) {
        # Any peptide counts for coverage
        length(which(colSums(!is.na(peptides_filtered)) > 0)) / ncol(peptides_filtered) >= presence_fraction
    }
    else {
        # Best peptide should match minimum coverage
        max(rowSums(!is.na(peptides_filtered))) / ncol(peptides_filtered) >= presence_fraction
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
        
        # Get rid of medians for sparse peptides, meaning these won't be weighted in the scaling calculations
        overlap_count <- rowSums(!is.na(ref_ratios))
        overlap_medians[which(overlap_count < min_overlap)] <- 0 # Should this not be assigned NA?
        
        # 3. Scale each peptide using the median differences in distance
        # Internal variations of peptides stay consistent, but focused around reference peptide
        x_scaled <- peptides + matrix(overlap_medians, nrow=num_peps, ncol=ncol(peptides))
        
        x_scaled_nooutlier <- remove_outliers(x_scaled, minPs=minPs, pvalue_thres=outlier_pval)
        
        # 4. Calculate sample-wise medians across reference and scaled peptides
        protein_val <- apply(x_scaled_nooutlier, 2, combine_func, na.rm=TRUE)
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
    
# Based on InfernoRDN
# Unable to target the reference peptide, is that correct?
# The outliers::grubbs.test produces NaN values if comparing a state with zero variance to a single value (i.e. c(1,1,2) or c(2,2,2,3))
remove_outliers = function(pep_mat, minPs=5, pvalue_thres=0.05) {
    
    xPeptideCount <- colSums(!is.na(pep_mat))
    for (sample_i in 1:ncol(pep_mat)) {

        if (xPeptideCount[sample_i] >= minPs) {
            
            no_outliers_found <- FALSE
            iterations <- 0
            
            while (!no_outliers_found) {
                iterations <- iterations + 1

                grubbs <- outliers::grubbs.test(pep_mat[, sample_i])
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

rollup.score = function(currPepSel, currProtSel, method) {

    N1 <- dim(currPepSel)[1]
    N2 <- dim(currPepSel)[2]
    pepCorr <- rep(numeric(0), N1)
    ws <- rep(numeric(0), N1)
    
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
    
    meanCorr <- stats::weighted.mean(pepCorr, ws, na.rm=TRUE) # Mean correlation value for each protein
    Penalty1 <- 1 - 1 / N1
    Score <- meanCorr * Penalty1
    
    out <- Score
    
    return(out)
}


# if (!interactive()) {
#     main()
# }
