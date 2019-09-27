library(outliers)
suppressPackageStartupMessages(library(tidyverse, warn.conflicts=FALSE))
library(matrixStats, warn.conflicts=FALSE)

# Inspired by DANTE and pmartR
# https://github.com/PNNL-Comp-Mass-Spec/InfernoRDN/blob/master/Rscripts/Rollup/RRollup.R

protein_rollup = function(protein_ids, pep_mat, rollup_func=rrollup, protein_col_name="Protein", get_debug_info=FALSE, one_hit_wonders=FALSE, min_presence=0, min_overlap=3) {
    # protein_rollup = function(pep_ids, protein_ids, pep_mat, rollup_func=self$rrollup, protein_col_name="Protein") {

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
    
    for (protein in unique_proteins) {
        
        row_inds <- which(protein_ids == protein)
        current_peps_only <- pep_mat[row_inds, , drop=FALSE]
        
        if (passes_presence_threshold(current_peps_only, min_presence)) {
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
    }
    
    message(low_presence_omitted, " proteins omitted due to low presence")

    annot_df <- data.frame(Protein=unique_proteins)
    colnames(annot_df) <- c(protein_col_name)
    
    # Generate the final annotated protein matrix
    out_df <- cbind(
        Protein=passing_proteins, 
        pep_count=pep_counts %>% unlist(), 
        rollup_score=rollup_scores  %>% unlist(),
        data.frame(do.call("rbind", prot_results), row.names = NULL)
    ) %>% arrange(Protein)
    
    if (!one_hit_wonders) {
        out_df <- out_df %>% filter(pep_count > 1)
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
   
passes_presence_threshold <- function(peptides, presence_fraction) {
    
    if (nrow(peptides) > 1) {
        message(paste(peptides, collapse=", "), " ", max(rowSums(!is.na(peptides))) / ncol(peptides))
    }
    
    max(rowSums(!is.na(peptides))) / ncol(peptides) >= presence_fraction
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
        
        x_scaled_nooutlier <- remove_outliers(x_scaled, minPs=minPs, pvalue=outlier_pval)
        
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
    
    if (argv$min_presence < 0 || argv$min_presence > 1) {
        stop("--min_presence expected to be in range 0 to 1, found: ", argv$min_presence)
    }
    
    if (!argv$one_column_mode) {
        ddf <- read_tsv(argv$ddf_fp, col_types=cols())
        rdf <- raw_rdf
        sdf <- rdf %>% dplyr::select(dplyr::one_of(ddf[[argv$sample_col]])) %>% as.matrix()
        protein_data <- rdf %>% dplyr::select(argv$protein_col) %>% unlist()
        protein_data[is.na(protein_data)] <- "NA"
    }
    else {
        sdf <- raw_rdf[, -1] %>% as.matrix()
        protein_data <- raw_rdf[, 1] %>% unlist()
        protein_data[is.na(protein_data)] <- "NA"
    }
    
    message("Performing rollup for ", nrow(raw_rdf), " peptides for ", length(unique(protein_data)), " unique protein IDs")
    
    prot_adf <- protein_rollup(
        protein_data,
        as.matrix(sdf),
        one_hit_wonders=argv$one_hit_wonders,
        min_presence=argv$min_presence,
        min_overlap=argv$min_overlap
    )

    if (argv$out_protein_name != "") {
        prot_adf <- cbind(prot_adf[["Protein"]], prot_adf)
        colnames(prot_adf) <- c(argv$out_protein_name, colnames(prot_adf)[-1])
    }

    write_tsv(prot_adf, path=argv$out_fp)
    message("Resulting file written to ", argv$out_fp)
}

parse_input_params <- function() {
    parser <- argparser::arg_parser("Protein rollup, standalone wrapper")
    parser <- argparser::add_argument(parser, "--ddf_fp", help="Design matrix path", type="character", default=NA)
    parser <- argparser::add_argument(parser, "--rdf_fp", help="Raw matrix path", type="character", nargs=1)
    parser <- argparser::add_argument(parser, "--out_fp", help="Output matrix path", type="character", nargs=1)
    parser <- argparser::add_argument(parser, "--sample_col", help="Design matrix sample column", type="character", default=NA)
    parser <- argparser::add_argument(parser, "--protein_col", help="Protein column in main data frame", type="character", default=NA)
    
    parser <- argparser::add_argument(parser, "--one_column_mode", help="If first column is the only annotation column, the --ddf_fp, --sample_col and --protein_col argument can be omitted", type="boolean", default=FALSE)
    
    parser <- argparser::add_argument(parser, "--min_overlap", help="Min. shared overlap between reference peptides and other", type="numeric", default=3)
    parser <- argparser::add_argument(parser, "--one_hit_wonders", help="Should proteins with a single peptide as support be included", type="boolean", default=FALSE)
    parser <- argparser::add_argument(parser, "--min_presence", help="Skip proteins where no peptide is present in at least this percentage", type="numeric", default=0)
    parser <- argparser::add_argument(parser, "--out_protein_name", help="Name of protein column in output", type="character", default="")
    # parser <- add_argument(parser, "--protein_rollup_path", help="CraftOmics protein tools path", type="character", default="ProteinRollup.R")

    parser <- argparser::add_argument(parser, "--show_warnings", help="Immediately print warnings", type="bool", default=FALSE)

    argv <- argparser::parse_args(parser)
    
    if (!argv$one_column_mode && (is.na(argv$sample_col) || is.na(argv$protein_col) || is.na(argv$ddf_fp))) {
        stop("If not running one_column_mode both --ddf_fp, --sample_col and --protein_col needs to be supplied")
    }
    
    if (argv$show_warnings) {
        options(warn=1)
        message("warn=1 assigned, warnings shown as they occur")
    }
    
    argv
}

if (!interactive()) {
    library(argparser)
    main()
}
