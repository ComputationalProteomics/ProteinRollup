library(R6)
library(outliers)
library(tidyverse)

# TODO
# Outlier detection
# One hit wonders

ProteinRollup <- R6Class(
    public = list(
        protein_rollup_on_matrix = function(design_fp, data_fp, protein_col, peptide_col, sample_col, out_fp, rollup_func="rrollup",
            protein_col_name="Protein") {
            
            ddf <- read_tsv(design_fp, col_types=cols())
            raw_rdf <- read_tsv(data_fp, col_types=cols())
            rdf <- raw_rdf %>% filter(!is.na(UQ(as.name(protein_col))))
            
            trimmed_count <- nrow(raw_rdf) - nrow(rdf)
            if (trimmed_count != 0) {
                message("Trimmed ", trimmed_count, " proteins with missing IDs")
            }
            
            sdf <- rdf %>% dplyr::select(dplyr::one_of(ddf[[sample_col]]))
            
            protein_data <- rdf %>% dplyr::select(protein_col) %>% unlist()
            peptide_data <- rdf %>% dplyr::select(peptide_col) %>% unlist()

            message("Performing rollup for ", length(peptide_data), " peptides to ", length(unique(protein_data)), " proteins")
            
            if (rollup_func == "rrollup") {
                rollup <- self$rrollup
            }
            else {
                stop("Unknown rollup func: ", rollup_func)
            }
            
            prot_sdf <- self$protein_rollup(peptide_data, protein_data, as.matrix(sdf), rollup_func=rollup, protein_col_name=protein_col_name)
            write_tsv(prot_sdf, path=out_fp)
        },
        
        # Inspired by DANTE and pmartR
        protein_rollup = function(protein_ids, pep_mat, rollup_func=self$rrollup, protein_col_name="Protein") {
            # protein_rollup = function(pep_ids, protein_ids, pep_mat, rollup_func=self$rrollup, protein_col_name="Protein") {
                
            warning("protein_rollup: Use on your own risk, mostly untested")
            
            # pep_expression <- cbind(pep_ids, pep_mat)
            # prot_expression <- cbind(protein_ids, pep_mat)
            # protein_map <- cbind(protein_ids, pep_ids)
            
            # browser()
            
            unique_proteins <- unique(protein_ids)
            pep_counts <- list()
            pep_results <- list()
            
            for (protein in unique_proteins) {
                
                row_inds <- which(protein_ids == protein)
                current_peps_only <- pep_mat[row_inds, , drop=FALSE]
                
                scaled_peps <- rollup_func(current_peps_only)
                pep_results[[protein]] <- scaled_peps
                pep_counts[[protein]] <- nrow(current_peps_only)
            }

            annot_df <- data.frame(Protein=unique_proteins)
            colnames(annot_df) <- c(protein_col_name)
            cbind(annot_df, pep_count=pep_counts %>% unlist(), data.frame(do.call("rbind", pep_results)))
        },
        
        
        rrollup = function(peptides, combine_func=median, min_overlap=1) {
            # warning("Relatively untested, and no Grubbs outlier test performed")
            
            num_peps <- nrow(peptides)
            if (num_peps == 1) {
                protein_val <- unlist(peptides)
            }
            else {
                # 1. Select peptide with least number missing values as reference
                na_counts <- apply(is.na(peptides), 1, sum)
                least_na_pepindex <- which(na_counts == min(na_counts))
                
                # If multiple peptides have same number non-missing, select with highest median abundance
                if (length(least_na_pepindex) > 1) {
                    mds <- matrixStats::rowMedians(peptides, na.rm=TRUE)[least_na_pepindex]
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
                
                # TODO Outlier removal
                minPs <- 5
                pvalue <- 0.05
                x_scaled_nooutlier <- self$remove_outliers(x_scaled, minPs=minPs, pvalue=pvalue)
                
                # TODO center
                
                # 4. Calculate sample-wise medians across reference and scaled peptides
                protein_val <- apply(x_scaled_nooutlier, 2, combine_func, na.rm=TRUE)
                
                # TODO calculate rollup score
            }
            
            protein_val
        },
        
        # Based on InfernoRDN
        # Unable to target the reference peptide, is that correct?
        # https://github.com/PNNL-Comp-Mass-Spec/InfernoRDN/blob/master/Rscripts/Rollup/RRollup.R
        remove_outliers = function(pep_mat, minPs=5, pvalue_thres=0.05) {
            
            warning("Completely untested")
            
            xPeptideCount <- colSums(!is.na(pep_mat))
            
            browser()
            
            print("Outside")
            
            for (sample_i in 1:ncol(pep_mat)) {
                if (xPeptideCount[sample_i] >= minPs) {
                    
                    no_outliers_found <- FALSE
                    iterations <- 0
                    print("In remove_outliers")
                    
                    while (!no_outliers_found && iterations < 1000) {
                        iterations <- iterations + 1
                        
                        grubbs <- outliers::grubbs.test(pep_mat[, sample_i])
                        if ((grubbs$p.value < pvalue_thres) && (!is.nan(grubbs$statistic[2])) && grubbs$statistic[2] != 0) {
                            pep_mat[, sample_i] <- self$rm.outlier.1(pep_mat[, sample_i], fill=TRUE, select_func=median)
                            print("Number iterations: ", iterations)
                        }
                        else {
                            print("No outliers found")
                            no_outliers_found <- TRUE
                        }
                    }
                }
            }
            pep_mat
        },
        
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
                    # if (median) {
                    #     res[which(x == outlier.1(x, opposite))] <- median(x[-which(x == outlier.1(x, opposite))], na.rm = na.rm) 
                    # }
                    # else {
                    #     res[which(x == outlier.1(x, opposite))] <- mean(x[-which(x == outlier.1(x, opposite))], na.rm = na.rm)
                    # }
                    res
                }
            }
        },
        
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
        },
        
        zrollup = function(peptides, combine_func=median) {
            warning("Relatively untested, and no Grubbs outlier test performed")
            
            num_peps <- nrow(peptides)
            #res <- matrix(NA, nrow=1, ncol=ncol(peptides))
            
            # Compute mean and sd of peptides
            mds <- matrixStats::rowMedians(peptides, na.rm=TRUE)
            sds <- matrixStats::rowSds(peptides, na.rm=TRUE)
            
            # Scale peptide data as pep_scaled = (pep - median) / sd
            medians_mat <- matrix(mds, nrow=num_peps, ncol=ncol(peptides), byrow=FALSE)
            standiv_mat <- matrix(sds, nrow=num_peps, ncol=ncol(peptides), byrow=FALSE)
            
            proteins_scaled <- apply((peptides - medians_mat) / standiv_mat, 2, combine_func, na.rm=TRUE)
            proteins_scaled
        },
        
        
        
        # qrollup_thres: 0 - 1 value, peptides above threshold used for rollup
        qrollup = function(peptides, qrollup_thres, combine_func=median) {
            
            warning("Relatively untested, and no Grubbs outlier test performed")
            num_peps <- nrow(peptides)
            if (num_peps == 1) {
                protein_val <- unlist(peptides)
                peps_used <- 1
            }
            else {
                # Step 1: Subset peptides with abundance >= qrollup_threshold
                means <- Matrix::rowMeans(peptides, na.rm=TRUE)
                quantil <- quantile(means, probs=qrollup_thres, na.rm=TRUE)
                
                quality_peps <- peptides[means >= quantil, ]
                peps_used <- nrow(quality_peps)
                
                # Step 1b: If only 1 peptide, set protein value to that
                if (nrow(quality_peps) == 1) {
                    protein_val <- unlist(quality_peps)
                }
                # Step 2: Set protein abundance to mean / median of peptide abundances
                else {
                    protein_val <- apply(quality_peps, 2, combine_func, na.rm=TRUE)
                }
                
                protein_val
            }
        }
    ),
    private = list(
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
        },
        
        dante_remove.outliers = function(Data, minPs=5, pvalue=0.05)
            # internal function used by normalize.proteins()
            # Calculates the protein abundances after removing outliers
            # Depends on package "outliers"
        {
            # library(outliers)
            # ColNames = colnames(Data)
            xPeptideCount <- colSums(!is.na(Data))
            for (i in 1:dim(Data)[2])
            {
                if (xPeptideCount[i] >= minPs)
                {
                    repeat
                    {
                        grubbs <- grubbs.test(Data[,i]) # Grubb's test
                        if ( (grubbs$p.value < pvalue) && (!is.nan(grubbs$statistic[2])) &&
                             (grubbs$statistic[2] != 0) ) # pass the p-value cutoff
                        {
                            Data[,i] <- rm.outlier(Data[,i], fill=TRUE, median=TRUE)
                            # fill the outlier with the median value
                        }
                        else { break }
                    }
                }
            }
            return(Data)
        },
        
        rm.outlier = function(data, fill, median) {
            
        }
    )
)


protroll <- ProteinRollup$new()
pr <- protroll
print("Loading module to 'protroll' and 'pr'")
