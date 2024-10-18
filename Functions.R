#install.packages("rstudioapi")
library(rstudioapi)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape")
library(reshape)
#install.packages("reshape2")
library(reshape2)
#install.packages("stringi")
library(stringi)

#Function that normalizes protein abundances and creates a small QC report
protein.Report.QC <- function(data, normalization = "median"){
  
  #INPUT
  #data: data frame of protein expressions.
  #      First column should be the protein accessions.
  #      Abundance columns should start after the 4th column (VIQoR output).
  #      Quantitative columns should be named properly.
  #      For example: 4th replicate of the 1st condition should be stated as C_1_R_4
  #normalization: method to use for protein expression normalization. Available options: meadian or mean
  
  #OUTPUT
  #QC plots/QC_proteins.pdf: a small QC report on protein level in .pdf format. It contains:
  #                          1) Bar plot of missing values per condition.
  #                          2) Density plot of abundances per condition after normalization.
  #                          3) Violin plot of samples after normalization.
  #                          4) PCA plot of samples after normalization.
  #data: a data frame of normalized protein expressions.
  
  if(missing(data)){
    
    stop("Error: Protein Expression data missing!")
    
  }
  
  #Create directory for QC plots.
  if(all(!list.files() == "QC plots")){ dir.create("QC plots")}
  
  #Count number of missing values per sample and create a barplot.
  NAs <- data.frame(Samples = colnames(data)[-c(1:4)], value = colSums(is.na(data[,c(5:ncol(data))])))
  NAs$Conditions <- substring(as.character(NAs$Samples), 3, 3)
  NAs$Samples <- factor(NAs$Samples, levels = NAs$Samples)
  
  p.NAs <- ggplot2::ggplot(NAs, ggplot2::aes(x = Samples, y = value, fill = Conditions)) + 
           ggplot2::geom_bar(stat = "identity", position = "dodge", ggplot2::aes(color = Samples))  +
           ggplot2::geom_text(ggplot2::aes(label = value, y = value), position = ggplot2::position_dodge(1), vjust = -1) +
           ggplot2::ylab("Missing values count") +
           ggplot2::xlab("Samples") +
           ggplot2::ggtitle("Number of missing values per sample") +
           ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
           ggplot2::guides(color = FALSE)
  

  #Normalization and scaled density plot per condition.
  if(normalization == "median"){
    
    data[ ,c(5:ncol(data))] <- data.frame(t(t(data[ ,c(5:ncol(data))]) - apply(data[,c(5:ncol(data))], 2, median, na.rm = TRUE)))
    
  } else if(normalization == "mean"){
    
    data[ ,c(5:ncol(data))] <- data.frame(t(t(data[ ,c(5:ncol(data))]) - apply(data[,c(5:ncol(data))], 2, mean, na.rm = TRUE)))
    
  } else {
    
    stop("Error: normalization method is not valid!")
    
  }
  
  exprs.long <- reshape::melt(data[ ,-c(1:4)], id.vars = NULL)
  exprs.long$Conditions <- substring(as.character(exprs.long$variable), 3,3)
  
  #Distribution plot per condition
  p.Dist <- ggplot2::ggplot(na.omit(exprs.long), ggplot2::aes(x = value, color = Conditions, ..scaled..)) + 
            ggplot2::geom_density() +
            ggplot2::ylab("Scaled density") +
            ggplot2::xlab("Protein expression (log2)") +
            ggplot2::ggtitle(paste0("Distributions of normalized (", normalization,") protein expressions per condition")) +
            ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  #Violin plot combined with box plot of protein expressions per sample.
  p.Violin <- ggplot2::ggplot(na.omit(exprs.long), ggplot2::aes(x = variable, y = value)) +
              ggplot2::geom_violin(width = 1, ggplot2::aes(fill = Conditions)) +
              geom_boxplot(width = 0.1, ggplot2::aes(color = variable)) +
              ggplot2::ylab("Protein expression (log2)") +
              ggplot2::xlab("Samples") +
              ggplot2::ggtitle(paste0("Violin plot of normalized (", normalization,") protein expressions per condition")) +
              ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
              ggplot2::guides(color = FALSE, fill = FALSE)
  
  #Protein feature CV per condition.
  
  cv.df <- 2^data[,5:ncol(data)]
  
  conditions <- unique(substring(colnames(cv.df), 1,3))
  conditions.B <- lapply(conditions, function(x) grepl(pattern = x, colnames(cv.df)))
  conditions.CV <- lapply(conditions.B, function(x) 100 * apply(cv.df[,x], 1, sd, na.rm = TRUE)/rowMeans(cv.df[,x], na.rm = T))
  conditions.CV <- as.data.frame(do.call(cbind, conditions.CV))
  colnames(conditions.CV) <- conditions
  
  cv.df <- reshape::melt(conditions.CV, id.vars = NULL)
  
  annotation <- rbind(round(colMeans(conditions.CV, na.rm = T), 2), round(apply(conditions.CV, 2, median, na.rm = T), 2))
  annotation <- cbind(c("Mean CV", "Median CV"), annotation)
  colnames(annotation)[1] <- ""
  
  p.CV <- ggplot2::ggplot(data = na.omit(cv.df), mapping = ggplot2::aes(x = variable, y = value)) + 
          ggplot2::geom_boxplot(ggplot2::aes(color = variable)) +
          ggplot2::xlab("Conditions") +
          ggplot2::ylab("Protein CV (%)") +
          ggplot2::ggtitle(paste0("Per condition CV boxplot (total mean: ", round(mean(cv.df$value, na.rm = T), 2), "%, total median: ", round(median(cv.df$value, na.rm = T), 2))) +
          ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(50, 50, 50, 50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
          ggplot2::ylim(0, max(cv.df$value, na.rm = T) + 20) +
          ggplot2::annotation_custom(gridExtra::tableGrob(annotation,theme = ttheme_default(base_size = 10), rows = NULL), xmin = -Inf, xmax = Inf, ymin = max(cv.df$value, na.rm = T))
  
  #PCA analysis and PCA plot of the two first components.
  pca <- stats::prcomp(~ ., data[ ,-c(1:4)], na.action = na.omit)
  
  pca.df <- data.frame(pca$rotation[, 1:2], rownames(pca$rotation), substr(rownames(pca$rotation), start = 1, stop = 3), stringsAsFactors = F)
  
  colnames(pca.df) <- c("x", "y", "Samples", "Conditions")
  
  p.PCA <- ggplot2::ggplot(data = pca.df, ggplot2::aes(x = x, y = y, color = Conditions, fill = Samples)) +
           ggplot2::geom_point(size = 3, shape = 21, stroke = 2) +
           ggplot2::xlab(paste0("Component ", 1, " - ", round(pca$sdev[1]^2/(sum(pca$sdev^2)), digits = 2)*100, " %" )) +
           ggplot2::ylab(paste0("Component ", 2, " - ", round(pca$sdev[2]^2/(sum(pca$sdev^2)), digits = 2)*100, " %" )) + 
           ggplot2::ggtitle(paste0("PCA plot of per sample relative protein expressions (Components ", 1, " and ", 2, ")")) +
           ggplot2::theme(plot.margin = ggplot2::margin(50, 50, 50, 50))
      
  
  #Save pdf file
  ggsave("QC plots/QC_proteins.pdf", plot = grid.arrange(p.NAs, p.Dist, p.Violin, p.CV, p.PCA, ncol = 2, nrow = 3, top = "Quality control plots for protein expression data"),
         width = 22, height = 30, limitsize = FALSE)
  
  #Return the normalized protein expressions.
  return(data)
  
}

#Function to filter regulated features.
regulated.Features <- function(data, FC.c, FC.thresholds, Stats.c, Stats.threshold, label){
  
  #INPUT
  #data: data frame of log2 FC ratios of protein expressions and a statistics for each comparison (PolySTest data)
  #FC.c: start and end column indices for the FC columns in the data frame.
  #FC.thresholds: lower and upper thresholds for FC.
  #Stats.c: start and end column indices for the statistical score columns in the data frame.
  #Stats.threshold: upper limit of the statistical score.
  #label: a label for the barplot describing the features in the data.
  
  #OUTPUT
  #data: filtered data according to the FC and statistical score to at least one comparison.
  #Figures/Significant_#label: a barplot to count regulated features per comparison
  
  arguments <- setdiff(ls(), names(match.call())[-1])
  
  if(length(arguments) > 0){

    stop(paste0("Error: missing argument(s) -> ", paste0(arguments, collapse = ", ")))

  }
  
  #Create directory for QC plots.
  if(all(!list.files() == "Figures")){ dir.create("Figures")}
  
  #Selection by FC
  FC.data <- data[, FC.c[1]:FC.c[2]]
  FC.data <- as.data.frame((FC.data <= FC.thresholds[1]) | (FC.data >= FC.thresholds[2]))
  FC.data$'In at least one condition (FC)' <- apply(FC.data, 1, any)
  FC.data[is.na(FC.data)] <- FALSE

  #Selection by q-value
  Stats.data <- data[, Stats.c[1]:Stats.c[2]]
  Stats.data <- as.data.frame(Stats.data <= Stats.threshold)
  Stats.data$'In at least one condition (q-value)' <- apply(Stats.data, 1, any)
  Stats.data[is.na(Stats.data)] <- FALSE
  
  #Merge and prepare to plot
  to.Plot <- cbind(FC.data[, 1:(ncol(FC.data)-1)], Stats.data[, 1:(ncol(Stats.data)-1)])
  to.Plot$'Selected by FC' <- FC.data[, ncol(FC.data)]
  to.Plot$'Selected by q-value' <- Stats.data[, ncol(Stats.data)]
  to.Plot$'Total selected' <- ((to.Plot$'Selected by FC' + to.Plot$'Selected by q-value') == 2)
  
  data$Selected <- to.Plot$'Total selected'
  
  to.Plot <- data.frame(variable = colnames(to.Plot), value = colSums(to.Plot), group = c(rep("FC", ncol(FC.data)-1), rep("q-value", ncol(Stats.data)-1), rep("Total", 3)))
  to.Plot$variable <- factor(to.Plot$variable, levels = to.Plot$variable)
  to.Plot$group <- factor(to.Plot$group, levels = unique(to.Plot$group))
  
  p <- ggplot2::ggplot(to.Plot, ggplot2::aes(x = variable, y = value, fill = variable)) +
       ggplot2::geom_bar(stat = "identity", position = "dodge") +
       ggplot2::geom_text(ggplot2::aes(label = value, y = value), position = ggplot2::position_dodge(1), vjust = -1) +
       ggplot2::facet_grid(~group, switch = "x", scales = "free_x", space = "free_x") +
       ggplot2::ylab(paste0("Number of ", label, " features")) +
       ggplot2::xlab("Labels") +
       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
       ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(50, 50, 50, 50)) +
       ggplot2::ggtitle(paste0("Selected ", label, " features by FC and q-value threshold"))

  #Save plot.
  ggsave(paste0("Figures/Significant_", label, "_Features.pdf"), plot = p , width = 10, height = 9, limitsize = FALSE)
  
  return(data)

}

#Function to normalize PTMs abundances as a ratio to protein expressions.
PTMs.Normalization <- function(peptides, proteins, fasta, quantitative.columns, normalization = "median"){
  
  #INPUT
  #peptides: data frame of log2 tranformed modified peptide abundances. First column should contrain the modified peptide sequence.
  #          The quantitative columns should be placed at the end.
  #proteins: data frame of log2 protein expressions. First column should contain the protein group accessions.
  #          The quantitative columns should be placed at the end.
  #fasta: a fasta file to provide protein sequences as a list (protr::readFASTA())
  #quantitative.columns: integer that determines the number of quantitative columns (should be the same for both peptide and protein data)
  #normalization: peptide normalization method before and after normalization by protein expression.(median or mean)
  
  #OUTPUT
  #A list that contains:
  #mapped: normalized modified peptides that can be mapped to the quantified protein sequences.
  #unmapped: modified peptides that could not be mapped on the the quantified protein sequences.

  arguments <- setdiff(ls(), names(match.call())[-1])
  
  if(length(arguments) > 0){
    
    stop(paste0("Error: missing argument(s) -> ", paste0(arguments, collapse = ", ")))
    
  }
  
  #Fasta to df
  fasta <- data.frame(protein = sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), sequence = unlist(fasta), stringsAsFactors = F)
  
  #Redude the proteins in Fasta.
  quantified.proteins <- unlist(strsplit(proteins[ ,1], "\\|"))
  fasta <- fasta[unlist(lapply(quantified.proteins, function(x) which(fasta$protein == x))), ]
  
  #Map counterpart sequences on the quantified protein sequences.
  modified.peptides.mapping <- lapply(peptides$Counterpart, function(x) fasta$protein[which(stringi::stri_detect_fixed(fasta$sequence, pattern = x))])
  
  
  #Normalization of the modified peptide abundances.
  if(normalization == "median"){
    
    peptides[ ,(ncol(peptides) - quantitative.columns + 1):ncol(peptides)] <- data.frame(t(t(peptides[ ,(ncol(peptides) - quantitative.columns + 1):ncol(peptides)]) - apply(peptides[ ,(ncol(peptides) - quantitative.columns + 1):ncol(peptides)], 2, median, na.rm = TRUE)))
    
  } else if(normalization == "mean"){
    
    peptides[ ,(ncol(peptides) - quantitative.columns + 1):ncol(peptides)] <- data.frame(t(t(peptides[ ,(ncol(peptides) - quantitative.columns + 1):ncol(peptides)]) - apply(peptides[ ,(ncol(peptides) - quantitative.columns + 1):ncol(peptides)], 2, mean, na.rm = TRUE)))
    
  } else {
    
    stop("Error: normalization method is not valid!")
    
  }
  
  #Median normalization of the modified peptide abundances
  
  #Separate the modified peptides that were not mapped.
  unmapped.modified.peptides <- peptides[which(lengths(modified.peptides.mapping) == 0), ]
  rownames(unmapped.modified.peptides) <- 1:nrow(unmapped.modified.peptides)

  #Separate the mapped modified peptides
  mapped.modified.peptides <- peptides[-which(lengths(modified.peptides.mapping) == 0), ]
  rownames(mapped.modified.peptides) <- 1:nrow(mapped.modified.peptides)

  #Assign modified peptides to a single protein group.
  which.proteins <- vector(mode = "list", length = length(mapped.modified.peptides))

  for (i in 1:length(modified.peptides.mapping)){

    prots <- modified.peptides.mapping[[i]]
    which <- unlist(lapply(prots, function(x) which(stringi::stri_detect_fixed(proteins$Protein, pattern = x))))

    if(length(which) > 1){

      if(length(unique(which)) == 1){

        which.proteins[[i]] <- which[1]

      } else {

        #If there are multiple protein groups, select the one with the most assigned peptides.
        id <- unique(which)
        size <- proteins$Total.quantified.peptides[id]
        which.proteins[[i]] <- id[which(size == max(size))[1]]

      }

    } else {

      which.proteins[[i]] <- unlist(lapply(prots, function(x) which(stringi::stri_detect_fixed(proteins$Protein, pattern = x))))

    }

  }
  
  which.proteins <- unlist(which.proteins)
  
  #Normalize mapped modified peptide abundances as a log2 ratio to the corresponding protein expression.
  normalized.abundances <- as.data.frame(t(sapply(1:nrow(mapped.modified.peptides), function(x) as.numeric(mapped.modified.peptides[x, (ncol(mapped.modified.peptides) - quantitative.columns + 1):ncol(mapped.modified.peptides)] - proteins[which.proteins[x], (ncol(proteins) - quantitative.columns + 1):ncol(proteins)])) ))
  
  #Normalize the result by median or mean.
  if(normalization == "median"){
    
    normalized.abundances <- data.frame(t(t(normalized.abundances) - apply(normalized.abundances, 2, median, na.rm = TRUE) ))
    
  } else if(normalization == "mean"){
    
    normalized.abundances <- data.frame(t(t(normalized.abundances) - apply(normalized.abundances, 2, mean, na.rm = TRUE) ))
    
  }
  
  mapped.modified.peptides <- cbind(mapped.modified.peptides[ ,1:2], proteins$Protein[which.proteins], normalized.abundances)
  colnames(mapped.modified.peptides) <- c("Peptide", "Counterpart", "Protein group", colnames(proteins[,(ncol(proteins) - quantitative.columns + 1):ncol(proteins)]))

  
  #QC report
  #Create directory for QC plots.
  if(all(!list.files() == "QC plots")){ dir.create("QC plots")}
  
  #Count number of missing values per sample and create a barplot.
  NAs <- data.frame(Samples = colnames(mapped.modified.peptides)[-c(1:3)], value = colSums(is.na(mapped.modified.peptides[,c(4:ncol(mapped.modified.peptides))])))
  NAs$Conditions <- substring(as.character(NAs$Samples), 3, 3)
  NAs$Samples <- factor(NAs$Samples, levels = NAs$Samples)

  p.NAs <- ggplot2::ggplot(NAs, ggplot2::aes(x = Samples, y = value, fill = Conditions)) + 
    ggplot2::geom_bar(stat = "identity", position = "dodge", ggplot2::aes(color = Samples))  +
    ggplot2::geom_text(ggplot2::aes(label = value, y = value), position = ggplot2::position_dodge(1), vjust = -1) +
    ggplot2::ylab("Missing values count") +
    ggplot2::xlab("Samples") +
    ggplot2::ggtitle("Number of missing values per sample") +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::guides(color = FALSE)
  
  exprs.long <- reshape::melt(mapped.modified.peptides[ ,-c(1:3)], id.vars = NULL)
  exprs.long$Conditions <- substring(as.character(exprs.long$variable), 3,3)
  
  #Distribution plot per condition
  p.Dist <- ggplot2::ggplot(na.omit(exprs.long), ggplot2::aes(x = value, color = Conditions, ..scaled..)) + 
    ggplot2::geom_density() +
    ggplot2::ylab("Scaled density") +
    ggplot2::xlab("PTM abundance (log2)") +
    ggplot2::ggtitle(paste0("Distributions of normalized (", normalization,") PTM abundances per condition")) +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  #Violin plot combined with box plot of PTMs abundances per sample.
  p.Violin <- ggplot2::ggplot(na.omit(exprs.long), ggplot2::aes(x = variable, y = value)) +
    ggplot2::geom_violin(width = 1, ggplot2::aes(fill = Conditions)) +
    geom_boxplot(width = 0.1, ggplot2::aes(color = variable)) +
    ggplot2::ylab("PTM abundance (log2)") +
    ggplot2::xlab("Samples") +
    ggplot2::ggtitle(paste0("Violin plot of normalized (", normalization,") PTMs abundances per condition")) +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::guides(color = FALSE, fill = FALSE)
  
  
  #PCA analysis and PCA plot of the two first components.
  pca <- stats::prcomp(~ . ,mapped.modified.peptides[ ,-c(1:3)] ,na.action = na.omit)
  
  pca.df <- data.frame(pca$rotation[, 1:2], rownames(pca$rotation), substr(rownames(pca$rotation), start = 1, stop = 3), stringsAsFactors = F)
  
  colnames(pca.df) <- c("x", "y", "Samples", "Conditions")
  
  p.PCA <- ggplot2::ggplot(data = pca.df, ggplot2::aes(x = x, y = y, color = Conditions, fill = Samples)) +
    ggplot2::geom_point(size = 3, shape = 21, stroke = 2) +
    ggplot2::xlab(paste0("Component ", 1, " - ", round(pca$sdev[1]^2/(sum(pca$sdev^2)), digits = 2)*100, " %" )) +
    ggplot2::ylab(paste0("Component ", 2, " - ", round(pca$sdev[2]^2/(sum(pca$sdev^2)), digits = 2)*100, " %" )) + 
    ggplot2::ggtitle(paste0("PCA plot of per sample PTM abundances (Components ", 1, " and ", 2, ")")) +
    ggplot2::theme(plot.margin = ggplot2::margin(50, 50, 50, 50))
  
  
  #Save pdf file
  ggsave("QC plots/QC_PTMs.pdf", plot = grid.arrange(p.NAs, p.Dist, p.Violin, p.PCA, ncol = 2, nrow = 2, top = "Quality control plots for PTMs data"),
         width = 22, height = 19, limitsize = FALSE)
  
  
  unmapped.modified.peptides <- unmapped.modified.peptides[,c(1:2, (ncol(unmapped.modified.peptides) - quantitative.columns + 1):ncol(unmapped.modified.peptides) )]
  
  
  return(list(mapped = mapped.modified.peptides, unmapped = unmapped.modified.peptides))
}


