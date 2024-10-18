#install.packages("rstudioapi")
library(rstudioapi)
#install.packages("readxl")
library(readxl)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("reshape")
library(reshape)
#install.packages("reshape2")
library(reshape2)
#install.packages("grid")
library(grid)
#install.packages("gridExtra")
library(gridExtra)
#install.packages("scales")
library(scales)

#Functions of meRgreeMS

#Function to merge PSM datasets from different MS runs.
datasets.Merge <- function(files, runIDs, pepC, ptmC, intensitiesC, referenceC, log2, statC,statisticalThreshold, localC,localizationThreshold, additionalC, label, mods, colnames, norm = "none"){
  
  #IN
  #files: PSM data in xlsx format (PD reports)
  #runIDs: run identifiers
  #pepC: index corresponding to peptide sequence column in the imported files
  #ptmC: index corresponding to modification information column in the imported files
  #intensitiesC: columns corresponding to peptide intensities in the imported files
  #referenceC: the sample index in intensitiesC that is a reference channel (used for rescaling)
  #log2: TRUE if data are in log2 scale, else FALSE 
  #statC: index of the identification score column in the imported files
  #statisticalThreshold: identification score threshold
  #localC: index of the PTM localization score column in the imported files
  #localizationThreshold: PTM localization score threshold
  #additionalC: column index of a string column in the imported files to filter additional tag attributes
  #label: tag value to be excluded 
  #mods: array of modification names
  #colnames: array of names that channels should get once the merging is done (based on experimental design)
  #norm: the type of normalization applied to data from different runs before they get merged (default: none, options: median, quantile)
  #OUT
  #merged_runs: data frame of the merged PSM reports
  
  #Create directory for QC plots.
  if(all(!list.files() == "QC plots")){ dir.create("QC plots")}
  
  
  #Preprocessing of every dataset. (Filtering, log2-transformation, normalization, rescale and QC per file)
  sets <- vector(mode = "list", length = length(files))

  for (i in 1:length(files)) {
    
    sets[[i]] <- dataset.Preprocess(file = files[i], pepC = pepC, ptmC = ptmC, intensitiesC = intensitiesC, 
                                    referenceC = referenceC, log2 = log2, statC = statC, statisticalThreshold =  statisticalThreshold, 
                                    localC = localC, localizationThreshold = localizationThreshold, additionalC = additionalC, 
                                    label = label, mods = mods)
  }
  
  #Aggerate files according to the run they belong.
  runs <- vector(mode = "list", length = length(runIDs))
  
  for (i in 1:length(runIDs)) {
    
    which.files <- which(grepl(runIDs[i], files))
    
    aggregate.set <- sets[[which.files[1]]]
    
    if(length(which.files) > 1){
      
      for (j in 2:length(which.files)) {
        
        aggregate.set <- rbind(aggregate.set, sets[[which.files[j]]])
        
      }
      
    }
    
    aggregate.set <- aggregate(x = aggregate.set[,2:ncol(aggregate.set)], by = list(aggregate.set[, 1]), FUN = mean, na.rm = T)
    
    if(norm == "median"){
      
      aggregate.set[,2:ncol(aggregate.set)] <- data.frame(t(t(aggregate.set[,2:ncol(aggregate.set)]) - apply(aggregate.set[,2:ncol(aggregate.set)], 2, median, na.rm = TRUE)))
      
    } else if(norm == "quantile"){
      
      aggregate.set[,2:ncol(aggregate.set)] <- preprocessCore::normalize.quantiles(as.matrix(aggregate.set[,2:ncol(aggregate.set)]))
      
      aggregate.set[,2:ncol(aggregate.set)] <- scale(aggregate.set[,2:ncol(aggregate.set)], scale = FALSE)
      
    }
      
    runs[[i]] <- aggregate.set
    
  }
  
  
  #Unique peptides in all runs.
  unique.peptides <- sort(unique(unlist(lapply(1:length(runIDs), function(x) runs[[x]][,1]))))
  
  #Merge the data from all runs.
  merged_runs <- data.frame(matrix(NA, nrow = length(unique.peptides), ncol = (sum(unlist(lapply(1:length(runIDs), function(x) ncol(runs[[x]]) - 1))) + 1)))
  merged_runs[ ,1] <- unique.peptides
  indices <- sapply(runs, function(x) sapply(unique.peptides, function(y) if(length(which(x[,1] == y)) == 1){return(which(x[,1] == y))} else {return(NA)}))

  for (i in 1:length(runIDs)) {
     
     starting.column <- if(i != 1){ sum(unlist(lapply(1:(i-1), function(x) ncol(runs[[x]]) - 1))) + 2 } else {2}
     merged_runs[, starting.column:(starting.column + ncol(runs[[i]]) - 2)] <- runs[[i]][indices[,i], 2:ncol(runs[[i]])]

  }
  
  colnames(merged_runs) <- colnames
  merged_runs <- merged_runs[ ,c(1, order(colnames(merged_runs))[1:(ncol(merged_runs)-1)])]
  
  #Save pdf file with all QC plots.
  ggsave("QC plots/QC_merged.pdf", plot = grid.arrange(peptides.Bar.Plot(runs = runIDs, files = files, data = c(unlist(lapply(sets, function(x) nrow(x))), unlist(lapply(runs, function(x) nrow(x))), nrow(merged_runs))),
                                                       NA.Bar.Plot(data = merged_runs[,2:ncol(merged_runs)]),
                                                       distribution.Plot(file = "Merged runs", label = "Peptide", data = merged_runs[,2:ncol(merged_runs)]),
                                                       box.Plot(file = "Merged runs", label = "Peptide", data = merged_runs[,2:ncol(merged_runs)]),
                                                       add.Density.Lines(p = ggplot2::ggplot() + ggplot2::ylab("Density") + ggplot2::xlab("Peptide intensity") + ggplot2::ggtitle("Density plot per sample") + ggplot2::theme(plot.margin = ggplot2::margin(50, 50, 50, 50)), data = na.omit(reshape::melt(merged_runs[,2:ncol(merged_runs)], id.vars = NULL)), merged = T),
                                                       cv.Box.Plot(data = merged_runs[,2:ncol(merged_runs)]),
                                                       pca.Plot(data = merged_runs[,2:ncol(merged_runs)], pc1 = 1, pc2 = 2),
                                                       pca.Plot(data = merged_runs[,2:ncol(merged_runs)], pc1 = 2, pc2 = 3),
                                                       cor.Heatmap(data = merged_runs),
                                                       ncol = 2, nrow = 5, top = "Quality control plots for merged run data"),
                                                       width = 22, height = 50, limitsize = FALSE)
          
  return(merged_runs)
  
}

#Function to filter, normalize and extract peptide sequences.
dataset.Preprocess <- function(file, pepC, ptmC, intensitiesC, referenceC, log2, statC,statisticalThreshold, localC,localizationThreshold, additionalC, label, mods){
  
  dataset <- as.data.frame(readxl::read_xlsx(file, guess_max = 50000))
  
  #A) Filtering.
  
  #A.1) Filter out PSMs with localization score (ptmRS) lower than localizationThreshold to at least one site.
  
  localization.scores <- lapply(dataset[,localC], function(x) as.numeric(paste0(unlist(regmatches(x, gregexpr("[: ][[:digit:]]+", x))))))
  
  localization.index <- unlist(lapply(localization.scores, function(x){
    
    if(length(x) == 0){ return(TRUE) }
    
    if(all(x > localizationThreshold)){ return(TRUE) } else { return(FALSE) }
    
  }))
  
  #Remove lozalization.index modified PSMs.
  if(!all(localization.index)){
    
    data_set <- dataset[localization.index,]
    rownames(data_set) <- 1:nrow(data_set)
    
  } else {
    
    data_set <- dataset
    
  }
  
  #A.2) Filter out PSMs under an additional filtering criteria (a label).
  if(!missing(additionalC)){
    
    if(length(which(data_set[, additionalC] == label)) > 0) {
      
      data_set <- data_set[-which(data_set[, additionalC] == label), ]
      rownames(data_set) <- 1:nrow(data_set)
      
    }
      
  }
  
  #A.3) Remove PSMs with statistical score lower than statisticalThreshold and keep only columns with pepC, ptmC and intesitiesC columns.
  data_set <- data_set[(data_set[ ,statC] < statisticalThreshold), c(pepC, ptmC, intensitiesC)]
  rownames(data_set) <- 1:nrow(data_set)

  #A.4) Remove PSMs that referenceC is missing or all intensities are missing.
  na.rows <- sort(unique(c(which(apply(is.na(data_set[ ,setdiff(3:ncol(data_set), referenceC + 2)]), 1, all)), which(is.na(data_set[, referenceC + 2])))))
  if(length(na.rows) > 0){
    
    data_set <- data_set[-na.rows,]
    rownames(data_set) <- 1:nrow(data_set)
    
  }
  
  #Bar plot - amount of remaining and filtered PSMs
  PSMs.Bar.Plot <- PSMs.Bar.Plot(file = file, data =  c(length(which((dataset[ ,statC] >= statisticalThreshold))), 
                                                                length(which(!localization.index)),
                                                                length(which(dataset[, additionalC] == label)),
                                                                length(unique(c(which(apply(is.na(dataset[setdiff(intensitiesC, (intensitiesC[1] + referenceC - 1))]), 1, all)), which(is.na(dataset[, (intensitiesC[1] + referenceC - 1)]))))),
                                                                nrow(dataset) - nrow(data_set),
                                                                nrow(data_set)) )
  
  #Distribution plot of filtered PSM intensities.
  PSMs.Distribution.Plot <- distribution.Plot(file = file, label = "PSM" ,data = data_set[,3:ncol(data_set)])
  
  #B) Turn to log2 scale if is not in log2 scale already.
  if(!log2){
    
    data_set[,3:ncol(data_set)] <- log2(data_set[,3:ncol(data_set)])
  
  }

  #Distribution plot of filtered PSM intensities in log2 scale.
  PSMs.Distribution.Plot.Log2 <- distribution.Plot(file = file, label = "Log2 PSM"  , data = data_set[,3:ncol(data_set)])

  #Box plot of PSM intensities in log2 scale for each sample.
  PSMs.Box.Plot.Log2 <- box.Plot(file = file, label = "Log2 PSM", data = data_set[,3:ncol(data_set)])
  
  #C) Normalize to median.
  data_set[,3:ncol(data_set)] <- data.frame(t(t(data_set[,3:ncol(data_set)]) - apply(data_set[,3:ncol(data_set)], 2, median, na.rm = TRUE)))
  
  #Distribution plot of log2 normalized by median PSM intensities.
  PSMs.Distribution.Plot.Normalized <- distribution.Plot(file = file, label = "Normalized (Median) Log2 PSM", data = data_set[,3:ncol(data_set)])
  
  #Box plot of normalized PSM intensities in log2 scale for each sample.
  PSMs.Box.Plot.Normalized <- box.Plot(file = file, label = "Normalized (Median) Log2 PSM", data = data_set[,3:ncol(data_set)])

  #D) Rescale to internal standard referenceC.
  data_set[,3:ncol(data_set)] <- data_set[,3:ncol(data_set)] - data_set[,(referenceC + 2)]
  data_set <- data_set[,-(referenceC + 2)]

  #Distribution plot of intensities rescaled to internal reference.
  PSMs.Distribution.Plot.Scaled <- distribution.Plot(file = file, label = "Scaled Log2 PSM", data = data_set[,3:ncol(data_set)])

  #Box plot of normalized PSM intensities in log2 scale for each sample.
  PSMs.Box.Plot.Scaled <- box.Plot(file = file, label = "Scaled Log2 PSM", data = data_set[,3:ncol(data_set)])
  
  #E) Split dataset to modified and unmodified partitions.
  #Search for all unique type of modifications in data.
  unique_mods <- sort(unique(unlist(lapply(1:nrow(data_set), function(x) regmatches(data_set[x,2], gregexpr("(?<=\\().*?(?=\\))", data_set[x,2], perl = T))[[1]]))))
  #Intersect that to argument modifications.
  unique_mods <- intersect(unique_mods, mods)
  mod_IDs <- sort(unlist(lapply(unique_mods, function(x) which(grepl(x, data_set[,2])))))
  
  #Modified PSMs.
  data_set_mods <- data_set[mod_IDs,]
  rownames(data_set_mods) <- 1:nrow(data_set_mods)
  
  #Unmodified PSMs.
  data_set <- data_set[-mod_IDs,]
  rownames(data_set) <- 1:nrow(data_set)
  
  #For unmodified PSMs (extract PSM sequences and aggregate on peptide level)
  unmod_sequences <- toupper(gsub(".*\\.(.*)\\..*", "\\1", data_set[,1]))
  data_set[,1] <- unmod_sequences
  data_set <- data_set[,-2]
  data_set <- aggregate(x = data_set[,2:ncol(data_set)], by = list(data_set[, 1]), FUN = mean, na.rm = T)
  
  #For modified PSMs (extract PSM sequences annotated with ptms and aggregate on peptide level)
  data_set_mods[,1] <- unlist(lapply(1:nrow(data_set_mods) , function(x) map.Modifications(mods = data_set_mods[x,2], sequence = data_set_mods[x,1], unique.mods = unique_mods)))
  data_set_mods <- data_set_mods[,-2]
  data_set_mods <- aggregate(x = data_set_mods[,2:ncol(data_set_mods)], by = list(data_set_mods[, 1]), FUN = mean, na.rm = T)

  #Merge the modified and unmodified peptides by row.
  merged_set <- rbind(data_set, data_set_mods)
  
  #Distribution plot of peptide intensities intensities, with additional distributions of modified and unmodified peptides.
  Peptide.Distribution.Plot <- distribution.Plot(file = file, label = "Peptide", data = merged_set[,2:ncol(merged_set)])
  Peptide.Distribution.Plot <- add.Density.Lines(p = Peptide.Distribution.Plot, 
                                                 data = data.frame(value = c(c(t(data_set[,2:ncol(data_set)])), c(t(data_set_mods[,2:ncol(data_set_mods)]))),
                                                                   variable = c(rep("Unmodified", nrow(data_set)), c(rep("Modified", nrow(data_set_mods)))) ))
  
  #Box plot of peptide intensities for each sample.
  Peptide.Box.Plot <- box.Plot(file = file, label = "Peptide", data = merged_set[,2:ncol(merged_set)])
  
  ggsave(filename = paste0("QC plots/", file, ".pdf") , plot = gridExtra::grid.arrange(grid::textGrob("Filtered PSMs", rot = 90), PSMs.Bar.Plot, PSMs.Distribution.Plot,
                                                                                       grid::textGrob("Log2 PSMs", rot = 90), PSMs.Distribution.Plot.Log2, PSMs.Box.Plot.Log2,
                                                                                       grid::textGrob("Normalized PSMs", rot = 90), PSMs.Distribution.Plot.Normalized, PSMs.Box.Plot.Normalized,
                                                                                       grid::textGrob("Rescaled PSMs", rot = 90), PSMs.Distribution.Plot.Scaled, PSMs.Box.Plot.Scaled,
                                                                                       grid::textGrob("Peptides", rot = 90), Peptide.Distribution.Plot, Peptide.Box.Plot,
                                                                                       widths = c(0.1, 0.45, 0.45),
                                                                                       ncol = 3, nrow = 5, top = paste0("QC - ", file)), width = 22, height = 50, limitsize = FALSE)
  return(merged_set)

}

#Function that extract annotated peptide sequences (in VIQoR format)
map.Modifications <- function(mods, sequence, unique.mods){
  
  mods <- unlist(strsplit(mods, "; "))
  sequence <- toupper(gsub(".*\\.(.*)\\..*", "\\1", sequence))
  mods <- mods[which(grepl(paste0("(", paste(unique.mods, collapse = ")|(", sep = ""), ")", sep = "", collapse = ""), mods))]
  
  index <- rep(NA, length(mods))
  type <- rep(NA, length(mods))
  for(i in 1:length(mods)){
    
    index[i] <- as.numeric(paste0(unlist(regmatches(mods[i], gregexpr("[[:digit:]]+", mods[i])))))
    type[i] <- regmatches(mods[i], gregexpr("(?<=\\().*?(?=\\))", mods[i], perl = T))[[1]]
    
  }
  
  extra.chars <- 0
  
  for (i in 1:length(index)) {
    
    expression <- paste0('^(.{', index[i]-1+extra.chars, '})(.*)$')
    add <- paste0('\\1(', type[i], ')\\2')
    extra.chars <- extra.chars + nchar(type[i]) + 2
    sequence <- gsub(expression, add, sequence)
    
  }
  
  return(sequence)
  
}

#Function to generate a PSMs filtering bar plot.
PSMs.Bar.Plot <- function(file, data){

  data <- data.frame(variable = c("Removed by q-value", "Removed by ptmRS", "Removed by Quan", "Missing values","Total removed", "Total remained"),
                     value = data, stringsAsFactors = F)
  data$variable <- factor(data$variable, levels = data$variable)
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = variable, y = value, fill = variable)) + 
                ggplot2::geom_bar(stat = "identity", position = "dodge") +
                ggplot2::geom_text(ggplot2::aes(label = value, y = value), position = ggplot2::position_dodge(1), vjust = -1) +
                ggplot2::ylab("PSMs count") + 
                ggplot2::xlab("Category") +
                ggplot2::ggtitle(label = paste0(file, " - PSMs filtering")) +
                ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(50, 50, 50, 50))

  return(p)  
  
}

#Function to generate a bar plot of unique peptides summarized per file/run/total.
peptides.Bar.Plot <- function(runs, files, data){
  
  files.per.run <- lapply(runs, function(x) which(grepl(x, files)))
  group <- c(unlist(lapply(1:length(runs), function(x) rep(runs[x], lengths(files.per.run)[x])))[order(unlist(files.per.run))], rep("Runs", length(runs)), "Merged")

  data <- data.frame(variable = c(files, runs, "Total"), value = data, group = group)
  data$variable <- factor(data$variable, levels = data$variable)
  data$group <- factor(data$group, levels = unique(data$group))

  p <- ggplot2::ggplot(data, ggplot2::aes(x = variable, y = value, fill = variable)) +
                ggplot2::geom_bar(stat = "identity", position = "dodge") +
                ggplot2::geom_text(ggplot2::aes(label = value, y = value), position = ggplot2::position_dodge(1), vjust = -1) +
                ggplot2::facet_grid(~group, switch = "x", scales = "free_x", space = "free_x") +
                ggplot2::ylab("Peptide count") +
                ggplot2::xlab("Labels") + 
                ggplot2::ggtitle("Unique peptide count for all files and runs") +
                ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(50, 50, 50, 50))
  
  maxChar <- max(sapply(as.character(data$variable), function(x) nchar(x)))
  
  if(maxChar > 12){
    
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    
  }
  
  return(p)
}

#Creates a bar plot for the amount of missing values to each sample.
NA.Bar.Plot <- function(data){

  data <- data.frame(variable = colnames(data), value = colSums(is.na(data)))
  data$variable <- factor(data$variable, levels = data$variable)

  p <- ggplot2::ggplot(data, ggplot2::aes(x = variable, y = value, fill = variable)) +
                ggplot2::geom_bar(stat = "identity", position = "dodge") +
                ggplot2::geom_text(ggplot2::aes(label = value, y = value), position = ggplot2::position_dodge(1), vjust = -1) +
                ggplot2::ylab("Missing values count") +
                ggplot2::xlab("Samples") +
                ggplot2::ggtitle("Amount of missing values per sample") +
                ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(50, 50 ,50 ,50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  return(p)
}

#Function to generate distribution plots for whole PSMs/peptides datasets.
distribution.Plot <- function(file, label, data){
  
  binwidth <- (max(data, na.rm = T) - min(data, na.rm = T))/100
  data <- reshape::melt(data, id.vars = NULL)
  
  p <- ggplot2::ggplot(data = na.omit(data), mapping = ggplot2::aes(x = value)) + 
                ggplot2::geom_histogram(ggplot2::aes(y = ..density..), binwidth = binwidth, colour = "black", fill = "white") + 
                ggplot2::geom_density(alpha = 0.1, fill = "blue", colour = "blue") +
                ggplot2::xlab(paste0(label, " intensities")) +
                ggplot2::ylab("Density") +
                ggplot2::ggtitle(paste0(file, " - ", label, " intensities distribution plot.")) +
                ggplot2::theme(plot.margin = ggplot2::margin(50, 50, 50, 50))
    
  return(p)
  
}

#Function to add additional densitry lines to a ggplot graph.
add.Density.Lines <- function(p, data, merged = FALSE){

  p <- p + ggplot2::geom_density(data = na.omit(data), ggplot2::aes(x = value, colour = variable), alpha = 0.4, show.legend = FALSE) + 
                    ggplot2::stat_density(data = na.omit(data), ggplot2::aes(x = value, colour = variable), geom = "line", position = "identity") +
                    ggplot2::theme(legend.position = c(0.1, 0.7), legend.background = ggplot2::element_rect(fill = scales::alpha(colour = "white", alpha = 0))) 
  
  if(!merged){
    
    p <- p + ggplot2::scale_colour_manual(name = "", values = c("green", "red"))
    
  } else {
    
    p <- p + ggplot2::labs(color = "Samples")
    
  }
  
  return(p)
  
}

#Function to create box plot for each sample.
box.Plot <- function(file, label, data){
  
  data <- reshape::melt(data, id.vars = NULL)
  
  p <- ggplot2::ggplot(data = na.omit(data), mapping = ggplot2::aes(x = variable, y = value)) + 
                ggplot2::geom_boxplot(ggplot2::aes(color = variable)) +
                ggplot2::xlab("Samples") +
                ggplot2::ylab(paste0(label, " intensities")) +
                ggplot2::ggtitle(paste0(file, " - ", label, " intensities boxplot.")) +
                ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(50, 50, 50, 50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  return(p)
  
}

#Function to create PCA plot for given components.
pca.Plot <- function(data, pc1, pc2){
  
  pca <- prcomp(~ ., data, na.action = na.omit)

  df <- data.frame(pca$rotation[, c(pc1, pc2)], rownames(pca$rotation), substr(rownames(pca$rotation), start = 1, stop = 3), stringsAsFactors = F)

  colnames(df) <- c("x", "y", "sample", "condition")
  
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y, color = condition, fill = sample)) +
                ggplot2::geom_point(size = 3, shape = 21, stroke = 2) +
                ggplot2::xlab(paste0("Component ", pc1, " - ", round(pca$sdev[pc1]^2/(sum(pca$sdev^2)), digits = 2)*100, " %" )) +
                ggplot2::ylab(paste0("Component ", pc2, " - ", round(pca$sdev[pc2]^2/(sum(pca$sdev^2)), digits = 2)*100, " %" )) + 
                ggplot2::ggtitle(paste0("PCA plot of peptide relative abundances (Components ", pc1, " and ", pc2, ")")) +
                ggplot2::theme(plot.margin = ggplot2::margin(50, 50, 50, 50))
    
  return(p)  
}

#Function to create sample-wise correlation heatmap.
cor.Heatmap <- function(data){
  
  data <- round(cor(data[,2:ncol(data)], use = "complete.obs"), digits = 2)
  data[lower.tri(data, diag = F)] <- NA
  data <- reshape::melt(data)
  colnames(data) <- c("Var1", "Var2", "value")
  
  p <- ggplot2::ggplot(na.omit(data), ggplot2::aes(x = Var2, y = Var1, fill = value)) +
                ggplot2::geom_tile(color = "white") +
                ggplot2::geom_text(ggplot2::aes(Var2, Var1, label = value)) +
                ggplot2::scale_fill_gradient2(low = "purple", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), name="Correlation") +
                ggplot2::theme_minimal() +
                ggplot2::ylab("Samples") +
                ggplot2::xlab("Samples") +
                ggplot2::ggtitle("Correlation heatmap of all samples") +
                ggplot2::theme(plot.margin = margin(50, 50 ,50 ,50), 
                          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), 
                          panel.grid.major = ggplot2::element_blank(), 
                          panel.background = ggplot2::element_blank(), 
                          legend.position = c(0.1, 0.9))
  
  return(p)
  
}

#Function to create per condition CV boxplot of peptide features.
cv.Box.Plot <- function(data){
  
  data <- 2^data
  
  conditions <- unique(substring(colnames(data), 1, 3))
  conditions.B <- lapply(conditions, function(x) grepl(pattern = x, colnames(data)))
  conditions.CV <- lapply(conditions.B, function(x) 100 * apply(data[,x], 1, sd, na.rm = TRUE)/rowMeans(data[,x], na.rm = T))
  conditions.CV <- as.data.frame(do.call(cbind, conditions.CV))
  colnames(conditions.CV) <- conditions
  
  #Data structure to plot.
  data <- reshape::melt(conditions.CV, id.vars = NULL)
  
  #Annotation table
  annotation <- rbind(round(colMeans(conditions.CV, na.rm = T), 2), round(apply(conditions.CV, 2, median, na.rm = T), 2))
  annotation <- cbind(c("Mean CV", "Median CV"), annotation)
  colnames(annotation)[1] <- ""
  
  p <- ggplot2::ggplot(data = na.omit(data), mapping = ggplot2::aes(x = variable, y = value)) + 
                ggplot2::geom_boxplot(ggplot2::aes(color = variable)) +
                ggplot2::xlab("Conditions") +
                ggplot2::ylab("Peptide CV (%)") +
                ggplot2::ggtitle(paste0("Per condition CV boxplot (total mean: ", round(mean(data$value, na.rm = T), 2), "%, total median: ", round(median(data$value, na.rm = T), 2), ")")) +
                ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(50, 50, 50, 50), axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
                ggplot2::ylim(0, max(data$value, na.rm = T) + 20) +
                ggplot2::annotation_custom(gridExtra::tableGrob(annotation, theme = ttheme_default(base_size = 10), rows = NULL), xmin = -Inf, xmax = Inf, ymin = max(data$value, na.rm = T))

  return(p)
  
}


