library(limma)
library(edgeR)

# Reading in the feature count file as "counts.df"
counts.df <- read.csv("featCounts_S_cere_20200331.csv")

# Printing the start of the counts.df object in R...
head(counts.df)

# Using the "Geneid" column to set the rownames
rownames(counts.df) <- counts.df$Geneid

# Removing the "Geneid" column
counts.df$Geneid <- NULL

# Printing the start of the counts.df object in R...
head(counts.df)

# Reading in the design table as "design.df"
design.df <- read.csv("design_table.csv", fileEncoding="UTF-8-BOM")

# Printing the start of the design.df object in R...
print(design.df)

# Subsetting gene counts according to experimental condition
counts_standard.df  <- counts.df[,c("SRR8933535", "SRR8933536", "SRR8933537")]
counts_anaerobic.df <- counts.df[,c("SRR8933506", "SRR8933511", "SRR8933512")]
counts_high_temp.df <- counts.df[,c("SRR8933532", "SRR8933533", "SRR8933534")]
counts_low_pH.df    <- counts.df[,c("SRR8933530", "SRR8933531", "SRR8933539")]
counts_pressure.df  <- counts.df[,c("SRR8933509", "SRR8933510", "SRR8933538")]

# Printing the structure of the gene counts set and subsets
str(counts.df)

str(counts_standard.df)
str(counts_anaerobic.df)
str(counts_high_temp.df)
str(counts_low_pH.df)
str(counts_pressure.df)

