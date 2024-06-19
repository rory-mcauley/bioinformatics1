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