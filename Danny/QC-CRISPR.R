# title: "QC of Danny's CRISPR library"
# subtitle: "Library from Joung, J., et al. (2016) Nature protocols  DOI: 10.1038/nprot.2017.016"

# Load libraries
library(readxl)
library(ggplot2)
library(skimr)

# Set working directory
setwd("/home/martinezvbs/Downloads/Danny")

# Load CRISPR library
CRISPR_library <- read_excel("tf_orf_library.xlsx")

# Change Full sequence column name to Full_sequence
names(CRISPR_library)[names(CRISPR_library) == "Full sequence"] <- "Full_sequence"

# Summary
CRISPR_before <- skim(CRISPR_library)

# Bar plot of the average sizes for the sequences
CRISPR_library$Length <- nchar(CRISPR_library$Length)
ggplot(CRISPR_library, aes(x = Length)) +
  geom_bar(stat = "count", fill = "blue") +
  labs(title = "Average sizes of the sequences",
       x = "Size",
       y = "Count")

# Save table as .csv file
write.csv(CRISPR_library, "CRISPR_library.csv", row.names = FALSE)