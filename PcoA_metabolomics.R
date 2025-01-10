#title: "PCoA clustering and PERMANOVA distance analysis of metabolomics"

#Data: NMR serum metabolomics content 

#Abbreviations used in data and scripts:
#LS = lynch syndrome group, MMR genes mutation carriers
#MMR genes: MLH1, MSH2, MSH6 and PMS2


###-----------------------------------------------------------------------------

#Step 1. Load libraries, metabolomics data and phenodata

library(readxl)
library(dplyr)
library(ape)
library(vegan)
library(hagis)
library(ggplot2)
library(stats)


# Import mets data
mets1 <- read.delim("boxcox_transf_data.txt")

# Import phenodata
targets<- read.delim("phenoData_filtered.txt")


# Select phenodata which you wan to use in grouping, take out all else
targets_sel = targets[-c(77),-c(2:11, 13:34)] #for variant comparison
targets_sel = targets[,-c(2:32)] #for future cancer

#Delete PMS2 variant from mets
#rowname_to_delete <- "LSME0131"

# Subset the dataframe to exclude the row with the specified rowname
#mets2 <- mets1[rownames(mets1) != rowname_to_delete, ] #variant
mets2<-mets1 #future cancer

##------------------------------------------------------------------------------

#Step 2.Principal coordinate analysis PCoA to mets data

#Merge and reorder df
mets3 <- cbind(Filename = rownames(mets2), mets2)
mets4 <- merge(targets_sel, mets3, by = "Filename", all = TRUE)
rownames(mets4) <- mets4$Filename # Restore row names

#custom_order <- c("MLH1", "MSH2", "MSH6") # Define the custom sorting order
#sorted_df <- mets4[order(factor(mets4$Variant, levels = custom_order)), ] # Reorder the rows based on the custom order

custom_order <- c("NO", "YES")
sorted_df <- mets4[order(factor(mets4$Cancer_during_surveillance, levels = custom_order)), ] # Reorder the rows based on the custom order

mets <- sorted_df[,-c(1,2,3)]

#Compare also each group
#mets5 <- mets[-c(99:115), ] #MLH1-MSH2
#mets5 <- mets[-c(83:98), ] #MLH1-MSH6
#mets5 <- mets[-c(1:82), ] #MSH2-MSH6

# This time we need to compute the distances between the products
distance_matrix <- dist(mets, method = "euclidean")
distance_matrix

# fit the principal coordinates analysis using cmdscale
my_pcoa <- stats:::cmdscale(distance_matrix)

# print the result
print(my_pcoa)

# make file from Pcoa cordinates
write.table(my_pcoa, file="pcoa_coordinates_mets", sep="\t", quote=F, col.names=NA)

# Doing the mapping
plot(my_pcoa[,1], my_pcoa[,2])
text(my_pcoa[,1], my_pcoa[,2], labels = row.names(mets), cex= 1)
dev.off()


##--------------------------------------------------------------------------

#Step 4.PCOA visualization


# Load packages
library(ellipse)
library(tibble)


# Import PcoA_type file
filexx<- read.table("pcoa_coordinates_mets")
filexxx <- tibble::rownames_to_column(filexx, "Filename")
filex  <- merge(filexxx, targets_sel, by ="Filename", all=TRUE) # merge by row names
filex$Cancer_during_surveillance <- as.character(filex$Cancer_during_surveillance)

#rename cancer incidences from NO-YES to Healthy-Cancer
filex <- filex %>%
  mutate(Cancer_during_surveillance = ifelse(Cancer_during_surveillance == "NO", "Healthy", "Future Cancer"))
#Rename column name to status
filex <- filex %>%
  rename(Status = Cancer_during_surveillance)

p <- ggplot(data = filex, aes(x = V1, y = V2, colour =Status)) +
  geom_point() +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(face = "bold", size = 10),
    legend.key.size = unit(1, 'lines')
  ) +
  stat_ellipse() +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_gradient(low = "blue", high = "red") +
  guides(colour = "none")  # Hides the colour legend
p + facet_wrap(~Status)+
  ggtitle("Euclidean Distances PCoA cMets by Status")

##------------------------------------------------------------------------------

#Step 5. Beta-dispersion tests if the dispersion, variance, of two or more groups are significantly different or not.

# Count the number of observations for each variant
variant_counts <- count(sorted_df, Variant)

# Print the counts
print(variant_counts)

# set up groups, check order and numbers from step 2
groups <- factor(c(rep("MLH1", 82), rep("MSH2", 16), rep("MSH6", 17))) #tähäm montako missäkin ryhmässä on
groups <- factor(c(rep("MLH1", 82), rep("MSH2", 16)))
groups <- factor(c(rep("MLH1", 82), rep("MSH6", 17)))
groups <- factor(c(rep("MSH2", 16), rep("MSH6", 17)))
groups <- factor(c(rep("Healthy", 99), rep("Cancer", 17)))
                 
# calculates the beta-dispersion for each group, when comparing 2 or more
pathotype.disp <- betadisper(distance_matrix, groups)

# plot showing the dispersion for each group
plot(pathotype.disp, hull = FALSE, ellipse = TRUE)


#PERMANOVA

pathotype.adonis <- adonis2(distance_matrix ~ groups)
pathotype.adonis


#ANOSIM 

pathotype.anosim <- anosim(distance_matrix, groups)
pathotype.anosim
