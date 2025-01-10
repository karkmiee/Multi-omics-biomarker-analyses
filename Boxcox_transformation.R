### Box-Cox transformation of metabolomics data.


#Step 1. Install packages
#BiocManager::install("MASS")
library(MASS)

#-------------------------------------------------------------------------------
##Step 2: Import and sort data. The test can't handle 0 values or NAs. 

#Here is an example how to identify columns with 0 or NA values. 
Dat1 <- read.table("mets_filtered.txt")

# Identify columns with 0 values
zero_columns <- sapply(Dat1, function(x) any(x == 0))

# Print the column names containing 0 values
print(names(zero_columns)[zero_columns])

# Check for NA values in each column
na_columns <- colSums(is.na(Dat1)) > 0

# Print column names containing NA values
print(names(na_columns)[na_columns])


#Extract columns with 0 and NAs as own df and erase observations with no values. 


#------------------------------------------------------------------------------
#We have now identified the columns with 0 or NAs and extracted them as separate dfs. 
#Make separate transformation for each df and combine them by common id.

# Import metabolomic dataframes
Dat <- read.table("mets1_filtered.txt")
Dat2 <- read.table("mets2_filtered.txt")
Dat3 <- read.table("mets3_filtered.txt")
Dat4 <- read.table("mets4_filtered.txt")
Dat5 <- read.table("mets5_filtered.txt")

#--------------------------------------------------------------------------------
#Step 3. Boxcox transformation and significance check.

Dat <- Dat4 #Repeat the steps below for each df.

#convert values from character to numeric
row_names <- rownames(Dat) # Save row names
Dat1 <- data.frame(lapply(Dat, function(x) gsub(",", ".", x)), row.names = row_names) # Replace commas with dots
Dat1 <- data.frame(lapply(Dat1, function(x) as.numeric(as.character(x))), row.names = row_names) # Convert values from character to numeric in all columns of 'Dat1'

sapply(Dat1, class)

# Data frame for transformed data
Dat_boxcox <- Dat1
#Dat_boxcox2 <- Dat1
#Dat_boxcox3 <- Dat1
#Dat_boxcox4 <- Dat1
#Dat_boxcox5 <- Dat1


# Data frame for lambda parameters and diagnostics (Shapiro-Wilk test)
lambda_df <- data.frame(metabolite = character(0), # metabolite label
                        lambda = numeric(0), # lambda parameter
                        SW_test_raw = numeric(0), # Shapiro-Wilk test for raw data
                        SW_test_boxcox = numeric(0)) # Shapiro-Wilk test for transformed data

# Metabolite labels (change this based on what df you use and also to command below)
Dat_cols <- colnames(Dat_boxcox)

# Box Cox transformation and diagnostics - works only for metabolite data that has no zero values or NAs
for (i in c(1:164)) { # insert number of metabolite columns
  b <- boxcox(lm(Dat1[[i]] ~ 1)) # identify lambda
  lambda <- b$x[which.max(b$y)] # extract lambda
  lambda_df[i, "metabolite"] <- Dat_cols[[i]] # metabolite name
  lambda_df[i, "lambda"] <- lambda # lambda
  
  Dat_boxcox[, i] <- (Dat_boxcox[[i]] ^ lambda_df[i, "lambda"] - 1) / lambda_df[i, "lambda"] # transform data
  
  lambda_df[i, "SW_test_raw"] <- shapiro.test(Dat1[[i]])$p.value # Shapiro-Wilk test for non-transformed data
  lambda_df[i, "SW_test_boxcox"] <- shapiro.test(Dat_boxcox[[i]])$p.value # Shapiro-Wilk test for transformed data
}


##------------------------------------------------------------------------------
#Step 4. Combine data

# Extract row names
filenames <- rownames(Dat_boxcox) #extract row names as vector
#filenames <- rownames(Dat_boxcox2)
#filenames <- rownames(Dat_boxcox3)
#filenames <- rownames(Dat_boxcox4)
#filenames <- rownames(Dat_boxcox5)

# Add row names as a new column
df <- cbind(filenames, Dat_boxcox)
df2 <- cbind(filenames, Dat_boxcox2)
df3 <- cbind(filenames, Dat_boxcox3)
df4 <- cbind(filenames, Dat_boxcox4)
df5 <- cbind(filenames, Dat_boxcox5)

# Perform the left join
combined_df <- merge(df, df2, by = "filenames", all.x = TRUE)
combined_df <- merge(combined_df, df3, by = "filenames", all.x = TRUE)
combined_df <- merge(combined_df, df4, by = "filenames", all.x = TRUE)
combined_df <- merge(combined_df, df5, by = "filenames", all.x = TRUE)

sapply(combined_df, class)

# Set the extracted filenames as row names
filenames <- combined_df$filenames
rownames(combined_df) <- filenames
combined_df <- combined_df[,-c(1)] #delete first column

# Save transformed values
write.table(combined_df, file="boxcox_transf_data.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(lambda_df, file="boxcox_transf_data_SW_test.txt", sep="\t", row.names = TRUE, col.names = TRUE)


