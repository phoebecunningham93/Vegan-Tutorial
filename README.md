### Vegan Diversity Workshop"

# Setting up in R

Open an R script by choosing new file R script, and save this script "BIO771P_VeganExercise".

In this document, everything in the grey "chunks" is code that you will need to copy over into your R script and run. There is a descriptions of what each piece of code is doing, but for now you just need to use the code to understand the data and the factors influencing microbial diversity in ants. 

R has lots of functions built in, in what is known as base R, but for this kind of analysis you will need some extra "packages" which will provide some more complex functions to digest and visualise the data for diversity analysis. Most importantly you will be using the package vegan and ggplot2. 

Let's load the packages using the code chunk below. Copy and paste this into your R script, highlight the code and click the "run" button on the top right hand part of the screen. *Continue to do this for each line of code below, being careful to read the text above each chunk.* Remember you can annotate you code in the R script using #.

```{r load libraries, results='hide', message=FALSE}
# List of libraries to load
libraries <- c("vegan", "dplyr", "ggplot2","ape", "lme4")

# Function to install and load a package
load_or_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("Package", pkg, "is not installed. Attempting to install..."))
    tryCatch({
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }, error = function(e) {
      message(paste("Failed to install", pkg, "using install.packages(). Trying remotes or BiocManager..."))
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      if (!require(pkg, character.only = TRUE)) {
        if (pkg %in% rownames(available.packages())) {
          remotes::install_cran(pkg)
        } else {
          BiocManager::install(pkg)
        }
      }
      library(pkg, character.only = TRUE)
    })
  } else {
    message(paste("Package", pkg, "is already installed."))
  }
}

# Loop through the list and load or install each package
for (lib in libraries) {
  load_or_install(lib)
}

```

## Bringing in the data

Now we will bring in the datasets, these have been saved as either comma separated files (.csv) or tab separated file (.txt) and need to be brought into the R workspace using either read.csv() or read.table(). We need two datasets for our analysis:

- Metadata which contains all the information we now about each individual sample
- Presence/absence data of the microbes in the form of a count matrix of *amplicon sequence variants (ASV's)* 

``` {r data}
meta_data <- read.csv("/Users/phoebecunningham/Desktop/QMUL/Presentations/ant-metadata.csv")
asv_table <- read.table("/Users/phoebecunningham/Desktop/QMUL/Presentations/asv-table.txt", header=T, check.names = F)
rownames(asv_table) <- asv_table$sample_id
asv_table$sample_id <- NULL
```

Have a look at the datasets set you have brought in. What do they look like?

```{r looking at the data, results='hide'}
head(asv_table)
head(meta_data)
```

## Alpha Diversity 

For alpha diversity you can use vegan to calculate a vector which captures the diversity of each sample in one number. Usually, the higher the number the greater the diversity of the sample. Here we will capture the alpha diversity using Shannon'd vecor and Simpson diversity index which can both be generated using the function diversity() from the package vegan.

``` {r alpha diversity}
## Get the alpha diversity matrices
simpson <- vegan::diversity(asv_table, "simpson")
shannon <- vegan::diversity(asv_table, "shannon")
shannon.df <- as.data.frame(shannon)
shannon.df$sample_id <- rownames(shannon.df)
shannon.df <- merge(shannon.df, meta_data, by="sample_id")
```

## Beta diversity

For beta diversity you need to digest the data into distance matrices. These give you and indication of the dissimilarity between samples measure as aditance. You can create these in vegan using the function vegdist().

``` {r beta diversity}
## Get the distance matrices
jaccard <- vegdist(asv_table, method = "jaccard")
bray_curtis <- vegdist(asv_table, method = "bray")

pcoa <- pcoa(jaccard) # replace with bray once you have done jaccard
pcoa.vectors <- as.data.frame(pcoa$vectors[,1:2])
pcoa.values <- as.data.frame(pcoa$values)
pcoa.df <- cbind(pcoa.vectors, pcoa.values)
pcoa.df$loggedeigen <- pcoa.df$Eigenvalues^(1/3)
colnames(pcoa.df)[which(names(pcoa.df) == "Axis.1")] <- "PC1"
colnames(pcoa.df)[which(names(pcoa.df) == "Axis.2")] <- "PC2"
pcoa.df$sample_id <- row.names(pcoa.df)
pcoa.df <- merge(pcoa.df, meta_data, by="sample_id")
dm.pcoa.df <- pcoa.df
```    
 
## Visualing the diversity data 

For *Alpha Diversity* we will using boxplots to visualise the diversity for certain groups. Look at the metadata and think about which factors may be interesting to group the data by to understand trends. 

For *Beta Diversity* we have digested the data again into a Principle Coordinate Analysis, this takes the data from being multi-dimensional to being two-dimensional so we can visualise the distance between each sample in terms of the dissimilarity of microbial partners present between two samples. 

``` {r plotting the data}   
ggplot(dm.pcoa.df, aes(x = PC2, y = PC1, colour = feeding_guild)) +
  geom_point(size = 4) +
  geom_polygon(stat = "ellipse", aes(fill = feeding_guild), alpha = 0.05) +
  theme_classic()

```
