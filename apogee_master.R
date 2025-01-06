################################################################################
# Script: Analysis of apogee pipeline results
# Created by: Ishraq Akbar
# Data: 09/16/2024
################################################################################


################################################################################
#                       I.) SET UP PACKAGES & DEPENDENCIES
################################################################################
# install.packages("iNEXT")
# install.packages("ggvenn")
# install.packages("plotly")
# install.packages("ggraph")
# install.packages("igraph")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")
# install.packages("pacman")
# devtools::install_github("xia-lab/NetworkAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
# install.packages("pheatmap")

library(tidyverse)
library(vegan)
library(iNEXT)
library(ggrepel)
library(phyloseq)
library(ggpubr)
library(ggvenn)
library(gridExtra)



#...............................................................................



#...............................................................................
## Custom colour palettes:
# For stacked bar plots:
c25 <- c("#fdaf65", "#d886f3", "#a52f08", "#91d7f1", "#e5e327", "#ba69c3",
         "#6ce83c", "#27066f", "#fe8736", "#1f7264", "#0d1e42", 
         "#123c4e", "#065913", "#8de6cb", "#84fc5c", "#fe0a30", 
         "#773b70", "#08f4d8", "#911675", "#eb7126", "#54ad73", 
         "#286998", "#6c2520","#6744ff", "#b5ef31")


## For diversity plots:
treatment_colours <- c("Original DNA" = "#D5Ad36", "WGA DNA" = "#4488bf")



# Column orders
col_order <- c("Mock", "H2O", "Cpos", "Cneg", "C1", "C2", "C3", "C4", "N1", 
               "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10", "N11", 
               "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19", "N20", 
               "N21", "N22", "N23", "N24", "N25", "N26", "N27", "N28")
col_order2 <- c("MockNF", "H2ONF", "OG1", "OG2", "OG3", "OG4", "OG5", "OG6",
                "OG7", "NF1", "NF2", "NF3", "NF5", "NF6", "NF7")
col_order3 <- c("Mock", "H2O", "Cpos", "Cneg", paste0("C", 1:4), paste0("NF", c(1:3, 5:7)), 
                  paste0("OG", 1:7), paste0("N", 1:28))

#...............................................................................



#...............................................................................
## Functions:

# Function to clean pgpb csv Bacterial_strain column:
pgpb_cleaner <- function(x) {
  x <- iconv(x, "UTF-8", "ASCII", sub = "") # 1st removes non-ASCII characters
  parts <- str_split(x, " ")[[1]]   # First splits the column values by spaces
  sp_index <- which(parts == "sp.")   # Finds the index of 'sp.' for filtering
  
  if (length(sp_index) > 0) {
    # If sp. present, keeps preceding it and the 1st string after:
    result <- paste(parts[1:(sp_index + 1)], collapse = " ") 
  } else {
    # Otherwise, keeps only the first 2 strings (Genus + species):
    result <- paste(parts[1:2], collapse = " ")
  }
  return(result)
}

# Sorter function to reorganize taxonomy columns from Emu database:
col_sorter <- function(df) {
  df |>
    select(tax_id, superkingdom, phylum, class, order, family, genus, species, everything())
}

# Sorter function to organize columns based on numerical order
num_sorter <- function(df){
  cols <- colnames(df)   # Gets the column names
  fixed_cols <- cols[1:8] # Assigns the fixed columns I want to keep the same way
  sort_cols <- cols[-(1:8)] # Assigns the columns I want to manipulate (after col 8)
  sort_cols_number <- as.numeric(gsub("\\D", "", sort_cols)) # Extracts the numbers from the column names
  sort_cols <- sort_cols[order(sort_cols_number)] # Orders columns based on extracted numbers
  cols <- c(fixed_cols, sort_cols)   # Combines the fixed and sorted columns
  df <- df[, cols]   # Reorder the dataframe
  return(df) # Returns the dataframe
}


# ** Function to get the first string before any space or special character
get_first_string <- function(x) {
  str_to_lower(str_extract(x, "^[^\\s]+"))
}

# Function to check if any substring in mm$Species matches any substring in b2$Species
match_species <- function(b2_species, mm_species) {
  any(sapply(mm_species, function(x) grepl(x, b2_species, ignore.case = TRUE)))
}

# # Simpson diversity index calculator:
# simpson <- function(x){
#   n <- sum(x)
#   return(1 - (sum((x*(x-1))/(n*(n-1)))))
# }



#...............................................................................





################################################################################
#                             II.) DATA IMPORT
################################################################################
# Read files:

otu_master <- read_delim("data/otu_table_PW.csv", delim="\t") # Read otu table
tax_master <- read_csv("data/phyloseq_taxonomy_PW.csv") # Read phyloseq taxonomy table
pgpb <- read_csv('data/pgpb_tidy.csv') # PGPB genera
biomob <- read_csv('data/biomob_genome.csv') # Biomob PGPB taxonomy list
mock <- read_csv('data/mock-community-data.csv') # Information on mock community
b2 <- read_csv('data/emu-combined-tax_id-counts.csv')  # Contains absolute counts for Newfoundland samples
b2_app <- read_delim('data/otu_table_PW_newfoundland.csv', delim="\t")
b2_app_tax <- read_csv('data/phyloseq_taxonomy_PW_newfoundland.csv')
sampleinfo <- read_csv('data/PW_Metadata_batch2_batch3.csv')

################################################################################
#                            III.) DATA TIDYING
################################################################################

#...............................................................................
#### i.) Cleaning apogee sequencing data

## Rename first column of taxonomy dataframe to be just 'OTU':
tax_master <- rename(tax_master, OTU = '#OTU ID')
b2_app_tax <- rename(b2_app_tax, OTU = '#OTU ID')

## Filter out Eukarya from taxonomy list:
tax_master <- tax_master |>
  filter(Domain == "k__Bacteria")

b2_app_tax <- b2_app_tax |>
  filter(Domain == "k__Bacteria")


## Filter for OTUs above 70% mapping confidence
otu_master <- otu_master |>
  filter(MappingConfidence >= 70)
b2_app <- b2_app |>
  filter(MappingConfidence >= 70)

## Filter for taxa that are only present in filtered OTU list:
tax_master <- tax_master |>
  filter(OTU %in% otu_master$OTU) |>
  mutate(across(2:8, ~ gsub(".__", "", .))) |>
  mutate(across(2:8, ~ gsub("_", " ", .))) |> 
  mutate(across(2:8, ~ gsub("-", " ", .))) |>
  mutate(across(2:8, ~ ifelse(grepl("^(metagenome|uncultured)", .), "Other", .)))

b2_app_tax <- b2_app_tax |>
  filter(OTU %in% b2_app$OTU) |>
  mutate(across(2:8, ~ gsub(".__", "", .))) |>
  mutate(across(2:8, ~ gsub("_", " ", .))) |> 
  mutate(across(2:8, ~ gsub("-", " ", .))) |>
  mutate(across(2:8, ~ ifelse(grepl("^(metagenome|uncultured)", .), "Other", .)))
  



## Create metadata combining OTU information with taxonomy information:
# Rename M in OTU column to be Mock:
otu_master <- otu_master |>
  rename(Mock="PW_M") |>
  rename_with(~ sub("^PW_", "", .), 2:37) |>
  select(1, all_of(col_order), everything()) |>
  select(-c('TotalCount', 'MappingConfidence')) 

b2_app <- b2_app |>
  rename_with(~ sub("^PW_", "", .), 2:16) |>
  select(1, all_of(col_order2), everything()) |>
  select(-c('TotalCount', 'MappingConfidence')) 

otu_master <- otu_master |>
  filter(OTU %in% tax_master$OTU)
b2_app <- b2_app |>
  filter(OTU %in% b2_app_tax$OTU)


# Merge taxonomy column to otu_master to form otu_mm:
mm <- cbind(tax_master, otu_master) |>
  select(-9) 

b2_app_mm <- cbind(b2_app_tax, b2_app) |>
  select(-9)


# Create master dataframe for both groups:
b3_tmp <- mm |>
  pivot_longer(-c(1:8),
               names_to = "Sample",
               values_to = "Value")

b2_tmp <- b2_app_mm |>
  pivot_longer(-c(1:8),
               names_to = "Sample",
               values_to = "Value")

mm_all <- rbind(b3_tmp, b2_tmp) |>
  pivot_wider(names_from = "Sample", values_from = "Value")  |>
  select(-MockNF, -H2ONF) |>
  filter(Domain != 'Eukaryota') |>
  select(1:8, col_order3)
  

mm_all[is.na(mm_all)] <- 0
rm(b3_tmp, b2_tmp)

## Metadata for soil samples:
mm_soil <- mm_all |>
  select(1:8, c(c(NF1:NF3, NF5:NF7),N11:N16)) |>
  filter(rowSums(across(where(is.numeric)))!=0)

test <- mm_all |>
  select(1:8, c(C1:C4,N1, N1:N10, N17:N28))
  

## Metadata for wart samples:
mm_wart <- mm_all |>
  select(1:8, c(C1:C4,N1, N1:N10, N17:N28)) |>
  filter(rowSums(across(where(is.numeric)))!=0)


test <- mm_soil

test |> distinct(Phylum) |> nrow()
test |> distinct(Class) |> nrow()
test |> distinct(Order) |> nrow()
test |> distinct(Family) |> nrow()
test |> distinct(Genus) |> nrow()
test |> distinct(Species) |> nrow()

rm(test)

sum(colSums(mm_soil[,9:20]))




#...............................................................................



#...............................................................................
#### ii.) Cleaning Emu sequencing data

## Tidy Newfoundland sample dataframe:

# Set aside a dataframe for rarecurve analysis with mock and H2O controls:
b2_rare <- b2 |>
  col_sorter() |>
  num_sorter() |>
  select(1:8, ends_with("fastq")) |>
  rename(
    OG1 = barcode01.fastq, OG2 = barcode02.fastq, OG3 = barcode03.fastq,
    OG4 = barcode04.fastq, OG5 = barcode05.fastq, OG6 = barcode06.fastq,
    OG7 = barcode07.fastq, NF1 = barcode08.fastq, NF2 = barcode09.fastq,
    NF3 = barcode10.fastq, NF4 = barcode11.fastq, NF5 = barcode12.fastq,
    NF6 = barcode13.fastq, NF7 = barcode14.fastq, Mock = barcode15.fastq, 
    H2O = barcode16.fastq, OTU = tax_id,Domain = superkingdom, Phylum = phylum, 
    Class = class, Order = order, Family = family, Genus = genus, 
    Species = species
  ) |>
  mutate(across(9:24, round)) |>
  mutate(across(9:24, ~ ifelse(is.na(.), 0, .)))


# Does not contain information on Mock and H2O controls
b2 <- b2 |>
  col_sorter() |>
  num_sorter() |>
  select(1:8, ends_with("fastq"), -barcode15.fastq, -barcode16.fastq) |>
  rename(
    OG1 = barcode01.fastq, OG2 = barcode02.fastq, OG3 = barcode03.fastq,
    OG4 = barcode04.fastq, OG5 = barcode05.fastq, OG6 = barcode06.fastq,
    OG7 = barcode07.fastq, NF1 = barcode08.fastq, NF2 = barcode09.fastq,
    NF3 = barcode10.fastq, NF4 = barcode11.fastq, NF5 = barcode12.fastq,
    NF6 = barcode13.fastq, NF7 = barcode14.fastq, OTU = tax_id,
    Domain = superkingdom, Phylum = phylum, Class = class, Order = order,
    Family = family, Genus = genus, Species = species
  ) |>
  mutate(across(9:22, round)) |>
  mutate(across(9:22, ~ ifelse(is.na(.), 0, .)))

#...............................................................................



#...............................................................................
#### iii.) Combine and create metadata w/ Apogee + Emu data

# Convert dataframes to long format:

b2_long <- b2 |>
  pivot_longer(cols = starts_with("OG") | starts_with("NF"), 
               names_to = "Sample", 
               values_to = "Value")

mm_long <- mm |> 
  pivot_longer(cols = -c(OTU, Domain, Phylum, Class, Order, Family, Genus, Species), 
               names_to = "Sample", 
               values_to = "Value")


mm_temp <- bind_rows(b2_long, mm_long)
rm(b2_long, mm_long) # Remove temporary dfs

# Create master metadata list
metadata <- mm_temp %>%
  pivot_wider(names_from = Sample, values_from = Value) 

# Reorder columns. There are 58 columns total:
metadata <- metadata |>
  select(1:8, 23:24, 25:30, 9:22, 31:58)
rm(mm_temp) # Remove temporary df
#...............................................................................



#...............................................................................
#### iv.) Create additional dataframes for downstream analysis

### Taxonomy information for all metadata:
taxonomy <- metadata |>
  select(1:8)

### WGA vs original DNA treatments:
comparisons <- metadata |>
  select(1:8, 17:30)
#...............................................................................



#...............................................................................
#### iii.) Cleaning Mock Data
mock <- mock |>
  select(1,3) |>
  na.omit() |> # Get rid of NA values
  rename(Theoretical = '16S_only')|>
  rename(Species = 'species') |>
  arrange(desc(Species))

#...............................................................................


#...............................................................................
#### iv.) Cleaning PGPB dataframe:

# Cleans pgpb dataframe based on filter criteria (see function pgpb_cleaner):
pgpb <- pgpb |>
  mutate(Species = str_replace(Bacterial_Strain, "^[^\\s]+", Genus)) |> # Replaces the first string in Bacterial_Strain with Genus name
  mutate(Species = sapply(Species, pgpb_cleaner)) |> # Further cleans and stores results in Species column
  select(Species, Bacterial_Strain, PGPB_property, Genus)

#...............................................................................







################################################################################
#                            IV.) DATA VISUALIZATION
################################################################################

#...............................................................................
#### i.) Comparing sequencing results to Mock Data 
seq_mock <- mm |>
  select(8:9) |>
  filter(Species %in% mock$Species) |>
  group_by(Species) |>
  summarize(Mock = sum(Mock, na.rm=T))

# Calculate relative abundance of Mock community (sum is 25679):
seq_mock <- seq_mock |>
  mutate(Nanopore = round((Mock / sum(seq_mock$Mock)) * 100, 1)) |>
  arrange(desc(Species)) |>
  select(1,3)

# Merge dataframes together:
mock_dat <- mock |>
  left_join(seq_mock, by="Species")

# Convert dataframe to long format:
mock_dat <- mock_dat |>
  pivot_longer(cols = c("Nanopore", "Theoretical"),
               names_to = "sample", values_to = "percentage")



# Visualize mock comparisons:
mock_plot <- ggplot(mock_dat, aes(x = Species, y = percentage, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  
  labs(
    x = "Species", 
    y = "Relative Abundance (%)", 
    fill = "Sample",
    title = "ZymoBIOMICS D6305 Mock community recovery (positive control") +
  scale_fill_discrete(labels = c("Theoretical values", "Sequencing values")) +
  scale_fill_manual(values = c("maroon", "steelblue")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 25),
        axis.text.y = element_text(size = 30),
        axis.title.x = element_text(size = 30, margin = margin(t=15), face="bold"),
        axis.title.y = element_text(size = 30, margin = margin(r=15), face="bold"),
        plot.title = element_text(size = 30, face="bold", margin = margin(b=10)),
        axis.text.x.bottom = element_text(face="italic"),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30, face="bold"),
  ) +
  scale_x_discrete(labels = label_wrap_gen(width = 10))
mock_plot


#...............................................................................


#...............................................................................
#### ii.) Create a rarecurve of Apogee samples


# First ensure our dataframe only includes OTU IDs and samples:
mm_rare <- mm_all |>
  select('OTU', c(9:22, 30:57))

# Next assign the row names as OTUs and column names as samples:
## Setting OTU column values to character values
mm_all$OTU <- as.character(mm_all$OTU)

## Creating rownames
mm_rare <- mm_rare |>
  tibble::column_to_rownames("OTU") |>
  t() # Transpose

## Determine lowest number of counts for each sample:
raremax <- min(rowSums(mm_rare)) # Shows 180


## First run rarecurve using vegan package and store results into variable called rare_dat
rare_dat <- rarecurve(mm_rare, step=50, sample=raremax)

# Convert rare_dat to a data frame and add row names as a column
rare_df <- map_dfr(rare_dat, ~ as.data.frame(t(.x))) |>
  mutate(sample = rownames(mm_rare))

# Transform the data into a tidy format
rare_plot <- rare_df %>%
  pivot_longer(-sample, names_to = "step", values_to = "value") %>%
  drop_na() %>%
  mutate(step = as.numeric(str_replace(step, "N", "")))

# Visualize rarecurve plot provided by vegan
rare_plot <- rare_df |>
  pivot_longer(-sample, names_to = "step", values_to = "value") |>
  drop_na() |>
  mutate(step = as.numeric(str_replace(step, "N", "")))


# Create a new variable for Treatment and Category
rare_plot <- rare_plot |>
  mutate(Treatment = case_when(
    str_detect(sample, "^C") ~ "CFIA",
    str_detect(sample, "^NF|^OG") ~ "Newfoundland",
    str_detect(sample, "^N") ~ "Netherlands",
    str_detect(sample, "^Mock") ~ "Mock",
    str_detect(sample, "^H2O") ~ "H2O"
  ),
  Category = case_when(
    sample %in% c(paste0("C", 1:4), "N1", paste0("N", 6:10), paste0("N", 17:28)) ~ "Wart samples",
    sample %in% c(paste0("OG", 1:7), paste0("NF", 1:3), paste0("NF", 5:7), paste0("N", 2:5), paste0("N", 11:16)) ~ "Soil samples",
    TRUE ~ "Other"
  ))

# Get the last point for each sample
last_points <- rare_plot %>%
  group_by(sample) %>%
  filter(step == max(step))

# Plot using ggplot2
mm_rare_plot <- ggplot(rare_plot, aes(x = step, y = value, color = Category, shape = Treatment, group = sample)) +
  geom_vline(xintercept = raremax, color = "maroon", linewidth = 1, linetype = "longdash") +
  geom_line(linewidth = 0.8) +
  geom_text_repel(data = last_points, aes(x = step, y = value, label = sample), 
                  hjust = 0, vjust = 0.5, size = 6, max.overlaps = 60) +
  labs(title = "Figure 2. Rarefaction curves for 36 sequenced samples.", 
       x = "Number of Sequences", y = "Species Richness", 
       color = "Sample Type", shape = "Sample Source") +
  scale_color_manual(values = c("Wart samples" = "maroon", "Soil samples" = "steelblue", 
                                "Mock" = "red", "H2O" = "purple")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25, margin = margin(t = 15)),
    axis.title.y = element_text(size = 25, margin = margin(r = 15)),
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.x.bottom = element_text(face = "italic"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )

# View graph
mm_rare_plot

#...............................................................................

################################################################################
#                       ACCUMULATED RARECURVE FOR LOCATION
################################################################################

# Ensure the dataframe only includes OTU IDs and samples
mm_rare <- mm_all %>%
  select('OTU', c(9:22, 30:57)) %>%
  pivot_longer(cols = -OTU,
               names_to = "Sample",
               values_to = "Value") %>%
  mutate(Location = case_when(
    Sample %in% paste0("C", 1:4) ~ "CFIA",
    Sample %in% paste0("N", 1:28) ~ "Netherlands",
    Sample %in% paste0("NF", 1:7) ~ "Newfoundland",
    Sample == "Mock" ~ "Mock",
    Sample == "H2O" ~ "H2O"
  )) %>%
  select(OTU, Value, Location) %>%
  group_by(OTU, Location) %>%
  summarize(Value = sum(Value), .groups = 'drop') %>%
  pivot_wider(names_from = Location, values_from = Value, values_fill = 0) %>%
  select(1:6)

# Assign the row names as OTUs and column names as samples
mm_rare$OTU <- as.character(mm_all$OTU)

# Creating rownames
mm_rare <- mm_rare %>%
  column_to_rownames("OTU") %>%
  t() %>%
  as.data.frame()

# Determine lowest number of counts for each sample
raremax <- min(rowSums(mm_rare)) # Shows 180

# First run rarecurve using vegan package and store results into variable called rare_dat
rare_dat <- rarecurve(mm_rare, step=50, sample=raremax)

# Convert rare_dat to a data frame and add row names as a column
rare_df <- map_dfr(rare_dat, ~ as.data.frame(t(.x))) %>%
  mutate(sample = rownames(mm_rare))

# Transform the data into a tidy format
rare_plot <- rare_df %>%
  pivot_longer(-sample, names_to = "step", values_to = "value") %>%
  drop_na() %>%
  mutate(step = as.numeric(str_replace(step, "N", "")))

# Create a new variable for Treatment and Category
rare_plot <- rare_plot %>%
  mutate(Treatment = case_when(
    sample %in% rownames(mm_rare)[grep("CFIA", rownames(mm_rare))] ~ "CFIA",
    sample %in% rownames(mm_rare)[grep("Newfoundland", rownames(mm_rare))] ~ "Newfoundland",
    sample %in% rownames(mm_rare)[grep("Netherlands", rownames(mm_rare))] ~ "Netherlands",
    sample == "Mock" ~ "Mock",
    sample == "H2O" ~ "H2O"
  ),
  Category = case_when(
    sample %in% rownames(mm_rare)[grep("CFIA", rownames(mm_rare))] ~ "CFIA samples",
    sample %in% rownames(mm_rare)[grep("Newfoundland", rownames(mm_rare))] ~ "Newfoundland samples",
    sample %in% rownames(mm_rare)[grep("Netherlands", rownames(mm_rare))] ~ "Netherlands samples",
    sample == "Mock" ~ "Mock samples",
    sample == "H2O" ~ "H2O samples",
    TRUE ~ "Other"
  ))

# Get the last point for each sample
last_points <- rare_plot %>%
  group_by(sample) %>%
  filter(step == max(step))

# Plot using ggplot2
mm_location_rare <- ggplot(rare_plot, aes(x = step, y = value, color = Treatment, shape = Treatment, group = sample)) +
  geom_vline(xintercept = raremax, color = "maroon", linewidth = 1, linetype = "longdash") +
  geom_line(linewidth = 1) +
  geom_text_repel(data = last_points, aes(x = step, y = value, label = sample), 
                  hjust = 0, vjust = 0.5, size = 6, max.overlaps = 60) +
  labs(title = "Figure 2. Rarefaction curves for 36 sequenced samples.", 
       x = "Number of Sequences", y = "Species Richness", 
       color = "Location", shape = "Sample Source") +
  scale_color_manual(values = c("Netherlands" = "maroon", "Newfoundland" = "steelblue", 
                                "CFIA" = "darkgreen", "Mock" = "red", "H2O" = "purple")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25, margin = margin(t = 15)),
    axis.title.y = element_text(size = 25, margin = margin(r = 15)),
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.x.bottom = element_text(face = "italic"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )

# View graph
mm_location_rare



################################################################################
#                       ACCUMULATED RARECURVE FOR TREATMENT
################################################################################

# Ensure the dataframe only includes OTU IDs and samples
# Ensure the dataframe only includes OTU IDs and samples
mm_rare <- mm_all %>%
  select('OTU', c(9:22, 30:57)) %>%
  pivot_longer(cols = -OTU,
               names_to = "Sample",
               values_to = "Value") %>%
  mutate(Location = case_when(
    Sample %in% c(paste0("C", 1:4), paste0("N", 1:10), paste0("N", 17:28)) ~ "Wart",
    Sample %in% c(paste0("NF", 1:7), paste0("N", 11:16)) ~ "Soil",
    Sample == "Mock" ~ "Mock",
    Sample == "H2O" ~ "H2O"
  )) %>%
  select(OTU, Value, Location) %>%
  group_by(OTU, Location) %>%
  summarize(Value = sum(Value), .groups = 'drop') %>%
  pivot_wider(names_from = Location, values_from = Value, values_fill = 0) %>%
  select(1:5)

# Assign the row names as OTUs and column names as samples
mm_rare <- mm_rare %>%
  column_to_rownames("OTU") %>%
  t() %>%
  as.data.frame()

# Determine lowest number of counts for each sample
raremax <- min(rowSums(mm_rare)) # Shows 180

# Run rarecurve using vegan package and store results into variable called rare_dat
rare_dat <- rarecurve(mm_rare, step=50, sample=raremax)

# Convert rare_dat to a data frame and add row names as a column
rare_df <- map_dfr(rare_dat, ~ as.data.frame(t(.x))) %>%
  mutate(sample = rownames(mm_rare))

# Transform the data into a tidy format
rare_plot <- rare_df %>%
  pivot_longer(-sample, names_to = "step", values_to = "value") %>%
  drop_na() %>%
  mutate(step = as.numeric(str_replace(step, "N", "")))

# Create a new variable for Treatment
rare_plot <- rare_plot %>%
  mutate(Treatment = case_when(
    str_detect(sample, "Wart") ~ "Wart",
    str_detect(sample, "Soil") ~ "Soil",
    sample == "Mock" ~ "Mock",
    sample == "H2O" ~ "H2O"
  ))

# Get the last point for each sample
last_points <- rare_plot %>%
  group_by(sample) %>%
  filter(step == max(step))

# Plot using ggplot2
mm_rare_type <- ggplot(rare_plot, aes(x = step, y = value, color = Treatment, shape = Treatment, group = sample)) +
  geom_vline(xintercept = raremax, color = "maroon", linewidth = 1, linetype = "longdash") +
  geom_line(linewidth = 2) +  # Increase line width
  geom_text_repel(data = last_points, aes(x = step, y = value, label = sample), 
                  hjust = 0, vjust = 0.5, size = 6, max.overlaps = 60) +
  labs(title = "Figure 2. Rarefaction curves for selected samples.", 
       x = "Number of Sequences", y = "Species Richness", 
       color = "Sample Type", shape = "Sample Type") +
  scale_color_manual(values = c("Wart" = "maroon", "Soil" = "steelblue", 
                                "Mock" = "red", "H2O" = "purple")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25, margin = margin(t = 15)),
    axis.title.y = element_text(size = 25, margin = margin(r = 15)),
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.x.bottom = element_text(face = "italic"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  )


# View graph
mm_rare_type



#...............................................................................
#### ii.) Create a rarecurve of Newfoundland samples

# First ensure our dataframe only includes OTU IDs and samples:
b2_rare <- mm_all |>
  select('OTU', c(Mock, H2O, OG1:OG7, NF1:NF3, NF5:NF7, N2:N5, N11:N16)) |>
  

# Next assign the row names as OTUs and column names as samples:
## Setting OTU column values to character values
b2$OTU <- as.character(b2$OTU)

## Creating rownames
rownames(b2_rare) <- b2$OTU
b2_rare <- b2_rare |>
  select(- 'OTU')

## Transposing dataframe:
b2_rare <- t(b2_rare)


## Determine lowest number of counts for each sample:
raremax <- min(rowSums(b2_rare)) # Shows 74,683


## First run rarecurve using vegan package and store results into variable called rare_dat
rare_dat <- rarecurve(b2_rare, step=50, sample=raremax)

# Convert rare_dat to a data frame and add row names as a column
rare_df <- map_dfr(rare_dat, ~ as.data.frame(t(.x))) |>
  mutate(sample = rownames(b2_rare))

# Transform the data into a tidy format
rare_plot <- rare_df %>%
  pivot_longer(-sample, names_to = "step", values_to = "value") %>%
  drop_na() %>%
  mutate(step = as.numeric(str_replace(step, "N", "")))

# Visualize rarecurve plot provided by vegan
rare_plot <- rare_df |>
  pivot_longer(-sample, names_to = "step", values_to = "value") |>
  drop_na() |>
  mutate(step = as.numeric(str_replace(step, "N", "")))

# Create a new variable for colour
rare_plot <- rare_plot |>
  mutate(Treatment = case_when(
    str_detect(sample, "^OG") ~ "Original",
    str_detect(sample, "^NF") ~ "WGA treated",
    str_detect(sample, "^Mock") ~ "Mock",
    str_detect(sample, "^H2O") ~ "H2O"
  ))

# Get the last point for each sample
last_points <- rare_plot |>
  group_by(sample) |>
  filter(step == max(step))

# Plot using ggplot2
b2_rare_plot <- ggplot(rare_plot, aes(x = step, y = value, color = Treatment, group = sample)) +
  geom_vline(xintercept = raremax, color = "maroon", linewidth = 1, linetype = "longdash") +
  geom_line(linewidth = 0.8) +
  geom_text_repel(data = last_points, aes(x = step, y = value, label = sample), 
                  hjust = 0, vjust = 0.5, size = 6, max.overlaps = 30) +
  labs(title = "",
       x = "Number of sequences",
       y = "Species richness") +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 25, margin = margin(t = 15)),
    axis.title.y = element_text(size = 25, margin = margin(r = 15)),
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.x.bottom = element_text(face = "italic"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25)
  ) +
  scale_x_continuous(labels = scales::comma)

b2_rare_plot
#...............................................................................




#...............................................................................
## ALPHA DIVERSITY ON DIFFERENT TREATMENTS
#...............................................................................


## From Riffomonias: Want column for sample, column for taxa, column for count
b2_div <- b2 |>
  select('OTU', 9:22)

## Manipulate dataframe for calculations:
b2_wide <- b2_div |>
  as.data.frame() |>
  select(2:15)
rownames(b2_wide) <- b2_div$OTU
b2_wide <- b2_wide[rowSums(b2_wide)>0,] |> t()

## Calculate indices:
# Calculate true ginni-simpson index:
treatment_index <- 1/(vegan::diversity(b2_wide, "simpson")) 
treatment_index <- data.frame(GinniSimpson = treatment_index) |>
  mutate(Sample = rownames(b2_wide)) |>
  select(Sample, GinniSimpson)

# Calculate chao1 index:
chao1_results <- estimateR(b2_wide) |> t() |> as.data.frame() |>
  rename("Sobs" = 1, "Chao1" = 2, "Chao1se" = 3, "ACE" = 4, "ACEse" = 5)

# Calculate shannon index & combine w/ chao1 index:
treatment_index <- treatment_index |>
  mutate(ShannonTD = exp(vegan::diversity(b2_wide, "shannon"))) |>
  bind_cols(chao1_results)

# Create new column for treatment:
treatment_index <- treatment_index |>
  mutate(Treatment = ifelse(str_starts(Sample, "NF"), "WGA DNA", "Original DNA")) |>
  select(Sample, Treatment, 2:8) |>
  group_by(Sample)

# Revert rownames to default:
rownames(treatment_index) <- NULL

## Visualize plots:

# Shannon TD violin plot for treatment
treatment_shannon <- ggplot(treatment_index, aes(x = Treatment, y = ShannonTD, fill = Treatment)) +
  geom_violin(trim = FALSE, col = "#021a38", linewidth = 1, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.5), fill = "#021a38", outlier.shape = NA, show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.05, color = "white", linewidth = 0.5, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_y_continuous(breaks = seq(round(min(treatment_index$ShannonTD, na.rm=T), digits = -2) - 100,
                                  round(max(treatment_index$ShannonTD, na.rm=T), digits = -2) + 100, length.out = 6)) +
  ggtitle("Shannon True Diversity") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)))

# Ginni Simpson violin plot for treatment
treatment_simpson <- ggplot(treatment_index, aes(x = Treatment, y = GinniSimpson, fill = Treatment)) +
  geom_violin(trim = FALSE, col = "#021a38", linewidth = 1, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.5), fill = "#021a38", outlier.shape = NA, show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.05, color = "white", linewidth = 0.5, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_y_continuous(breaks = seq(round(min(treatment_index$GinniSimpson, na.rm=T), digits=1) - 0.05,
                                  round(max(treatment_index$GinniSimpson, na.rm=T), digits=1) + 0.05, length.out = 6)) +
  ggtitle("Ginni Simpson") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)))

# Chao1 violin plot for treatment
treatment_chao1 <- ggplot(treatment_index, aes(x = Treatment, y = Chao1, fill = Treatment)) +
  geom_violin(trim = FALSE, col = "#021a38", linewidth = 1, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.5), fill = "#021a38", outlier.shape = NA, show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.05, color = "white", linewidth = 0.5, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_y_continuous(breaks = seq(round(min(treatment_index$Chao1, na.rm=T), digits = -2) - 100,
                                  round(max(treatment_index$Chao1, na.rm=T), digits = -2) + 100, length.out = 6)) +
  ggtitle("Chao1") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)))



# ACE (Abundance based-coverage estimator) violin plot for treatment
treatment_ACE <- ggplot(treatment_index, aes(x = Treatment, y = ACE, fill = Treatment)) +
  geom_violin(trim = FALSE, col = "#021a38", linewidth = 1, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  geom_boxplot(width = 0.05, position = position_dodge(width = 0.5), fill = "#021a38", outlier.shape = NA, show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.05, color = "white", linewidth = 0.5, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_y_continuous(breaks = seq(round(min(treatment_index$ACE, na.rm=T), digits = -2) - 100,
                                  round(max(treatment_index$ACE, na.rm=T), digits = -2) + 100, length.out = 6)) +
  ggtitle("Abundance based-coverage estimator") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)))

# Faceted graph for all of them:
treatment_plots <- gridExtra::grid.arrange(treatment_shannon, treatment_simpson, treatment_chao1, treatment_ACE, ncol = 2)


################################################################################
####TESTING ONLY#############
# Shannon TD boxplot for treatment
# Shannon TD boxplot for treatment
treatment_shannon <- ggplot(treatment_index, aes(x = Treatment, y = ShannonTD, fill = Treatment)) +
  geom_boxplot(outlier.shape = 16, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.05, color = "white", linewidth = 0.5, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  geom_signif(comparisons = list(c("Treatment1", "Treatment2")), map_signif_level = TRUE, test = "t.test") +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_y_continuous(breaks = seq(round(min(treatment_index$ShannonTD, na.rm=T), digits = -2) - 100,
                                  round(max(treatment_index$ShannonTD, na.rm=T), digits = -2) + 100, length.out = 6)) +
  ggtitle("Shannon True Diversity") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)))

# Ginni Simpson boxplot for treatment
treatment_simpson <- ggplot(treatment_index, aes(x = Treatment, y = GinniSimpson, fill = Treatment)) +
  geom_boxplot(outlier.shape = 16, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.05, color = "white", linewidth = 0.5, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  geom_signif(comparisons = list(c("Treatment1", "Treatment2")), map_signif_level = TRUE, test = "t.test") +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_y_continuous(breaks = seq(round(min(treatment_index$GinniSimpson, na.rm=T), digits=1) - 0.05,
                                  round(max(treatment_index$GinniSimpson, na.rm=T), digits=1) + 0.05, length.out = 6)) +
  ggtitle("Ginni Simpson") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)))

# Chao1 boxplot for treatment
treatment_chao1 <- ggplot(treatment_index, aes(x = Treatment, y = Chao1, fill = Treatment)) +
  geom_boxplot(outlier.shape = 16, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.05, color = "white", linewidth = 0.5, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  geom_signif(comparisons = list(c("Treatment1", "Treatment2")), map_signif_level = TRUE, test = "t.test") +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_y_continuous(breaks = seq(round(min(treatment_index$Chao1, na.rm=T), digits = -2) - 100,
                                  round(max(treatment_index$Chao1, na.rm=T), digits = -2) + 100, length.out = 6)) +
  ggtitle("Chao1") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)))

# ACE (Abundance based-coverage estimator) boxplot for treatment
treatment_ACE <- ggplot(treatment_index, aes(x = Treatment, y = ACE, fill = Treatment)) +
  geom_boxplot(outlier.shape = 16, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.05, color = "white", linewidth = 0.5, show.legend = FALSE, position = position_dodge(width = 0.5)) +
  geom_signif(comparisons = list(c("Treatment1", "Treatment2")), map_signif_level = TRUE, test = "t.test") +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_y_continuous(breaks = seq(round(min(treatment_index$ACE, na.rm=T), digits = -2) - 100,
                                  round(max(treatment_index$ACE, na.rm=T), digits = -2) + 100, length.out = 6)) +
  ggtitle("Abundance based-coverage estimator") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(r = 10)))

# Faceted graph for all of them:
treatment_plots <- gridExtra::grid.arrange(treatment_shannon, treatment_simpson, treatment_chao1, treatment_ACE, ncol = 2)





#---------------


## Reformat dataframes:
alpha_NF <- as.data.frame(b2_wide) |>
  slice(8:14)

alpha_CF <- mm |>
  select(OTU, C1:C4) |>
  tibble::column_to_rownames("OTU") |>
  filter(rowSums(across(where(is.numeric)))!=0) |>
  t()

alpha_NL <- mm |>
  select(OTU, N1:N28) |>
  tibble::column_to_rownames("OTU") |>
  filter(rowSums(across(where(is.numeric)))!=0) |>
  t()

alpha_soil <- mm_soil |>
  select(OTU, 9:20) |>
  tibble::column_to_rownames("OTU") |>
  t()

alpha_wart <- mm_wart |>
  select(OTU, 9:34) |>
  tibble::column_to_rownames("OTU") |>
  t()

## Calculate indicies for soil:
shannon_soil <- exp(vegan::diversity(alpha_soil, index = "shannon"))
chao1_soil <- estimateR(alpha_soil)["S.chao1", ]
simpson_soil <- 1/(vegan::diversity(alpha_soil, "simpson"))
alpha_soil_indices <- data.frame(
  ShannonTD = shannon_soil,
  Chao1 = chao1_soil,
  SimpsonTD = simpson_soil) |>
  mutate(Type = "Soil")

## Calculate indicies for wart:
shannon_wart <- exp(vegan::diversity(alpha_wart, index = "shannon"))
chao1_wart <- estimateR(alpha_wart)["S.chao1", ]
simpson_wart <- 1/(vegan::diversity(alpha_wart, "simpson"))

alpha_wart_indices <- data.frame(
  ShannonTD = shannon_wart,
  Chao1 = chao1_wart,
  SimpsonTD = simpson_wart) |>
  mutate(Type = "Wart")


# Merge dataframe into one:
div_types <- rbind(alpha_soil_indices, alpha_wart_indices)


# ShannonTD plot
ggplot(div_types, aes(x = Type, y = ShannonTD, fill = Type)) +
  geom_boxplot(alpha = 0.5, size = 1) +  # Box plot with adjusted alpha and outline thickness
  geom_jitter(position = position_jitter(width = 0.2), alpha = 1, color = "black", size = 2) +  # Increase point size
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = position_dodge(0.75), size = 1) +
  theme_classic(base_size = 16) +
  labs(title = "ShannonTD", x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40, face = "bold"),
        axis.text.y = element_text(size = 30)) +
  scale_fill_manual(values = c("Soil" = "darkgreen", "Wart" = "maroon"))


# Chao1 plot
ggplot(div_types, aes(x = Type, y = Chao1, fill = Type)) +
  geom_boxplot(alpha = 0.5, size = 1) +  # Box plot with adjusted alpha and outline thickness
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, color = "black", size = 2) +  # Jitter to show individual points
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = position_dodge(0.75), size = 1) +
  theme_classic(base_size = 16) +
  labs(title = "Chao1", x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(size = 30, face = "bold"),
        axis.text.y = element_text(size = 30),
        legend.position = 'none') +
  scale_fill_manual(values = c("Soil" = "darkgreen", "Wart" = "maroon"))






#------------



# Remove non-finite values:
div_samples_clean <- div_samples %>%
  filter(is.finite(ShannonTD) & is.finite(GinniSimpson) & is.finite(Chao1) & is.finite(ACE))

# Calculate indices for CF
CF_shannon <- exp(vegan::diversity(alpha_CF, index = "shannon"))
CF_simpson <- 1/(vegan::diversity(alpha_CF, index = "simpson"))
CF_chao1 <- estimateR(alpha_CF)["S.chao1", ]
CF_ace <- estimateR(alpha_CF)["S.ACE", ]


alpha_CF_indices <- data.frame(
  ShannonTD = CF_shannon,
  GinniSimpson = CF_simpson,
  Chao1 = CF_chao1,
  ACE = CF_ace) |>
  mutate(Location = "CFIA")


# Calculate indices for NF
NF_shannon <- exp(vegan::diversity(alpha_NF, index = "shannon"))
NF_simpson <- 1/(vegan::diversity(alpha_NF, index = "simpson"))
NF_chao1 <- estimateR(alpha_NF)["S.chao1", ]
NF_ace <- estimateR(alpha_NF)["S.ACE", ]

alpha_NF_indices <- data.frame(
  ShannonTD = NF_shannon,
  GinniSimpson = NF_simpson,
  Chao1 = NF_chao1,
  ACE = NF_ace) |>
  mutate(Location = "Newfoundland") 

# Calculate indices for NL
NL_shannon <- exp(vegan::diversity(alpha_NL, index = "shannon"))
NL_simpson <- 1/(vegan::diversity(alpha_NL, index = "simpson"))
NL_chao1 <- estimateR(alpha_NL)["S.chao1", ]
NL_ace <- estimateR(alpha_NL)["S.ACE", ]

alpha_NL_indices <- data.frame(
  ShannonTD = NL_shannon,
  GinniSimpson = NL_simpson,
  Chao1 = NL_chao1,
  ACE = NL_ace) |>
  mutate(Location = "Netherlands")


# Merge dataframe into one:

div_samples <- rbind(alpha_CF_indices, alpha_NF_indices, alpha_NL_indices)

# Remove non-finite values:
div_samples_clean <- div_samples %>%
  filter(is.finite(ShannonTD) & is.finite(GinniSimpson) & is.finite(Chao1) & is.finite(ACE))


# ShannonTD
p1 <- ggplot(div_samples_clean, aes(x = Location, y = ShannonTD, fill = Location)) +
  geom_violin(trim = FALSE) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = position_dodge(0.75), size = 1) +
  theme_classic(base_size = 16) +
  labs(title = "ShannonTD", x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# GinniSimpson
p2 <- ggplot(div_samples_clean, aes(x = Location, y = GinniSimpson, fill = Location)) +
  geom_violin(trim = FALSE) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = position_dodge(0.75), size = 1) +
  theme_classic(base_size = 16) +
  labs(title = "GinniSimpson", x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Chao1
p3 <- ggplot(div_samples_clean, aes(x = Location, y = Chao1, fill = Location)) +
  geom_violin(trim = FALSE) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = position_dodge(0.75), size = 1) +
  theme_classic(base_size = 16) +
  labs(title = "Chao1",  x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ACE
p4 <- ggplot(div_samples_clean, aes(x = Location, y = ACE, fill = Location)) +
  geom_violin(trim = FALSE) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.5, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = position_dodge(0.75), size = 1) +
  theme_classic(base_size = 16) +
  labs(title = "Abundance-based Coverage Estimator (ACE)",  x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Arrange the plots in a grid with larger size
grid.arrange(p1, p2, p3, p4, ncol = 2, top = element_text("Alpha Diversity Indices by Location", size=20), heights = c(1, 1), widths = c(1, 1))

#-------------------------------------------------------------------------------


#...............................................................................


## Determine if the data displays normality (p>0.05 or p>0.01):

# Shannon TD:
shapiro.test(treatment_index$ShannonTD[treatment_index$Treatment == "Original DNA"]) # p=0.341
shapiro.test(treatment_index$ShannonTD[treatment_index$Treatment == "WGA DNA"]) # p=0.5066

# Ginni Simpson:
shapiro.test(treatment_index$GinniSimpson[treatment_index$Treatment == "Original DNA"]) # p=0.05289
shapiro.test(treatment_index$GinniSimpson[treatment_index$Treatment == "WGA DNA"]) # p=0.1193

# Chao1:
shapiro.test(treatment_index$Chao1[treatment_index$Treatment == "Original DNA"]) # p=0.5965
shapiro.test(treatment_index$Chao1[treatment_index$Treatment == "WGA DNA"]) # p=0.6894


# ACE:
shapiro.test(treatment_index$ACE[treatment_index$Treatment == "Original DNA"]) # p=0.1791
shapiro.test(treatment_index$ACE[treatment_index$Treatment == "WGA DNA"]) # p=0.9304

# Soil:
shapiro.test(div_types$ShannonTD[div_types$Type == "Soil"]) 
shapiro.test(div_types$Chao1[div_types$Type == "Soil"])

# Wart:
shapiro.test(div_types$ShannonTD[div_types$Type == "Wart"]) 
shapiro.test(div_types$Chao1[div_types$Type == "Wart"])

## Determine if statistical significance between the two groups using t-test:
t.test(data = treatment_index, ShannonTD ~ Treatment) # p=0.003377 (p < 0.05), there is a difference
t.test(data = treatment_index, GinniSimpson ~ Treatment) # p = 0.176, there is no difference
t.test(data = treatment_index, Chao1 ~ Treatment) # p =  p-value = 0.0004332, there is a difference
t.test(data = treatment_index, ACE ~ Treatment) # p-value = 0.002545, there is a difference

t.test(data = div_types, ShannonTD ~ Type) # p-value = 0.0106, no difference
t.test(data = div_types, SimpsonTD ~ Type) # p-value = 0.6117, no difference


# Results indicate a difference in species richness (Chao1 + ACE) and evenness (Shannon TD) between the two groups
# However, the dominant species present in both treatments remain largely unchanged in identification of taxa (Ginni Simpson)


#-------------------------------------------------------------------------------
# BETA DIVERSITY ANALYSIS:
#-------------------------------------------------------------------------------

# Tidy data and create new dataframe:
treatment_dat <- b2_div |>
  select(-OTU) |>
  as.data.frame()
row.names(treatment_dat) <- b2_div$OTU

# Transpose to have rows as samples, columns as taxa
treatment_dat <- t(treatment_dat) 

# Perform NMDS:
set.seed(19990115)
nmds_treatment <- metaMDS(treatment_dat, distance="robust.aitchison", k=2, trymax=999)

# Extract NMDS scores:
nmds_treatment_scores <- as.data.frame(scores(nmds_treatment, display = "sites"))
nmds_treatment_scores$Sample <- rownames(nmds_treatment_scores)

# Create a new dataframe for NMDS graph
meta_treatment <- data.frame(
  Sample = rownames(treatment_dat),
  Treatment = ifelse(grepl("OG", rownames(treatment_dat)), "Original DNA", "WGA DNA")
)

# Merge NMDS scores with metadata:
nmds_treatment_scores <- merge(nmds_treatment_scores, meta_treatment, by = "Sample")
nmds_treatment_scores <- nmds_treatment_scores |>
  mutate(Sample = str_replace(Sample, "^NF", "WGA"))


# Code for plot:
nmds_plot_treatment <- ggplot(data = nmds_treatment_scores, aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(size = 5) +
  geom_text(aes(label = Sample), vjust = 1.5, hjust = 0.5, size = 5, show.legend = F) +
  stat_ellipse(aes(fill = Treatment), type = "norm", alpha = 0.2, geom = "polygon", color = "black", show.legend = F) +
  theme_classic() +
  labs(title = "NMDS analysis of WGA vs original treated DNA microbiome profile",
       x = "NMDS1", y = "NMDS2") +
  scale_color_manual(values = c("purple", "steelblue")) +
  scale_fill_manual(values = c("purple", "steelblue")) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20)
  ) +
  annotate("text", x = Inf, y = Inf, label = "Stress = 0.128", hjust = 1, vjust = 2, size = 7, fontface = "bold")


# Visualize:
nmds_plot_treatment
#...............................................................................
# DETERMINING SHARED SPECIES BETWEEN 3 LOCATIONS
#...............................................................................

## First filter data sets:
# Convert to long format:
mm_spec <- mm |>
  mutate(Species = ifelse(Species == "", "Unclassified", Species)) |>
  select(Species, 13:44) |>
  group_by(Species) |>
  summarise(across(everything(), sum, na.rm = TRUE)) |>
  pivot_longer(
    cols = -Species,
    names_to = "Sample",
    values_to = "Counts"
  ) |>
  mutate(Species = str_remove_all(Species, "\\[|\\]")) |>
  mutate(Genus = str_extract(Species, "^[^ ]+")) |>
  select(Genus, Species, Sample, Counts)


b2_spec <- b2 |>
  mutate(Species = ifelse(Species == "", "Unclassified", Species)) |>
  select(Species, 9:22) |>
  group_by(Species) |>
  summarise(across(everything(), sum, na.rm = TRUE)) |>
  pivot_longer(
    cols = -Species,
    names_to = "Sample",
    values_to = "Counts"
  ) |>
  mutate(Species = str_remove_all(Species, "\\[|\\]")) |> # Removes square brackets, IMPORTANT
  mutate(Genus = str_extract(Species, "^[^ ]+")) |>
  select(Genus, Species, Sample, Counts)

# Join dataframes together:
all_spec <- bind_rows(mm_spec, b2_spec) |>
  filter(!Species %in% c(NA, "Other", "Unclassified", "bacterium", 
                         "Phage", "unidentified")) |>
  filter(!Genus %in% c(NA, "Other", "Unclassified", "bacterium", 
                       "Phage", "unidentified")) |>
  filter(!grepl("^[a-z]", Genus))

# Filter species for each sample type:
cmntax_C <- all_spec |> 
  filter(grepl("^C[0-9]", Sample))  |> filter(Counts > 0) |>
  select(Genus, Species) |> distinct()
cmntax_NF <- all_spec |> 
  filter(grepl("^NF[0-9]", Sample)) |> filter(Counts > 0) |>
  select(Genus, Species) |> distinct()
cmntax_N <- all_spec |> 
  filter(grepl("^N[0-9]", Sample)) |> filter(Counts > 0) |>
  select(Genus, Species) |> distinct()

# Create a list of genus from each site:
gen_C <- cmntax_C$Genus|> unique() 
gen_NF <- cmntax_NF$Genus |> unique() 
gen_N <- cmntax_N$Genus |> unique() 

# Reshape cmntax files to just have species now:
cmntax_C <- cmntax_C$Species
cmntax_NF <- cmntax_NF$Species
cmntax_N <- cmntax_N$Species

# Create lists:
spec_list <- list(CFIA=cmntax_C, Newfoundland=cmntax_NF, Netherlands=cmntax_N)
gen_list <- list(CFIA=gen_C, Newfoundland=gen_NF, Netherlands=gen_N)


## Visualize plots:
# Species plot:
venn_spec <- ggvenn(spec_list,
                    fill_color = c("red", "blue", "green"),
                    fill_alpha = 0.35,
                    stroke_size = 1, stroke_linetype = "solid",
                    set_name_size = 10, text_size = 7.1)
venn_spec

# Genus plot:
venn_gen <- ggvenn(gen_list,
                   fill_color = c("red", "blue", "green"),
                   fill_alpha = 0.35,
                   stroke_size = 1, stroke_linetype = "solid",
                   set_name_size = 10, text_size = 7.1)
venn_gen

#...............................................................................



#...............................................................................
### CREATING STACKED BAR PLOTS

### Individual bar plots for each location

## Calculate top 20 genera:
gen20_C <- all_spec |>
  filter(grepl("^C[0-9]", Sample)) |> # Grepling for desired location
  filter(Counts > 0) |> # Filtering out counts that are 0 or lower
  group_by(Genus) |>
  summarize(TotalCounts = sum(Counts, na.rm=T), .groups="drop") |>
  top_n(20, TotalCounts) # Grabbing top 20 genera for the desired location

gen20_NF <- all_spec |>
  filter(grepl("^NF", Sample)) |> # Grepling for desired location
  filter(Counts > 0) |> # Filtering out counts that are 0 or lower
  group_by(Genus) |>
  summarize(TotalCounts = sum(Counts, na.rm=T), .groups="drop") |>
  top_n(20, TotalCounts) # Grabbing top 20 genera for the desired location

gen20_N <- all_spec |>
  filter(grepl("^N[0-9]", Sample)) |> # Grepling for desired location
  filter(Counts > 0) |> # Filtering out counts that are 0 or lower
  group_by(Genus) |>
  summarize(TotalCounts = sum(Counts, na.rm=T), .groups="drop") |>
  top_n(20, TotalCounts) # Grabbing top 20 genera for the desired location


## Calculate total counts per sample:
n_genC <- all_spec |>
  filter(grepl("^C[0-9]", Sample)) |> # Grepling for desired location
  group_by(Sample) |>
  summarize(n = sum(Counts, na.rm=T), .groups="drop") |>
  as.data.frame()

n_genNF <- all_spec |>
  filter(grepl("^NF", Sample)) |> # Grepling for desired location
  group_by(Sample) |>
  summarize(n = sum(Counts, na.rm=T), .groups="drop") |>
  as.data.frame()

n_genN <- all_spec |>
  filter(grepl("^N[0-9]", Sample)) |> # Grepling for desired location
  group_by(Sample) |>
  summarize(n = sum(Counts, na.rm=T), .groups="drop") |>
  as.data.frame()


## Rearrange dataframes to have those top 20 genera:
gendat_C <- all_spec |>
  filter(Genus %in% gen20_C$Genus,
         grepl("^C[0-9]", Sample),
         Counts > 0) |>
  as.data.frame()

gendat_NF <- all_spec |>
  filter(Genus %in% gen20_NF$Genus,
         grepl("^NF", Sample),
         Counts > 0) |>
  as.data.frame()

gendat_N <- all_spec |>
  filter(Genus %in% gen20_N$Genus,
         grepl("^N[0-9]", Sample),
         Counts > 0) |>
  as.data.frame()





## Merge dataframe with total counts, and calculate relative abundance:
gendat_C <- gendat_C |>
  left_join(n_genC, by="Sample") |>
  mutate(RelAbund = (Counts/n) * 100) |>
  select(-n)

gendat_NF <- gendat_NF |>
  left_join(n_genNF, by="Sample") |>
  mutate(RelAbund = (Counts/n) * 100) |>
  select(-n)

gendat_N <- gendat_N |>
  left_join(n_genN, by="Sample") |>
  mutate(RelAbund = (Counts/n) * 100) |>
  select(-n)




# Calculate the total relative abundance for each genus across all samples
relabund_gen_C <- gendat_C |>
  group_by(Genus) |>
  summarize(TotalRelAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop")
relabund_gen_NF <- gendat_NF |>
  group_by(Genus) |>
  summarize(TotalRelAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop")
relabund_gen_N <- gendat_N |>
  group_by(Genus) |>
  summarize(TotalRelAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop")

# Reorder the Genus factor levels based on total relative abundance
gendat_C <- gendat_C |>
  mutate(Genus = factor(Genus, levels = relabund_gen_C$Genus[order(-relabund_gen_C$TotalRelAbund)]))
gendat_NF <- gendat_NF |>
  mutate(Genus = factor(Genus, levels = relabund_gen_NF$Genus[order(-relabund_gen_NF$TotalRelAbund)]))
gendat_N <- gendat_N |>
  mutate(Genus = factor(Genus, levels = relabund_gen_N$Genus[order(-relabund_gen_N$TotalRelAbund)])) |>
  mutate(Sample_num = as.numeric(sub("N", "", Sample))) |>
  arrange(Sample_num) |>
  mutate(Sample = factor(Sample, levels = unique(Sample)))


## Barplot codes:
genbar_C <- ggplot(gendat_C, aes(x = Sample, y = RelAbund, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative Abundance of the top 20 genera in PEI samples",
       x = "Sample",
       y = "Relative Abundance (%)") +
  scale_fill_manual(values = c25) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size = 20, hjust = 1),
        legend.position = "right",
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold", margin = margin(b = 10)),
        axis.title = element_text(size = 20, face = "bold", margin = margin(b = 10)),
        legend.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(face = "italic", size = 15))

genbar_NF <- ggplot(gendat_NF, aes(x = Sample, y = RelAbund, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative Abundance of the top 20 genera in Newfoundland samples",
       x = "Sample",
       y = "Relative Abundance (%)") +
  scale_fill_manual(values = c25) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size=20, hjust = 1),
        legend.position = "right",
        axis.text.y = element_text(size=20),
        plot.title = element_text(size=20, face="bold", margin=margin(b=10)),
        axis.title = element_text(size=20, face="bold", margin=margin(b=10)),
        legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(face = "italic", size=15))

genbar_N <- ggplot(gendat_N, aes(x = Sample, y = RelAbund, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative Abundance of the top 20 genera in Netherland samples",
       x = "Sample",
       y = "Relative Abundance (%)") +
  scale_fill_manual(values = c25) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, size=20, hjust = 1),
        legend.position = "right",
        axis.text.y = element_text(size=20),
        plot.title = element_text(size=20, face="bold", margin=margin(b=10)),
        axis.title = element_text(size=20, face="bold", margin=margin(b=10)),
        legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(face = "italic", size=15))


###############################################################################

## First calculate total genera per sample (not including Other or NA):

## First determine genera found in all samples
gen_samples <- all_spec |>
  group_by(Genus, Sample) |>
  summarize(Count = sum(Counts, na.rm = TRUE), .groups="drop") |>
  ungroup() |>
  group_by(Genus) |>
  summarize(SampleCount = n(), .groups = "drop") |>
  filter(SampleCount == n_distinct(all_spec$Sample)) |>
  pull(Genus)

## Calculate total counts for genera present in all samples
gen_counts <- all_spec |>
  filter(Genus %in% gen_samples) |>
  group_by(Genus) |>
  summarize(TotalCounts = sum(Counts, na.rm = TRUE), .groups = "drop")

## Calculate top 20 genera in the samples:
gen_20 <- gen_counts |>
  top_n(20, TotalCounts) |>
  pull(Genus)

## Filter the original dataframe to include only the top 20 genera
all_gen_20 <- all_spec |>
  filter(Genus %in% gen_20)

## Calculate total counts per sample:
total_counts <- all_gen_20 |>
  group_by(Sample) |>
  summarize(Total = sum(Counts, na.rm = TRUE), .groups="drop")


## Calculate relative abundance
all_gen_20 <- all_gen_20 |>
  left_join(total_counts, by="Sample") |>
  mutate(RelAbund = Counts/Total) |>
  select(-Total)


library(dplyr)
library(ggplot2)

# Create the Location and Facet columns based on the Sample values
all_gen_20 <- all_gen_20 %>%
  mutate(Location = case_when(
    Sample %in% paste0("C", 1:4) ~ "CFIA",
    Sample %in% paste0("N", 1:28) ~ "Netherlands",
    Sample %in% paste0("NF", 1:7) ~ "Newfoundland",
    TRUE ~ "Other"
  ),
  Facet = case_when(
    Sample %in% c(paste0("C", 1:4), paste0("N", 1:10), paste0("N", 17:28)) ~ "Wart",
    Sample %in% c(paste0("N", 11:16), paste0("NF", 1:7)) ~ "Soil",
    TRUE ~ "Exclude"
  ))

# Filter out samples that should not be displayed
all_gen_20 <- all_gen_20 %>%
  filter(Facet != "Exclude")

# Arrange the Sample factor levels to organize the columns and Locations
all_gen_20 <- all_gen_20 %>%
  mutate(Sample = factor(Sample, levels = c(paste0("C", 1:4), paste0("N", 1:10), paste0("N", 17:28), paste0("N", 11:16), paste0("NF", 1:7))),
         Location = factor(Location, levels = c("CFIA", "Netherlands", "Newfoundland")))

# Visualize stacked bar plot with faceting and larger facet fonts
ggplot(all_gen_20, aes(x = Sample, y = RelAbund, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative abundance of the top 20 genera per sample",
       x = "Sample",
       y = "Relative Abundance") +
  scale_fill_manual(values = c25) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 15, face = "italic"),
        axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 20, face = "bold")) +  # Increase facet font size
  facet_grid(~ Location + Facet, scales = "free_x", space = "free_x")







## Visualize stacked bar plot:
ggplot(all_gen_20, aes(x = Sample, y = RelAbund, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Relative abundance of the top 20 genera per sample",
       x = "Sample",
       y = "Relative Abundance") +
  scale_fill_manual(values = c25) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right",
        legend.text = element_text(size=15, face="italic"),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.title = element_text(size=20, face="bold"),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 15))

#...............................................................................

test <- unique(all_gen_20$Genus) |>
  as.data.frame() |>
  rename(Genus = "unique(all_gen_20$Genus)") |>
  filter(Genus %in% pgpb$Genus)


#...............................................................................
# DETERMINING CORE MICROBIOME
#...............................................................................

## Define core microbiome per location:
core_C <- mm_all |>
  select(1:8, C1:C4) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(C1:C4) > 0) / length(c_across(C1:C4))) |>
  #filter(Ratio > 0.7) |> # Adjust based on what ratio of species present in each site you desire, do this for all areas
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus))

core_soil <- mm_soil |>
  rowwise() |>
  mutate(Ratio = sum(c_across(9:20) > 0) / length(c_across(9:20))) |>
  filter(Ratio > 0.4) |> 
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus))

core_wart <- mm_wart |>
  rowwise() |>
  mutate(Ratio = sum(c_across(9:34) > 0) / length(c_across(9:34))) |>
  filter(Ratio > 0.7) |> 
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus))

core_N_wart <- mm_wart |>
  select(1:8, starts_with("N")) |>
  mutate(Ratio = sum(c_across(9:30) > 0) / length(c_across(9:30))) |>
  filter(Ratio > 0) |> 
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus))



core_NF <- mm_all |>
  select(1:8, c(NF1:NF3, NF5:NF7)) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(c(NF1:NF3, NF5:NF7)) > 0) / length(c_across(c(NF1:NF3, NF5:NF7)))) |>
  filter(Ratio > 0) 

core_N <- mm_all |>
  select(1:8, N1:N28) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(N1:N28) > 0) / length(c_across(N1:N28))) |>
  filter(Ratio > 0.7) |>
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus)) 

core_N_soil <- mm_all |>
  select(1:8, N11:N16) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(9:14) > 0) / length(c_across(9:14))) |>
  filter(Ratio > 0) |> 
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus)) 

core_b2app <- b2_app |>
  select(1, 11:16) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(NF1:NF7) > 0) / length(c_across(NF1:NF7))) |>
  filter(Ratio > 0.7) |>
  select(1:7) |>
  as.data.frame()


core_OG <- b2 |>
  select(1:8, OG1:OG7) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(OG1:OG7) > 0) / length(c_across(OG1:OG7))) |>
  filter(Ratio > 0.05) |>
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus)) 

core_WGA <- b2 |>
  select(1:8, NF1:NF7) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(NF1:NF7) > 0) / length(c_across(NF1:NF7))) |>
  filter(Ratio > 0.05) |>
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus)) 


core_CAD <- mm_all |>
  select(1:8, C1:C4, NF1:NF7) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(9:18) > 0) / length(c_across(9:18))) |>
  filter(Ratio > 0.05) |>
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus)) 

core_NL <- mm_all |>
  select(1:8, N1:N28) |>
  rowwise() |>
  mutate(Ratio = sum(c_across(9:36) > 0) / length(c_across(9:36))) |>
  filter(Ratio > 0.05) |>
  filter(!Species %in% c("unidentified", "Other", NA),
         grepl("\\S", Species),
         grepl("\\S", Genus)) 


  

## Create a list of common species across sites:
core_spec_C <- unique(core_C$OTU)
core_spec_NF <- unique(core_NF$OTU)
core_spec_OG <- unique(core_OG$OTU)
core_spec_WGA <- unique(core_WGA$OTU)
core_spec_N <- unique(core_N$OTU)
core_spec_N_wart <- unique(core_N_wart$OTU)
core_spec_soil_N <- unique(core_N_soil$OTU)
core_asv_wart <- unique(core_wart$OTU)
core_asv_soil <- unique(core_soil$OTU)
core_spec_CAD <- unique(core_CAD$OTU) 
core_spec_NL <- unique(core_NL$OTU) # For continent

core_spec <- list(CFIA=core_spec_C, Newfoundland=core_spec_NF, Netherlands=core_spec_N)
core_type <- list(Soil=core_asv_soil, Wart=core_asv_wart)
core_treatment <- list(Original = core_spec_OG, WGA = core_spec_WGA)
core_continent <- list(Canada = core_spec_CAD, Netherlands = core_spec_NL)
core_wart <- list(Canada = core_spec_C, Netherlands = core_spec_N_wart)
core_soil_spec <- list(Canada=core_spec_NF, Netherlands=core_spec_soil_N)

## List of shared ASVs between soil and wart:
shared_type <- unique(intersect(core_spec_NF, core_spec_soil_N))

# Determine top PGPB per OTU
top_pgpb <- mm_soil |> filter(OTU %in% shared_type) |>
  select(OTU, Species, 9:20) |> # 9:30 for wart, 9:20 for soil
  filter(unique(OTU %in% shared_type)) |>
  filter(Species %in% pgpb$Species) |>
  select(-OTU) |>
  pivot_longer(cols = -Genus,
               names_to = "Sample",
               values_to = "Value") |>
  group_by(Genus) |>
  summarize(Value = sum(Value, na.rm=T)) |>
  arrange(desc(Value))


## Deterimne PGPB species:
soil_pgpb <- mm_soil |>
  filter(Genus == 'Bacillus') |>
  filter(Species %in% pgpb$Species)
print(unique(soil_pgpb$Species))

wart_pgpb <- mm_wart |>
  filter(Genus == 'Staphylococcus') |>
  filter(Species %in% pgpb$Species)
print(unique(wart_pgpb$Species))

intersect(unique(wart_pgpb$Species), unique(soil_pgpb$Species))

shared_wart <- unique(intersect(core_spec_C, core_spec_N_wart))
shared_soil <- unique(intersect(core_spec_NF, core_spec_soil_N))


# Load necessary libraries
library(dplyr)

# Load necessary libraries
library(dplyr)

# QPCR info
qpcr <- sampleinfo |>
  select(Sample, Delta_CQ) 


pgpb_wart <- mm_wart |>
  filter(OTU %in% shared_wart) |>
  filter(Species %in% pgpb$Species) |>
  filter(rowSums(across(where(is.numeric))) != 0)  |>
  select(c(1, 13:34)) |>
  tibble::column_to_rownames("OTU") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Sample")

# Merge pgpb_wart and qpcr based on the Sample column
pgpb_wart <- pgpb_wart %>%
  left_join(qpcr, by = "Sample") %>%
  mutate(Sample = ifelse(!is.na(Delta_CQ), paste0(Sample, " (", Delta_CQ, ")"), Sample)) %>%
  select(-Delta_CQ) |>  # Remove Delta_CQ column after merging 
  tibble::column_to_rownames("Sample") |>
  as.matrix()

# Ensure pgpb_wart is a data frame
pgpb_wart <- as.data.frame(pgpb_wart)

# Convert absolute counts to relative abundance
relative_abundance <- pgpb_wart %>%
  mutate(across(everything(), ~ . / sum(.)))

# Add a small pseudocount before log transformation
pseudocount <- 1e-6
relative_abundance <- relative_abundance + pseudocount

# Log-transform the data
log_transformed <- log1p(relative_abundance)  # log1p(x) is equivalent to log(x + 1)

# Ensure all columns are numeric and handle problematic values
log_transformed <- log_transformed %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ ifelse(is.nan(.), 0, .))) %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.), 0, .))) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  as.matrix()

library(pheatmap)

# Generate heatmap with clustering and scaling
pheatmap(log_transformed, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         scale = "row",  # Scale data across rows to make patterns more visible
         fontsize_number = 10,
         main = "Relative Abundance Heatmap of potential PGPB species in Netherland wart samples")






# Subset pgpb_wart dataset to include only the filtered samples
filtered_pgpb_wart <- pgpb_wart %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  filter(Sample %in% filtered_samples) %>%
  column_to_rownames("Sample")

# Ensure all columns are numeric and handle problematic values
filtered_pgpb_wart <- filtered_pgpb_wart %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ ifelse(is.infinite(.) | is.na(.), 0, .)))

# Convert absolute counts to relative abundance
relative_abundance <- filtered_pgpb_wart %>%
  rowwise() %>%
  mutate(across(everything(), ~ . / sum(c_across(everything())))) %>%
  ungroup() %>%
  as.matrix()

# Label samples with their respective Delta_CQ values
sample_labels <- qpcr %>%
  filter(Sample %in% rownames(relative_abundance)) %>%
  arrange(match(Sample, rownames(relative_abundance))) %>%
  pull(Delta_CQ)

# Generate heatmap
heatmap(relative_abundance, 
        Rowv = NA, 
        Colv = NA, 
        col = heat.colors(256), 
        scale = "row", 
        margins = c(5, 10), 
        xlab = "OTUs", 
        ylab = "Samples", 
        labRow = sample_labels)














unique(test$Species) 
test <- test |> filter(Genus %in% pgpb$Genus)







## Create a list of common genera across sites:
core_gen_C <- unique(core_C$Genus)
core_gen_NF <- unique(core_NF$Genus)
core_gen_N <- unique(core_N$Genus)
core_gen <- list(CFIA=core_gen_C, Newfoundland=core_gen_NF, Netherlands=core_gen_N)

## Create a list of common species between soil samples:
core_gen_soil_N <- unique(core_N_soil$Genus)

core_soil_gen <- list(Canada=core_gen_NF, Netherlands=core_gen_soil_N)


## Visualize coommon core species:
core_spec_venn <- ggvenn(core_spec,
                         fill_color = c("red", "blue", "green"),
                         fill_alpha = 0.35,
                         stroke_size = 1, stroke_linetype = "solid",
                         set_name_size = 10, text_size = 7.1)
core_spec_venn
# Shows only 3 species common across the core microbiome for all 3 locations
# All 3 species are baciullus
# Similarities in core microbiome between CFIA and Netherlands could be bc they're mainly wart samples instead of soil?


## For soil vs wart comparison
core_type_venn <- ggvenn(core_type,
                         fill_color = c("darkgreen", "maroon"),
                         fill_alpha = 0.35,
                         stroke_size = 1, stroke_linetype = "solid",
                         set_name_size = 15, text_size = 10)
core_type_venn


## For wart comparison between locations
ggvenn(core_wart,
       fill_color = c("maroon", "magenta"),
       fill_alpha = 0.35,
       stroke_size = 1, stroke_linetype = "solid",
       set_name_size = 15, text_size = 10)

## For continent
core_continent_venn <- ggvenn(core_continent,
                         fill_color = c("maroon", "orange"),
                         fill_alpha = 0.35,
                         stroke_size = 1, stroke_linetype = "solid",
                         set_name_size = 15, text_size = 10)
core_continent_venn



## For WGA vs OG comparison
core_treatment_venn <- ggvenn(core_treatment,
                         fill_color = c("darkgreen", "maroon"),
                         fill_alpha = 0.35,
                         stroke_size = 1, stroke_linetype = "solid",
                         set_name_size = 15, text_size = 10)
core_treatment_venn

## Visualize common core species in soil:
soil_venn_spec <- ggvenn(core_soil_spec,
                         fill_color = c("darkgreen", "cyan"),
                         fill_alpha = 0.35,
                         stroke_size = 1, stroke_linetype = "solid",
                         set_name_size = 15, text_size = 10)
soil_venn_spec # Seems to be Bacillus subtilis. Bacillus flexus, Bacillus cereus


## Visualize common core genera in soil:
soil_venn_gen <- ggvenn(core_soil_gen,
                        fill_color = c("red", "green"),
                        fill_alpha = 0.35,
                        stroke_size = 1, stroke_linetype = "solid",
                        set_name_size = 10, text_size = 7.1)
soil_venn_gen # Bacilus and Pseudomonas




## Visualize coommon core species:
core_gen_venn <- ggvenn(core_gen,
                        fill_color = c("red", "blue", "green"),
                        fill_alpha = 0.35,
                        stroke_size = 1, stroke_linetype = "solid",
                        set_name_size = 10, text_size = 7.1)
core_gen_venn





################################################################################
#                               CONTINENT COMPARISON
################################################################################

alpha_CAD <- core_CAD |>
  select(-c(2:8, Ratio)) |>
  tibble::column_to_rownames("OTU") |>
  t()

alpha_NL <- core_NL |>
  select(-c(2:8, Ratio)) |>
  tibble::column_to_rownames("OTU") |>
  t()


shannon_CAD <- exp(vegan::diversity(alpha_CAD, index = "shannon"))
simpson_CAD <- 1/(vegan::diversity(alpha_CAD, "simpson"))
chao1_CAD <- estimateR(alpha_CAD)["S.chao1", ]

shannon_NL <- exp(vegan::diversity(alpha_NL, index = "shannon"))
simpson_NL <- 1/(vegan::diversity(alpha_NL, "simpson"))
chao1_NL <- estimateR(alpha_NL)["S.chao1", ]

alpha_CAD <- data.frame(
  ShannonTD = shannon_CAD,
  Chao1 = chao1_CAD,
  SimpsonTD = simpson_CAD
) |>
  mutate(Location = "Canada")

alpha_NL <- data.frame(
  ShannonTD = shannon_NL,
  Chao1 = chao1_NL,
  SimpsonTD = simpson_NL
) |>
  mutate(Location = "Netherlands")

div_continent <- rbind(alpha_CAD, alpha_NL)


ggplot(div_continent, aes(x = Location, y = Chao1, fill = Location)) +
  geom_boxplot(alpha = 0.5, size = 1) +  # Box plot with adjusted alpha and outline thickness
  geom_jitter(position = position_jitter(width = 0.2), alpha = 1, color = "black", size = 2) +  # Increase point size
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = position_dodge(0.75), size = 1) +
  theme_classic(base_size = 16) +
  labs(title = "Chao1", x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40, face = "bold"),
        axis.text.y = element_text(size = 30)) +
  scale_fill_manual(values = c("Canada" = "maroon", "Netherlands" = "orange"))


t.test(data = div_continent, SimpsonTD ~ Location) 


## List of shared ASVs between soil and wart:
shared_type <- unique(intersect(core_spec_CAD, core_spec_NL))
test <- mm_soil |> filter(OTU %in% shared_type)
test <- test |> filter(Species %in% pgpb$Species)




################################################################################
#                               WGA VS OG CORE
################################################################################

## Calculate indicies for soil:

core_OG_tmp <- core_OG |>
  select(1, 9:15) |>
  tibble::column_to_rownames("OTU") |>
  t()

core_WGA_tmp <- core_WGA |>
  select(1, 9:15) |>
  tibble::column_to_rownames("OTU") |>
  t()


shannon_core_OG <- exp(vegan::diversity(core_OG_tmp, index = "shannon"))
shannon_core_WGA <- exp(vegan::diversity(core_WGA_tmp, index = "shannon"))
simpson_core_OG <- 1/(vegan::diversity(core_OG_tmp, "simpson"))
simpson_core_WGA <- 1/(vegan::diversity(core_WGA_tmp, "simpson"))
chao1_core_OG <- estimateR(core_OG_tmp)["S.chao1", ]
chao1_core_WGA <- estimateR(core_WGA_tmp)["S.chao1", ]


alpha_OG_core <- data.frame(
  ShannonTD = shannon_core_OG,
  Chao1 = chao1_core_OG,
  SimpsonTD = simpson_core_OG
) |>
  mutate(Treatment = "Original")

alpha_WGA_core <- data.frame(
  ShannonTD = shannon_core_WGA,
  Chao1 = chao1_core_WGA,
  SimpsonTD = simpson_core_WGA
) |>
  mutate(Treatment = "WGA")


alpha_treat_core <- rbind(alpha_OG_core, alpha_WGA_core)


# Visualize:
ggplot(alpha_treat_core, aes(x = Treatment, y = SimpsonTD, fill = Treatment)) +
  geom_boxplot(alpha = 0.5, size = 1) +  # Box plot with adjusted alpha and outline thickness
  geom_jitter(position = position_jitter(width = 0.2), alpha = 1, color = "black", size = 2) +  # Increase point size
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, fill = "white") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = position_dodge(0.75), size = 1) +
  theme_classic(base_size = 16) +
  labs(title = "SimpsonTD", x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40, face = "bold"),
        axis.text.y = element_text(size = 30)) +
  scale_fill_manual(values = c("Original" = "darkgreen", "WGA" = "maroon"))


## Determine significance:
t.test(data = alpha_treat_core, SimpsonTD ~ Treatment) 



################################################################################
#                               PGPB RECOVERY
################################################################################
all_CF <- mm |>
  select(Genus, Species, C1:C4) |>
  filter(rowSums(across(where(is.numeric)))!=0)
pgpb_spec_CF <- all_CF |>
  filter(Species %in% pgpb$Species)
pgpb_core_CF <- core_C |>
  filter(Species %in% pgpb$Species)

all_NF <- b2 |>
  select(Genus, Species, NF1:NF7) |>
  filter(rowSums(across(where(is.numeric)))!=0)
pgpb_spec_NF <- all_NF |>
  filter(Species %in% pgpb$Species)
pgpb_core_NF <- core_NF |>
  filter(Species %in% pgpb$Species)

all_NL <- mm |>
  select(Genus, Species, N1:N28) |>
  filter(rowSums(across(where(is.numeric)))!=0)
pgpb_spec_NL <- all_NL |>
  filter(Species %in% pgpb$Species)
pgpb_core_NL <- core_N |>
  filter(Species %in% pgpb$Species)

pgpb_perc_CF <- (nrow(pgpb_spec_CF) / nrow(all_CF)) * 100
pgpb_perc_NF <- (nrow(pgpb_spec_NF) / nrow(all_NF)) * 100
pgpb_perc_NL <- (nrow(pgpb_spec_NL) / nrow(all_NL)) * 100


pgpb_perc <- data.frame(
  Location = c("CF", "NF", "NL"),
  Total = c(nrow(all_CF), nrow(all_NF), nrow(all_NL)),
  PGPB = c(nrow(pgpb_spec_CF), nrow(pgpb_spec_NF), nrow(pgpb_spec_NL))
)


pgpb_bar_all <- ggplot(pgpb_perc, aes(y = reorder(Location, Total))) +
  geom_bar(aes(x = Total, fill = "Total"), stat = "identity", width = 0.5, color = "#111111") +
  geom_bar(aes(x = PGPB, fill = "PGPB"), stat = "identity", width = 0.5, color = "#111111") +
  geom_text(aes(x = Total, label = Total), vjust = -3.5, color = "black", size = 7) +
  geom_text(aes(x = PGPB, label = PGPB), vjust = -3.5, color = "black", size = 7) +
  theme_classic(base_size = 16) +
  labs(title = "Number of total species vs PGPB species from each Location",
       x = "Number of Species",
       y = "Location",
       fill = "Species Type") +
  scale_y_discrete(labels = c("NF" = "Newfoundland", "CF" = "CFIA", "NL" = "Netherlands")) +
  scale_fill_manual(values = c("Total" = "#58a36c", "PGPB" = "#e3d044")) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 30, face="bold"),
        axis.title.y = element_text(size = 30, face="bold"),
        axis.text.y = element_text(size = 20, face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))
pgpb_bar_all




pgpb_perc <- data.frame(
  Location = c("CF", "NF", "NL"),
  Total = c(nrow(core_C), nrow(core_NF), nrow(core_N)),
  PGPB = c(nrow(pgpb_core_CF), nrow(pgpb_core_NF), nrow(pgpb_core_NL))
)


pgpb_bar_core <- ggplot(pgpb_perc, aes(y = reorder(Location, Total))) +
  geom_bar(aes(x = Total, fill = "Total"), stat = "identity", width = 0.5, color = "#111111") +
  geom_bar(aes(x = PGPB, fill = "PGPB"), stat = "identity", width = 0.5, color = "#111111") +
  geom_text(aes(x = Total, label = Total), vjust = -3.5, color = "black", size = 7) +
  geom_text(aes(x = PGPB, label = PGPB), vjust = -3.5, color = "black", size = 7) +
  theme_classic(base_size = 16) +
  labs(title = "Number of PGPB species found in core microbiomes",
       x = "Number of Species",
       y = "Location",
       fill = "Species Type") +
  scale_y_discrete(labels = c("NF" = "Newfoundland", "CF" = "CFIA", "NL" = "Netherlands")) +
  scale_fill_manual(values = c("Total" = "#58a36c", "PGPB" = "#e3d044")) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 30, face = "bold"),
        axis.title.x = element_text(size = 30, face="bold"),
        axis.title.y = element_text(size = 30, face="bold"),
        axis.text.y = element_text(size = 20, face="bold"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20))
pgpb_bar_core



################################################################################
#                         PGPB ABUNDANCE LEVEL TO SENDO
################################################################################













################################################################################
#                               BETA DIVERSITY ANALYSIS
################################################################################
beta_CF <- core_C |>
  select(OTU, C1:C4) |>
  pivot_longer(cols = -OTU,
               names_to = "Sample",
               values_to = "Value")

beta_NF <- core_b2app |>
  select(OTU, NF1:NF7) |>
  pivot_longer(cols = -OTU,
               names_to = "Sample",
               values_to = "Value")

beta_NL <- core_N |>
  select(OTU, N1:N28) |>
  pivot_longer(cols = -OTU,
               names_to = "Sample",
               values_to = "Value")

beta_mm <- rbind(beta_CF, beta_NF, beta_NL) |>
  pivot_wider(names_from = "Sample",
              values_from = "Value") |>
  mutate(across(everything(), ~ replace_na(.x, 0))) |>
  tibble::column_to_rownames("OTU") |>
  t() |>
  as.matrix()


## Do permanova:
# Distance matrix:
beta_mm_dist <- vegdist(beta_mm, method = "robust.aitchison")

# Define the groups (locations) as a factor
groups <- factor(c(rep("CFIA", 4), rep("Newfoundland1", 6), rep("Newfoundland2", 28)))


# Run the PERMANOVA test
result <- adonis2(beta_mm_dist ~ groups)
# Results indicate highly significant difference (0.001). R^2 suggests 46% of the variation is attributed to the difference between locations.

# View the results
print(result)



set.seed(11122024)
nmds_all <- metaMDS(beta_mm, distance = "robust.aitchison", k = 2, trymax = 100)

# Extract NMDS site scores and add sample identifiers
nmds_scores_all <- as.data.frame(scores(nmds_all)$sites) |>
  tibble::rownames_to_column("Sample") |>
  mutate(group = case_when( startsWith(Sample, "C") ~ "CFIA", 
                            startsWith(Sample, "NF") ~ "Newfoundland", 
                            startsWith(Sample, "N") ~ "Netherlands", 
                            TRUE ~ "Other" ))


# Plot the NMDS with the updated group names and stress value

# Calculate the stress value if it's not already in your nmds object
stress_value <- nmds_all$stress  # Replace 'nmds' with your NMDS object

ggplot(nmds_scores_all, aes(x = NMDS1, y = NMDS2, color = group, shape = group, label = Sample, group = group)) +
  geom_point(size = 4) +  # Increase point size
  stat_ellipse(aes(fill = group), type = "norm", alpha = 0.5, geom = "polygon", size = 5) +
  geom_text_repel(size = 5, max.overlaps = Inf) +  # Increase text size and allow for all overlaps
  theme_classic(base_size = 16) +
  labs(title = "NMDS plot comparing core microbiomes",
       x = "NMDS1",
       y = "NMDS2",
       color = "Location",
       shape = "Location",
       fill = "Location") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold")) +
  annotate("text", x = Inf, y = Inf, label = paste("Stress =", round(stress_value, 3)), 
           hjust = 1, vjust = 1, size = 15, color = "black", fontface = "bold")



################################################################################
#                               NETWORK ANALYSIS
################################################################################

## Read data:
WGA_nodes <- read_delim("data/MENA/WGA_node.tsv", delim="\t")
OG_nodes <- read_delim("data/MENA/OG_node.tsv", delim="\t")


## Determine which taxa are unique to each treatment
shared_nodes <- intersect(WGA_nodes$Name, OG_nodes$Name) |>
  as.data.frame() |>
  rename(Name = "intersect(WGA_nodes$Name, OG_nodes$Name)")

## Determine unique nodes for both groups:
unique_WGA <- setdiff(WGA_nodes$Name, OG_nodes$Name) 
unique_OG <- setdiff(OG_nodes$Name, WGA_nodes$Name)

## Add a column called 'Type' indicating shared/unique for those groups:
WGA_nodes$Type <- ifelse(WGA_nodes$Name %in% shared_nodes$Name, "Shared", "Unique")
OG_nodes$Type <- ifelse(OG_nodes$Name %in% shared_nodes$Name, "Shared", "Unique")

## Merge taxonomy information:
b2_tax <- b2 |>
  select(1:8)

WGA_nodes <- merge(WGA_nodes, b2_tax, by.x = "Name", by.y = "OTU", all.x = TRUE)
OG_nodes <- merge(OG_nodes, b2_tax, by.x = "Name", by.y = "OTU", all.x = TRUE)


## Create new column for PGPB:
WGA_nodes <- WGA_nodes |>
  mutate(PGPB = ifelse(Genus %in% pgpb$Genus, "PGPB", "Non-PGPB"))

OG_nodes <- OG_nodes |>
  mutate(PGPB = ifelse(Genus %in% pgpb$Genus, "PGPB", "Non-PGPB"))


write_csv(WGA_nodes, 'data/MENA/WGA_node_master.csv')
write_csv(OG_nodes, 'data/MENA/OG_node_master.csv')



