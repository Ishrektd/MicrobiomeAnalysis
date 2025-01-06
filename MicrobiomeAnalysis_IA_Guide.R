#-------------------------------------------------------------------------------
#                                 Ishraq Akbar
#                             ishraqakbar@gmail.com
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
####                            INTRODUCTION
#-------------------------------------------------------------------------------

# Welcome to my EDA project! In this script, we will analyze 16S metabarcoding information to explore bacterial community compositions in potato-wart infested soil samples. Sequencing of the data was done using the Emu pipeline, you may find more information about it on their github (https://gitlab.com/treangenlab/emu) if you're interested.

# In this project, I hope to discover genera of bacteria that are extensively studied as Biocontrol Agents (BCAs) in hopes of discovering unique species that could serve as potential defense organisms against Synchytrium endobioticum - a fungal pathogen responsible for inducing potato wart disease. In doing so, it can pave the way to additional research to investigate gene expression pathways associated with those given species to determine if they can confer resistances against the plant pathogen. I will also be assessing the difference in bacterial community composition between original 16S genomic DNA, and whole genome amplified (WGA) 16S genomic DNA to assess if laboratory preparation can affect the sequencing steps.

# Here is an outline of the rest of my project:


# Part 1
# I will install and import all the necessary packages, and I will import my data and define functions that will be used later on in the data tidying process.


# Part 2 - Data Tidying
# I will attempt to organize and filter my dataframes so that I am able to work with the functions used in this analysis.

# Part 3 - Visualization
# I will attempt to visualize my dataframes after data tidying and try to make sense of it using various visualization methods such as NMDS plots, barplots, PCoA plots, etc.

# Part 4 - Discussion
# I will wrap up my findings of this EDA in a concluding paragraph



#-------------------------------------------------------------------------------
####                         PART 1 - DATA IMPORT
#-------------------------------------------------------------------------------
# Install following libraries if you haven't.
# Commented out to avoid breaking code:

# install.packages("tidyverse")
# install.packages("vegan")
# install.packages("ggpubr")
# install.packages("tidytext")
# install.packages("forcats")


# Libraries used:
library(tidyverse)
library(vegan)
library(ggpubr)
library(tidytext)
library(forcats)


# Read files
mm <- read_csv('data/batch2.csv') # Full meta data on relative abundance
mock <- read_csv('data/mock-community-data.csv') # Mock community data on relative abundance
b2 <- read_csv('data/batch2.csv') # Full meta data on relative abundance, will be tidied
b2_abs <- read_csv('data/emu-combined-tax_id-counts.csv') # Meta data on absolute counts


# Functions used in this EDA:

# Create Sorter function to organize taxonomic ranks in order of highest to lowest after tax_id:
col_order <- function(df) {
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
#...............................................................................



#-------------------------------------------------------------------------------
####                         PART 2 - DATA TIDYING
#-------------------------------------------------------------------------------



#...............................................................................
# Seeing if the mock community positive control has the expected species:

# Filtering for the top 10 abundant species for the positive control (barcode15.fastq), and arranging them in descending order:
mm <- mm |>
  select(tax_id, species, barcode15.fastq) |>
  filter(barcode15.fastq > 0) |>
  mutate(barcode15.fastq = barcode15.fastq * 100) |>
  arrange(desc(barcode15.fastq)) |>
  slice(1:10)

# Seeing which species are common between the expected mock community list and our positive controls, and creating a list of species that are common
common_species <- inner_join(mm, mock, by = "species") |>
  select(species) |>
  distinct()

# Overrides the mm dataframe to combine the common species list with the values that we received:
mm <- mm |>
  filter(species %in% common_species$species) |>
  select(species, barcode15.fastq) |>
  arrange(species)

# Combining the mock and mm dfs with inner_join() so that we can visualize the similarities and differences using ggplot():

# First combine dataframes into a new dataframe called data_mock:

# To do this, filter mock data down to the same species as in the data df:
mock <- mock |>
  filter(species %in% common_species$species) |>
  select(species, "16S_only") |>
  arrange(species)

# Now bind the dataframes together:
data_mock <- inner_join(mm, mock, by = "species")
data_mock # Seeing if it worked

# Reshape the data to a long format for visualization:
data_mock <- data_mock |>
  pivot_longer(cols = c(barcode15.fastq, "16S_only"), names_to = "sample", values_to = "percentage")

# Changing the column values for visualization:
data_mock$sample[data_mock$sample == "barcode15.fastq"] <- "actual abundance"
data_mock$sample[data_mock$sample == "16S_only"] <- "theoretical abundance"
#...............................................................................




#...............................................................................

# Quickly insepect number of columns and rows:
tibble(b2) # 3135 rows x 37 columns
tibble(b2_abs) # 3135 rows x 37 columns


# Run the functions on the dataframes:
b2 <- col_order(b2)
b2 <- num_sorter(b2)
b2_abs <- col_order(b2_abs)
b2_abs <- num_sorter(b2_abs)

# Inspect:
b2
b2_abs

# Now we want to remove the columns containing relative abundance information and absolute counts for only thresholds over 0.0001 because we are interested in analyzing all the samples. As such, we'll go ahead and remove those columns. We'll also remove the barcode15 & barcode16 columns as they represent positive and negative controls respectively. Because the goal of this project is to compare between the original treatments and WGA treatments, it is not necessary to keep them:

# Note that the first 8 columns contain taxonomy information, so we'll select those columns and columns ending with fastq:
b2 <- b2 |>
  select(1:8, ends_with("fastq"), -barcode15.fastq, -barcode16.fastq)

b2_abs <- b2_abs |>
  select(1:8, ends_with("fastq"), -barcode15.fastq, -barcode16.fastq)


# Now let's rename these columns to something that's easier to understand. barcodes01-07 represent 7 soil samples with original DNA, whereas barcodes08-14 represent the whole genome amplified (WGA) treatment of those original DNA samples to increase the available DNA concentration to work with during library preparation. barcode15 & barcode16 represent positive and negative samples respectively. We are interested in visualizing the microbial community characteristics between these two treatments using various models.

# Rename the columns to something easier to understand:
b2 <- b2 |>
  rename(
    og_rel_1 = barcode01.fastq, og_rel_2 = barcode02.fastq, og_rel_3 = barcode03.fastq,
    og_rel_4 = barcode04.fastq, og_rel_5 = barcode05.fastq, og_rel_6 = barcode06.fastq,
    og_rel_7 = barcode07.fastq, wga_rel_1 = barcode08.fastq, wga_rel_2 = barcode09.fastq,
    wga_rel_3 = barcode10.fastq, wga_rel_4 = barcode11.fastq, wga_rel_5 = barcode12.fastq,
    wga_rel_6 = barcode13.fastq, wga_rel_7 = barcode14.fastq
  )

b2_abs <- b2_abs |>
  rename(
    og_abs_1 = barcode01.fastq, og_abs_2 = barcode02.fastq, og_abs_3 = barcode03.fastq, 
    og_abs_4 = barcode04.fastq, og_abs_5 = barcode05.fastq, og_abs_6 = barcode06.fastq,
    og_abs_7 = barcode07.fastq, wga_abs_1 = barcode08.fastq, wga_abs_2 = barcode09.fastq,
    wga_abs_3 = barcode10.fastq, wga_abs_4 = barcode11.fastq, wga_abs_5 = barcode12.fastq,
    wga_abs_6 = barcode13.fastq, wga_abs_7 = barcode14.fastq
  )

# See if it worked:
colnames(b2)
colnames(b2_abs)
#...............................................................................





#...............................................................................

# Now it will be better if we could organize the dataframes such that one dataframe will contain all the metadata information on original samples, and the other will contain all the metadata information on wga samples. Let's go ahead and do this:

# Select the taxonomy and 'og' columns from both b2 and b2_abs dataframes:
b2_og <- select(b2, tax_id:species, starts_with("og_rel"))
og_b2_abs <- select(b2_abs, starts_with("og_abs"))

# Combine the selected columns into a new dataframe
og <- bind_cols(b2_og, og_b2_abs)

# Select the taxonomy and 'wga' columns from both b2 and b2_abs dataframes
b2_wga <- select(b2, tax_id:species, starts_with("wga_rel"))
wga_b2_abs <- select(b2_abs, starts_with("wga_abs"))

# Combine the selected columns into a new dataframe
wga <- bind_cols(b2_wga, wga_b2_abs)

# Check to see if it worked:
tibble(og)
tibble(wga)


# Let's create a new column called taxonomy for the b2 column. This column will be a concatenated format of the other taxonomy columns in a way such that it will be easy to extract taxonomy information if needed. We will do this by combining the taxonomy columns into one string where each hierarchical classification is separated the following way: k__x; p__x; c__x; o__x; f__x; g__x; s__x; Note that x represents the hierarchical classification for each letter (k = superkingdom, p = phylum, c = class, o = order, f = family, s = species):

b2 <- b2 |>
  mutate(taxonomy = paste0("k__", superkingdom, "; ",
                           "p__", phylum, "; ",
                           "c__", class, "; ",
                           "o__", order, "; ",
                           "f__", family, "; ",
                           "g__", genus, "; ",
                           "s__", species))

b2_abs <- b2_abs |>
  mutate(taxonomy = paste0("k__", superkingdom, "; ",
                           "p__", phylum, "; ",
                           "c__", class, "; ",
                           "o__", order, "; ",
                           "f__", family, "; ",
                           "g__", genus, "; ",
                           "s__", species))
# See if it worked:
head(b2$taxonomy)
head(b2_abs$taxonomy)
#...............................................................................



#...............................................................................
# After that, we will separate the asv dataframe columns for further analysis. We will create a taxonomy dataframe only containing taxonomy information for each tax_id, and an asv dataframe only containing the absolute counts for each tax_id. This will be for the b2_abs dataframe as that is what we will largely be doing our asv analyses with to compare between treatment groups.

# Start with taxonomy dataframes:
b2_tax <- b2 |>
  select(tax_id, 2:8, taxonomy) |>
  relocate(taxonomy, .after = tax_id)

# Now make the ASV tables. We will only filter for the tax_id column, and the sample columns. We will also convert any NA values to 0:
b2_asv <- b2 |>
  select(tax_id, starts_with(c("og", "wga"))) |>
  mutate(across(starts_with(c("og", "wga")), ~replace_na(., 0))) |>
  mutate(across(starts_with(c("og", "wga")), ~ . * 100)) # Multiplying by 100 here because the values are currently a fraction and total to 1 for each sample


# Let's quickly check to see if all the sample columns sum to 100:
all(round(colSums(b2_asv[, 2:ncol(b2_asv)])) == 100)


# Do the same for the absolute count data. Also note that for the analysis steps, functions such as rarefy will only accept integer count values, so while we're at it, let's convert the values to whole numbers as well.
b2_abs_asv <- b2_abs |>
  select(tax_id, starts_with(c("og", "wga"))) |>
  mutate(across(starts_with(c("og", "wga")), ~replace_na(., 0))) |>
  mutate(across(starts_with(c("og", "wga")), ~as.integer(round(.))))


# Check to see if it worked:
head(b2_abs_asv)

#...............................................................................



#-------------------------------------------------------------------------------
####                         PART 3 - VISUALIZATION
#-------------------------------------------------------------------------------

# Let's attempt to visualize the mock community data and compare it to actual seq data:
ggplot(data_mock, aes(x = species, y = percentage, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  
  labs(title = "Figure 1. Comparing the Relative Abundance Composition of Positive Controls to Theoretical ZymoBIOMICS™ Microbial Community DNA Standards",
       x = "Species", 
       y = "Percentage (%)", 
       fill = "Sample") +
  scale_fill_discrete(labels = c("Theoretical values", "Sequencing values")) +
  scale_fill_manual(values = c("steelblue", "maroon")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.title = element_text(size=15),
        axis.text.x.bottom = element_text(face="italic"))

# Shows results are similar and fall within 15% of each other which is good
# Indicates sequencing quality is good and can now proceed to next steps

#...............................................................................




#...............................................................................

# Now we will attempt to rarefy our data with the vegan R package. To first do so, we need to convert our dataframe containing absolute counts into wide format, and make our columns are tax_id while the rows will be our samples. Let's go ahead and do this:


# First temporarily convert the dataframe to long format
b2_abs_wide <- pivot_longer(b2_abs_asv, cols = -tax_id, names_to = "sample", 
                            values_to = "value") 

# Revert the dataframe back to wide format
b2_abs_wide <- pivot_wider(b2_abs_wide, 
                           names_from = tax_id, 
                           values_from = value) |> as.data.frame()

# Assigning the row names as the sample:
rownames(b2_abs_wide) <- b2_abs_wide$sample


# Removing the column containing sample information because we need this in matrix format to be able to work with it. But before we do that, let's pull the names of the sample columns as it will help us with subsequent analysis. We'll also create a new column titled group which indicates whether it is og or wga for the treatment type for a given sample ID:

b2_names <- b2_abs_wide |>
  select(sample) |>
  mutate(group = case_when(
    startsWith(sample, "og") ~ "og",
    startsWith(sample, "wga") ~ "wga"
  ))



# Check to see if it worked:
b2_names

# Now let's go ahead and remove the sample column
b2_abs_wide <- b2_abs_wide[,-1]
rownames(b2_abs_wide)[0] <- "samples"


# Before we rarefy our data, we need to see what the lowest number of counts is for each sample. Let's go ahead and see that:
raremax <- min(rowSums(b2_abs_wide))
raremax # Shows the lowest number of counts is 95,460

# Now let's visualize our rarefy data. The initial output by the vegan package will look messy, so we'll extract the output from the rarecurve() function and visualize it in a neater way using ggplot:

# First run the rarecurve() function and store the output into a variable called rarecurve_data
rarecurve_data <- rarecurve(b2_abs_wide, step = 100, sample = raremax)

# This is taking each element of the list and binds them together to make a dataframe
map_dfr(rarecurve_data, bind_rows) |>
  bind_cols(sample = rownames(b2_abs_wide),) |>
  pivot_longer(-sample) |>
  drop_na() |>
  mutate(counts = as.numeric(str_replace(name, "N",""))) |>
  
  # Now plotting the graph:
  ggplot(aes(x=counts, y=value, group=sample, color=sample)) +
  geom_vline(xintercept = raremax, color ="gray") +
  geom_line(size=0.8) +
  labs(x="Sequencing read counts", y="Number of observed species", 
       title="Figure 2. Rarefaction visualization of sequencing samples.") +
  scale_color_discrete(name = "Sample ID") +
  theme_classic2() +
  scale_x_continuous(labels = scales::comma)  


# Based on the results, all our samples plateau as the sequencing read counts approach a constant number of species for each sample, thereby indicating there was enough sequencing to cover all species within the samples. The rarefraction curves allow us to quickly visualize the number of observed species recorded for the number of sequences. 

#...............................................................................





#...............................................................................

# For data analysis to assess community composition, we'll go ahead and use the og_b2 and wga_b2 dataframes.

# Reshape the dataframes into long format
og_long <- b2_og |>
  pivot_longer(cols = starts_with("og_rel"), names_to = "sample", values_to = "value") %>%
  drop_na(value) |>
  group_by(sample) |>
  mutate(percentage = value / sum(value) * 100) |>
  as.data.frame()

wga_long <- b2_wga |>
  pivot_longer(cols = starts_with("wga_rel"), names_to = "sample", values_to = "value") %>%
  drop_na(value) |>
  group_by(sample) |>
  mutate(percentage = value / sum(value) * 100) |>
  as.data.frame()


# Now let's visualize the most abundant phyla within our groups. To do so, we'll first create new dataframes titled og_phy and wga_phy to store our information so that we'll be able to compare them side by side using ggarrange().We are summing the total relative abundance for each phylum across each sample, and then multiplying it by 100 to get a percentage as the relative abundance information is stored as a fraction. We are then replacing any NA values as Unclassified, and plotting our results using ggplot():

# Tidying data to select phyla and replacing NA values with Unclassified:
og_phy <- og_long |>
  group_by(sample, tax_id, phylum, value) |>
  summarize(value = sum(value), .groups="drop") |>
  group_by(sample, phylum) |>
  summarize(value = value*100, .groups="drop") |>
  mutate(phylum = replace_na(phylum, "Unclassified")) 

wga_phy <- wga_long |>
  group_by(sample, tax_id, phylum, value) |>
  summarize(value = sum(value), .groups="drop") |>
  group_by(sample, phylum) |>
  summarize(value = value*100, .groups="drop") |>
  mutate(phylum = replace_na(phylum, "Unclassified"))

# Reordering levels of phylum before to descending order before creating the plots
og_phy$phylum <- fct_reorder(og_phy$phylum, og_phy$value, .desc = TRUE)
wga_phy$phylum <- fct_reorder(wga_phy$phylum, wga_phy$value, .desc = TRUE)


# Code for plots:
og_phy_bar <- ggplot(data=og_phy, aes(x=sample, y=value, fill=phylum)) +
  geom_col() +
  labs(y="Percentage (%)", x="Sample ID", fill="Phyla", title = "Figure 3A. Visualizing most abundant phyla for original DNA") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8, face="italic"))

wga_phy_bar <- ggplot(data=wga_phy, aes(x=sample, y=value, fill=phylum)) +
  geom_col() +
  labs(y="Percentage (%)", x="Sample ID", fill="Phyla", 
       title = "Figure 3B. Visualizing most abundant phyla for WGA treatment") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8, face="italic"))

# Visualize:

# Visualize:
ggarrange(og_phy_bar, wga_phy_bar)

# Based on this, we can see that the composition in phyla clearly differ between the original and WGA treatments. The most abundant treatment for the original samples appears to be Chloroflexi, and for the WGA samples it appears to be Firmicutes. Regardless, both are common BCAs that have been determined in soil samples based on past studies (https://www.jstage.jst.go.jp/article/jsme2/32/4/32_ME17131/_article). This is a good sign based on this preliminary data.

#...............................................................................





#...............................................................................

# Now let's do the same but with the genera. However note that if trying to visualize the percentage of genera across all samples, there will be a lot. As such, we would need to further filter our data down to only visualize genera across a certain threshold. Let's go ahead and do that


# Creating temporary variables to tmp1 and tmp2 to store filtered og and wga data respectively containing information on the relative abundance percentages of each genus per sample:

tmp1 <- og_long |>
  group_by(sample, tax_id, genus, value) |>
  summarise(value = sum(value), .groups="drop") |>
  group_by(sample, genus) |>
  summarise(value = value*100, .groups="drop") |>
  mutate(genus = replace_na(genus, "Unclassified"))

tmp2 <- wga_long |>
  group_by(sample, tax_id, genus, value) |>
  summarise(value = sum(value), .groups="drop") |>
  group_by(sample, genus) |>
  summarise(value = value*100, .groups="drop") |>
  mutate(genus = replace_na(genus, "Unclassified"))


# Before we pool for rare taxa, let's first quickly assess the top 10 genera percentages for each of the treatment groups, and use that as the respective cut-off points:


# Shows the lowest percentage would be 2.29%, we'll use that as our cut-off point for the og samples:
tmp1_min <- tmp1 |>
  group_by(genus) |>
  summarize(max = max(value)) |>
  arrange(desc(max)) |>
  slice(10) |>
  pull(max)

# Shows the lowest percentage would be 6.46%, we'll use that as our cut-off point for the wga samples:
tmp2_min <- tmp2 |>
  group_by(genus) |>
  summarize(max = max(value)) |>
  arrange(desc(max)) |>
  slice(10) |>
  pull(max)

# Pooling rare taxa for percentages <= 2% and storing them into respective variables for og and wga:
og_pool <- tmp1|>
  group_by(genus) |>
  summarize(pool = max(value) <= tmp1_min, .groups="drop") 

wga_pool <- tmp2 |>
  group_by(genus) |>
  summarize(pool = max(value) <= tmp2_min, .groups="drop") 


# Now visualizing the dataframes to filter by genus <= 2% for both the original and wga treatments:

# Combining the pooled data information with the temporary dataframes, and then graphing the results using ggplot():
og_gen <- inner_join(tmp1, og_pool, by="genus") |>
  mutate(genus = if_else(pool, "Other", genus)) |>
  group_by(sample, genus) |>
  summarize(value = sum(value))


og_gen_bar <- og_gen |>
  mutate(genus = fct_reorder(genus, value, .desc = TRUE)) |>
  ggplot(aes(x=sample, y=value, fill=genus)) +
  geom_col() +
  labs(y="Percentage (%)", x="Sample ID", fill="Genera",
       title = "Figure 3A. Visualizing the percentage of genera present across all original treatment samples") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8, face="italic"))

wga_gen <- inner_join(tmp2, wga_pool, by="genus") |>
  mutate(genus = if_else(pool, "Other", genus)) |>
  group_by(sample, genus) |>
  summarize(value = sum(value)) 

wga_gen_bar <- wga_gen |>
  mutate(genus = fct_reorder(genus, value, .desc=TRUE)) |>
  ggplot(aes(x=sample, y=value, fill=genus)) +
  geom_col() +
  labs(y="Percentage (%)", x="Sample ID", fill="Genera", 
       title = "Figure 3B. Visualizing the percentage of genera present across all wga treatment samples") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8, face="italic"))


# Visualize:
ggarrange(og_gen_bar, wga_gen_bar)

# Based on this, we see that the most abundant genera for both groups falls within the other category, which somewhat makes sense since there are a lot of genera in our sample. Regardless, we see that genera such as Bacillus are present, which has been noted to have BCA properties against other pathogens in the past as well (https://www.jstage.jst.go.jp/article/jsme2/32/4/32_ME17131/_article). 

#...............................................................................



#...............................................................................

# Now let's go ahead and determine if there is a significant difference between the original treatment and wga treatment groups. To do this, we will use various functions within the vegan R package. We will attempt to visualize our data using MetaMDS, PERMANOVA, and dispersion plots.

# Normalize the data using decostand(). There are multiple methods to choose from, but I have decided to go with hellinger transformation. 

b2_norm <- decostand(b2_abs_wide, method = "hell")


# Run the NMDS model to visualize composition data
set.seed(100)
nmds_results <- metaMDS(b2_norm, distance="robust.aitchison")

# Extract the NMDS scores
nmds_scores <- as.data.frame(scores(nmds_results)$sites)

# Determine the centroids for the NMDS graph:
group_centroids <- data.frame(
  Treatment = c("og", "wga"),
  Centroid_X = c(mean(nmds_scores$NMDS1[b2_names$group == "og"]),
                 mean(nmds_scores$NMDS1[b2_names$group == "wga"])),
  Centroid_Y = c(mean(nmds_scores$NMDS2[b2_names$group == "og"]),
                 mean(nmds_scores$NMDS2[b2_names$group == "wga"]))
)

# Now create a dataframe containing the nmds data that we want to visualize with ggplot:
nmds_data <- data.frame(
  Treatment = b2_names$group,
  NMDS1 = nmds_scores$NMDS1,  
  NMDS2 = nmds_scores$NMDS2,  
  xend = c(rep(group_centroids[1,2], 1), rep(group_centroids[2,2], 1)),
  yend = c(rep(group_centroids[1,3], 1), rep(group_centroids[2,3], 1))
)

# Now plot the data:
ggplot(nmds_data, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = Treatment), size=2) + 
  stat_ellipse(geom = "polygon", alpha = 0.04, aes(group = Treatment),
               colour = "black", fill = "steelblue") +
  geom_point(data = group_centroids, aes(x = Centroid_X, y = Centroid_Y),
             colour = "black", size = 2, shape = 7) +
  geom_segment(data = nmds_data, aes(x = NMDS1, y = NMDS2, xend = xend, 
                                     yend = yend, colour = Treatment), alpha = 0.5) +
  scale_colour_manual(name = "Treatment type", labels = unique(nmds_data$Treatment),
                      values = c("cyan", "magenta")) +
  labs(title = "Figure 5. NMDS plot comparing the two soil sample treatment groups") +
  theme_classic()
#...............................................................................




#...............................................................................

# Now let's do PERMANOVA to assess whether there are significant differences between the original and wga treatment groups or not


# Calculating dissimilarity indices:
perm_dist <- vegdist(b2_norm, method = "robust.aitchison")
perm_dist

# Let's now check to see if the two treatment groups have similar statistically similar multivariate dispersion. 

dispersion <- betadisper(perm_dist, group = b2_names$group, type = "centroid")
plot(dispersion, main = "Figure 6. Multivariate dispersion plot")
# Visually, we can see the groups are distinct from each other.


perm_results <- adonis2(perm_dist~as.factor(nmds_data$Treatment), data=perm_dist,
                        permutations=9999)
perm_results

# Running a PERMANOVA test, we can see that the p-value as indicated by Pr(>F) is 0.0004, which indicates that the difference between the two treatments is very significant.

#...............................................................................




#-------------------------------------------------------------------------------
####                         PART 4 - DISCUSSION
#-------------------------------------------------------------------------------

# In conclusion, this EDA has shown me that there is a significant difference between the original and WGA treatments when it comes to bacterial community composition. This suggests to me that if possible, I should avoid using WGA DNA for my library preparation as it could lead to biases and amplification of different species as highlighted in Figures 4 and 5. However, I did find bacterial taxa within both treatments that have been previously found to have BCA properties. As such, these preliminary results seem promising in suggesting that upon pathogen infection to S. endobioticum, soil microbiomes may undergo a shift in their community composition to favour bacterial taxa having BCA properties to protect them from environmental stressors. However, more studies are required to validate this.

#...............................................................................

