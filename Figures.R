# data viz for all-hands meeting: 

# Intro: map of all sequences, map of E-V159A, map of NS4A-A85T.

# overview of mutational landscape: genome with EVERY AA mutation. could just be a big barplot
# with a bar for each AA in orer, marked up in Designer to show gene locations and key mutations
# barplot of all mutations occuring in more than 1% of genomes, all in grey except E-V159A and NS4A-A85T
# to show how NS4A-A85T may not be particularly unique in its ability to gain prominence
#mutation co-occurance maps for NS4A-A85T and NS5-K314R, NS2A-R188K and NS4B-I240M, state level

# rest of presentation is walking through questions to answer: are there regional clusters? 
# what are there ranges and what climactic or environmental variables predict distribution?
# When and where did each of the mutations come from? 
#

library(choroplethr)
library(choroplethrMaps)
library(ggplot2)
library(RColorBrewer)
#import and process data
WNV_metadata <- read.csv("data_summary/Summary_13_feb_2024.csv", header = TRUE)
accessions_of_good_genomes <- read.csv("data_summary/AA_polymorphisms_sparse_filtered.csv", header = TRUE)
filtered_metadata <- merge(WNV_metadata, accessions_of_good_genomes, by = "accession")
filtered_metadata <- filtered_metadata[filtered_metadata$county != "", ]
filtered_metadata <- filtered_metadata[filtered_metadata$FIPS != "", ]

#PLOT 1: US distribution of all high quality genomes ----

my.cols <- brewer.pal(6, "OrRd")

# Count the occurrences of each county_FIPS value
counts <- table(filtered_metadata$county_FIPS)

# Convert the result to a dataframe
us_distribution_df <- as.data.frame(counts)

# Rename the columns
names(us_distribution_df) <- c("region", "value")
us_distribution_df$region <- as.numeric(as.character(us_distribution_df$region))
all_genomes_map <- choroplethr::county_choropleth(us_distribution_df, title = "Distribution of WNV genomes", legend = "# of genomes")
all_genomes_map + scale_fill_manual(values = my.cols, na.value = "white")

#PLOT 2: distribution of NW02 ----
my.cols <- brewer.pal(6, "OrRd")
E.V159A_data <- aggregate(E.V159A ~ county_FIPS, data = filtered_metadata, sum)
names(E.V159A_data) <- c("region", "value")
E.V159A_data$region <- as.numeric(as.character(E.V159A_data$region))
E.V159A_map <- choroplethr::county_choropleth(E.V159A_data, title = "Distribution of WN02", legend = "Presence of E-V159A")
E.V159A_map + scale_fill_manual(values = my.cols, na.value = "white")

#PLOT 3: distribution of SW03 ----
my.cols <- brewer.pal(6, "OrRd")
my.cols[1] <- '#ffffff'
NS4A.A85T_data <- aggregate(NS4A.A85T ~ county_FIPS, data = filtered_metadata, sum)
names(NS4A.A85T_data) <- c("region", "value")
NS4A.A85T_data$region <- as.numeric(as.character(NS4A.A85T_data$region))
NS4A.A85T_map <- choroplethr::county_choropleth(NS4A.A85T_data, title = "Distribution of NS4A.A85T", legend = "Presence of NS4A-A85T")
NS4A.A85T_map + scale_fill_manual(values = my.cols, na.value = "white")

#PLOT 4: distribution of NS5-K314R ----
my.cols <- brewer.pal(4, "OrRd")
my.cols[1] <- '#ffffff'
NS5.K314R_data <- aggregate(NS5.K314R ~ county_FIPS, data = filtered_metadata, sum)
names(NS5.K314R_data) <- c("region", "value")
NS5.K314R_data$region <- as.numeric(as.character(NS5.K314R_data$region))
NS5.K314R_map <- choroplethr::county_choropleth(NS5.K314R_data, title = "Distribution of NS5.K314R", legend = "Presence of NS5.K314R")
NS5.K314R_map + scale_fill_manual(values = my.cols, na.value = "white")
#PLOT 5: distribution of NS5-K314R AND NS4A.A85T together ----

# Create a new column that indicates whether both NS5.K314R and NS4A.A85T are present
filtered_metadata$NS5.K314R_NS4A.A85T <- ifelse(filtered_metadata$NS5.K314R == 1 & filtered_metadata$NS4A.A85T == 1, 1, 0)

# Aggregate this new column by county_FIPS
both_present_data <- aggregate(NS5.K314R_NS4A.A85T ~ county_FIPS, data = filtered_metadata, sum)

# Rename the columns and convert the region column to numeric
names(both_present_data) <- c("region", "value")
both_present_data$region <- as.numeric(as.character(both_present_data$region))

NS5.K314R_NS4A.A85T_map <- choroplethr::county_choropleth(both_present_data, title = "Distribution of NS5.K314R and NS4A.A85T", legend = "Presence of both NS5.K314R and NS4A.A85T")
NS5.K314R_NS4A.A85T_map + scale_fill_manual(values = my.cols, na.value = "white")

#PLOT 6: distribution of NS4B-I240M ----
my.cols <- brewer.pal(6, "OrRd")
my.cols[1] <- '#ffffff'
NS4B.I240M_data <- aggregate(NS4B.I240M ~ county_FIPS, data = filtered_metadata, sum)
names(NS4B.I240M_data) <- c("region", "value")
NS4B.I240M_data$region <- as.numeric(as.character(NS4B.I240M_data$region))
NS4B.I240M_map <- choroplethr::county_choropleth(NS4B.I240M_data, title = "Distribution of NS4B.I240M", legend = "Presence of NS4B.I240M")
NS4B.I240M_map + scale_fill_manual(values = my.cols, na.value = "white")

#PLOT 7: distribution of NS5-A860T ----
my.cols <- brewer.pal(6, "OrRd")
my.cols[1] <- '#ffffff'
NS5.A860T_data <- aggregate(NS5.A860T ~ county_FIPS, data = filtered_metadata, sum)
names(NS5.A860T_data) <- c("region", "value")
NS5.A860T_data$region <- as.numeric(as.character(NS5.A860T_data$region))
NS5.A860T_map <- choroplethr::county_choropleth(NS5.A860T_data, title = "Distribution of NS5.A860T", legend = "Presence of NS5.A860T")
NS5.A860T_map + scale_fill_manual(values = my.cols, na.value = "white")

#PLOT 7: distribution of NS2A-R188K ----
NS2A.R188K_data <- aggregate(NS2A.R188K ~ county_FIPS, data = filtered_metadata, sum)
names(NS2A.R188K_data) <- c("region", "value")
NS2A.R188K_data$region <- as.numeric(as.character(NS2A.R188K_data$region))
NS2A.R188K_map <- choroplethr::county_choropleth(NS2A.R188K_data, title = "Distribution of NS2A.R188K", legend = "Presence of NS2A.R188K")
NS2A.R188K_map + scale_fill_manual(values = my.cols, na.value = "white")


##---------------------------- summary plot of mutations -------------------##
library(tidyverse)
library(dplyr)
polymorphisms <- read.csv("data_summary/AA_polymorphisms_dense.csv", header=TRUE)

# Function to replace the most common amino acid with NA, except for E.V159
replace_most_common <- function(x, name) {
  if (name == "E.V159") return(x)
  
  most_common <- names(which.max(table(x)))
  x[x == most_common] <- NA
  return(x)
}

# Apply the function to each column (except the first one)
polymorphisms_filtered <- polymorphisms
polymorphisms_filtered[-1] <- lapply(names(polymorphisms_filtered[-1]), function(name) replace_most_common(polymorphisms_filtered[[name]], name))

polymorphisms_long <- polymorphisms_filtered %>% pivot_longer(cols = -accession, names_to = "residue", values_to = "amino_acid")


# Remove rows with NA in the 'amino_acid' column
polymorphisms_long <- na.omit(polymorphisms_long)
# Remove rows with '-' in the 'amino_acid' column
polymorphisms_long <- polymorphisms_long %>% filter(amino_acid != '-')

# Convert 'residue' to a factor and specify the levels
polymorphisms_long$residue <- factor(polymorphisms_long$residue, levels = names(polymorphisms_filtered)[-1])

ggplot(polymorphisms_long, aes(x = residue, y = ..count../1211, fill = amino_acid)) +
  geom_bar(width = 0.7) +
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 20)]) +
  coord_cartesian(ylim = c(0, 0.25)) +  # Adjust these values as needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Residue", y = "Frequency")

##----------------------- time series for common mutations -------------#
library(lubridate)
library(gridExtra)

# Convert 'collection_date' to a date format and extract the year
filtered_metadata$year <- year(dmy(filtered_metadata$collection_date))

# Calculate the sum of each mutation by year
mutations_by_year <- filtered_metadata %>%
  group_by(year) %>%
  summarise(NS4A.A85T = sum(NS4A.A85T) / n(),
            NS2A.R188K = sum(NS2A.R188K) / n(),
            NS4B.I240M = sum(NS4B.I240M) / n(),
            NS5.K314R = sum(NS5.K314R) / n())

# Reshape the data to a long format
mutations_long <- mutations_by_year %>%
  gather(key = "mutation", value = "frequency", -year)

# Calculate the total number of samples per year
total_samples_per_year <- filtered_metadata %>%
  group_by(year) %>%
  summarise(total = n())

# Calculate the count of samples by region and year
region_by_year <- filtered_metadata %>%
  group_by(year, region_name) %>%
  summarise(count = n())

# Join the two dataframes and calculate the frequency
region_by_year <- region_by_year %>%
  left_join(total_samples_per_year, by = "year") %>%
  mutate(frequency = count / total)

# Create the first plot for NS5.K314R and NS4A.A85T
plot1 <- ggplot(mutations_long %>% filter(mutation %in% c("NS5.K314R", "NS4A.A85T")), 
                aes(x = year, y = frequency, color = mutation)) +
  geom_line() +
  labs(x = "Year", y = "Frequency of mutations", color = "Mutation") +
  theme_minimal()

# Create the second plot for NS2A.R188K and NS4B.I240M
plot2 <- ggplot(mutations_long %>% filter(mutation %in% c("NS2A.R188K", "NS4B.I240M")), 
                aes(x = year, y = frequency, color = mutation)) +
  geom_line() +
  labs(x = "Year", y = "Frequency of mutations", color = "Mutation") +
  theme_minimal()

# Create the third plot for the sampling of regions over time
plot3 <- ggplot(region_by_year, aes(x = year, y = frequency, color = region_name)) +
  geom_line() +
  labs(x = "Year", y = "Frequency of samples", color = "Region") +
  theme_minimal()

# Combine the plots
grid.arrange(plot1, plot2, plot3, ncol=1)