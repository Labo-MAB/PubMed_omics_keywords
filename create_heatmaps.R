#!/usr/bin/Rscript

# From the output of pubmed_trends_analysis, this script creates heatmaps of 
# keywords frequency over time and co-occurrence in pubmed articles. 

# Note:
# This script is an adaptation of the trends_plot.r which can be found at  
# https://github.com/lab42open-team/pubmed_trend_analysis/blob/master/scripts/trends_plots.r

# Original Author: Savvas Paragkamian (s.paragkamian@hcmr.gr)
#                  Institute of Marine Biology Biotechnology and Aquaculture (IMBBC)
#                  Hellenic Centre for Marine Research (HCMR)
# Original Date of creation: 2020-06-05


# Packages version
# tidyverse=1.3.0, readr=1.3.1, ggplot2=3.3.0, dplyr=0.8.5, Matrix=1.2-15

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(igraph)
  library(ggraph)
  library(pheatmap)
  library(tidygraph)
  library(svglite)
  library(viridis)
})

# Number of papers in pubmed
n_papers_pubmed <- 38096619

input_file <- "results/omics_2024-02-14_16-34_dig_analysis.tsv"
user_prefix <- "omics"
keywords_file"omics_keywords.txt"


trends_pubmed <- read_delim(input_file, delim = "\t", col_names = FALSE, 
                            col_types = cols())
colnames(trends_pubmed) <- c("year", "PMID", "synonym", "abstract_keywords")

trends_categories <- read_delim(keywords_file, delim = "\t", col_names = FALSE,
                                col_types = cols()) |>
                     arrange(X3)

colnames(trends_categories) <- c("synonym", "keyword", "category")

trends_categories_only <- trends_categories |> distinct(keyword,category)
## filter only the keywords that are listed in the trends_categories and then
## join them to keep the general categories. Also remove the the synonyms to
## keep only the unique number of PMIDs per keyword.

trends_pubmed <- trends_pubmed |>
  filter(synonym %in% trends_categories$synonym) |>
  dplyr::left_join(trends_categories, by = c("synonym" = "synonym")) |>
  dplyr::distinct(year, PMID, keyword, category)

# ------------------------------------------------------------------------------
## trends per year

keywords_per_year <- trends_pubmed |>
  distinct(PMID, keyword, category, year) |>
  group_by(year, keyword, category) |>
  summarize(counts = n(), .groups = "keep") |>
  ungroup() |>
  arrange(year) |>
  group_by(keyword, category) |>
  mutate(cumulative_count = cumsum(counts)) |>
  ungroup() |>
  mutate(count_bin = cut(counts, breaks = c(0, 50, 250, 500, 1000, 2500, 5000,
                                            max(counts, na.rm = TRUE)),
                         labels = c("1-50", "50-250", "250-500", "500-1000",
                                    "1000-2500", "2500-5000", "5000+")))
# original breaks: 0, 10, 50, 100, 500, 2000, max(counts, na.rm = TRUE)

# change the order to descreasing to appear with alphabetical order
keywords_per_year$keyword <- factor(keywords_per_year$keyword,
                                    levels = unique(keywords_per_year$keyword[order(keywords_per_year$category,
                                                                                  keywords_per_year$keyword,
                                                                                  decreasing = TRUE)]))

# ------------------------------------------------------------------------------
# timeline heatmap

# Filter to keep only match from 2010 to last year
filtered_yearly_keywords <- keywords_per_year |>
  filter(2009 < year & year < 2024)

pubmed_keyword_per_year_heatmap <- ggplot() +
  geom_tile(data = filtered_yearly_keywords,
            aes(x = year, y = keyword,
                fill = count_bin,
                height = 1,
                width = 1),
            color = "white",
            show.legend = TRUE) +
  ggtitle("Occurence of keywords over years") +
  scale_x_continuous(breaks = seq(2010, 2022, 2), n.breaks = 10) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  scale_fill_manual(values = c("#dadaeb", "#bcbddc", "#9e9ac8", "#807dba",
                               "#6a51a3", "#4a1f7c", "#2b0955")) +
  ylab("") +
  xlab("") +
  theme_bw() +
  guides(fill = guide_legend(title = "Number of abstracts")) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 13),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.8, "cm"),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.title = element_text(hjust = 0.5, size = 19),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())

ggsave(paste0("results/",
              user_prefix,
              "words_over_time.svg"),
       plot = pubmed_keyword_per_year_heatmap,
       width = 45,
       height = 15,
       units = "cm",
       device = "svg",
       dpi = 300)

# ------------------------------------------------------------------------------
# Data preparation for plots

# create the edglist of keywords and PMID's
papers_keywords_network <- trends_pubmed |>
  group_by(PMID, keyword) |>
  distinct(PMID, keyword) |>
  ungroup() |>
  left_join(trends_categories_only, by = c("keyword" = "keyword"))

keywords_n_papers <- papers_keywords_network |>
  group_by(keyword) |>
  summarise(n_papers = n()) |>
  mutate(freq = n_papers / n_papers_pubmed)

# very important for the correct order of the keywords based on categories.
papers_keywords_network$keyword <- factor(papers_keywords_network$keyword,
                                          levels = unique(papers_keywords_network$keyword[order(papers_keywords_network$category,
                                                                                              papers_keywords_network$keyword)]))

# create a matrix class spMatrix (handles better sparse matrices) to do inverse table multiplication
papers_keywords_matrix <- spMatrix(nrow = length(unique(papers_keywords_network$PMID)),
                                   ncol = length(unique(papers_keywords_network$keyword)),
                                   i = as.numeric(factor(papers_keywords_network$PMID)),
                                   j = as.numeric(factor(papers_keywords_network$keyword)),
                                   x = rep(1, length(as.numeric(papers_keywords_network$PMID))))

row.names(papers_keywords_matrix) <- levels(factor(papers_keywords_network$PMID))
colnames(papers_keywords_matrix) <- levels(factor(papers_keywords_network$keyword))

# with the inverse cross product we do the projection of the edgelist to
# keywords in order to calculate how many times keyword pairs appear together in
# abstracts.
keywords_heatmap <- tcrossprod(t(papers_keywords_matrix))

# becaue the matrix is summetric we keep the triangle
keywords_heatmap[upper.tri(keywords_heatmap)] <- 0

keywords_heatmap <- as.data.frame(as.matrix(keywords_heatmap))
#write_delim(keywords_heatmap,"keywords_heatmap.tsv",delim="\t")

# transform to long format for plotting and remove zero's and NA's and assign -1
# to loops (self occurrence)
keywords_heatmap_long <- as.data.frame(as.matrix(keywords_heatmap)) |>
    rownames_to_column() |>
    pivot_longer(-rowname, names_to = "colname", values_to = "count") |>
    # filter(count != 0, colname != rowname) |>
    filter(count != 0) |>
    na.omit()

colnames(keywords_heatmap_long) <- c("from", "to", "count")

# original breaks 0, 10, 50, 100, 500, 800, 1000
keywords_heatmap_long$diag <- keywords_heatmap_long$from == keywords_heatmap_long$to
keywords_heatmap_long$count_ref <- keywords_heatmap_long$count * ifelse(keywords_heatmap_long$diag, -1, 1)
keywords_heatmap_long$count_bin <- cut(keywords_heatmap_long$count_ref,
                                       breaks = c(min(keywords_heatmap_long$count_ref, na.rm = TRUE)-1,
                                                  -25000, -10000, -5000, -1000, 
                                                  0, 100, 250, 500, 1000, 2500,
                                                  max(keywords_heatmap_long$count_ref, na.rm = TRUE)),
                                       labels = c("25000+", "10000-25000", "5000-1000", "1-1000", "0",
                                                  "1-100", "100-250", "250-500",
                                                  "500-1000", "1000-2500", "2500+"))
# assign the order levels of the count_bin
keywords_heatmap_long$count_bin <- factor(as.character(keywords_heatmap_long$count_bin),
                                          levels=rev(levels(keywords_heatmap_long$count_bin)))

## add the categories for the keywords
keywords_heatmap_long <- keywords_heatmap_long %>%
    left_join(trends_categories_only, by = c("from" = "keyword")) %>%
    left_join(trends_categories_only, by = c("to" = "keyword"))

# elements of the plot
keywords <- trends_categories |>
  filter(keyword %in% unique(c(keywords_heatmap_long$from,
                               keywords_heatmap_long$to))) |>
  mutate(from = factor(keyword,
                       levels = as.character(unique(keyword)))) |>
  mutate(to = factor(keyword,
                     levels = as.character(unique(unique(keyword)))),
         count = 0,
         count_bin = "0",
         jaccard = 0) |>
  dplyr::select(from, to, count, count_bin)

diagonal <- tibble(from = factor(keywords$from,
                                 levels = unique(keywords_heatmap_long$from[order(keywords_heatmap_long$category.x,keywords_heatmap_long$from)])),
                   to = factor(keywords$to,
                               levels = unique(keywords_heatmap_long$to[order(keywords_heatmap_long$category.y,
                                                                                  keywords_heatmap_long$to)])),
                   count = -1,
                   jaccard = 0)

# we defined here the diagonal because the raw values don't include them.
# In addition we need the diagonal seperate from the raw data because we will paint it differently

## summaries to dynamically set the break points and limits of the plot

summary <- summary(keywords_heatmap_long$count)

# For the heatmap we need breaks to define the specific points of the legend
# and limits to ensure that all values will be included in the plot.
# To create breaks and limits and make them scalable we used the base R
# functions summary and quantile.
# also quantiles because the raw counts are far apart

# big probabilities because of order of magnitude difference of values
quantile <- as.vector(quantile(keywords_heatmap_long$count,
                               probs = c(50, 90, 95, 98) / 100))

breaks <- c(floor(min(summary)),
            round(quantile[1]),
            round(quantile[2]),
            round(quantile[3]),
            round(quantile[4]),
            ceiling(max(summary)))

limits <- c(min(breaks), max(breaks))

# add the order based on the categories so they appear in that order
keywords_heatmap_long$from <- factor(keywords_heatmap_long$from,
                                     levels = unique(keywords_heatmap_long$from[order(keywords_heatmap_long$category.x,
                                                                                      keywords_heatmap_long$from)]))
keywords_heatmap_long$to <- factor(keywords_heatmap_long$to, 
                                   levels = unique(keywords_heatmap_long$to[order(keywords_heatmap_long$category.y,
                                                                                  keywords_heatmap_long$to)]))


# ------------------------------------------------------------------------------
# Plotting Concurrence heatmap

# legend title and combination of different colors and shapes into one legend
# g <- guide_legend("no of abstracts")

keywords_heatmap_long <- keywords_heatmap_long %>% filter(!is.na(count_bin))

pubmed_keyword_coocurrence_heatmap <- ggplot() +
  geom_tile(data = keywords_heatmap_long,
            aes(x = from, y = to, fill = count_bin, width = .98, height = .98),
            alpha = 0.75,
            show.legend = TRUE) +
  scale_fill_manual(breaks = c("0", "1-1000", "5000-1000", "10000-25000", "25000+", 
                                                  "1-100", "100-250", "250-500",
                                                  "500-1000", "1000-2500", "2500+"),
                    values = c("#fcfcfc", "#c7c7c7", "#949494", "#656565", "#393939",   
                             "#fafa6e", "#9cdf7c", "#4abd8c",
                               "#00968e", "#106e7c", "#395b6e"),
                    drop = FALSE) +
  coord_fixed() +
  scale_x_discrete(position = "top",
                   limits = rev(levels(keywords_heatmap_long$from))) +
  scale_y_discrete(limits = rev(levels(keywords_heatmap_long$to))) +
  guides(fill = guide_legend("Number of abstracts", ncol = 2, reverse = TRUE), 
         color = FALSE) +
  xlab("") +
  ylab("") +
  theme_bw() +
  ggtitle("Co-occurrence of keywords in abstract") +
  theme(plot.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15), 
        axis.text.x = element_text(angle = 45, hjust = 0),
        legend.position = c(1.1, .2),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 19),
        plot.margin = unit(c(0.2, 5, 1, 1), "cm"))

  ggsave(paste0("results/",
              user_prefix,
              "_coocurences.svg"),
       plot = pubmed_keyword_coocurrence_heatmap,
       device = "svg",
       dpi = 150)


#write_delim(keywords_heatmap_long,"heatmap_data.txt", delim="\t")


# ------------------------------------------------------------------------------
# running the log2 heatmap
### log2 transformation was the best way to gather together counts. Log10 was
# too aggresive and sqrt too soft.

# Remove the diagonal and calculate the log2 of the counts
keywords_heatmap_long$log2 <- log2(keywords_heatmap_long$count)

## summaries to dynamically set the break points and limits of the plot
summary_log <- summary(keywords_heatmap_long$log2)
breaks_log <- unname(c(floor(summary_log[1]),
                       round(summary_log[2]),
                       round(summary_log[4]),
                       round(summary_log[5]),
                       ceiling(max(summary_log))))

# Try hardcoded breaks values
breaks_log <- c(5, 8, 11, 14, 17)
limits_log <- c(min(breaks_log), max(breaks_log))
# legend title and combination of different colors and shapes into one legend
g_log <- guide_legend("log2(no of abstracts)")

pubmed_keyword_coocurrence_heatmap <- ggplot() +
  geom_tile(data = keywords_heatmap_long,
            aes(x = from,
                y = to,
                fill = log2),
            alpha = 0,
            show.legend = FALSE) +
  geom_point(data = keywords_heatmap_long |> filter(diag == FALSE),
             aes(x = from,
                 y = to,
                 colour = log2,
                 size = log2))  +
  scale_size(name = "co-occurrence",
             range = c(7, 25),
             breaks = breaks_log,
             limits = limits_log) +
  scale_colour_gradientn(colours = c("steelblue1", "yellowgreen", "yellow",
                                    "goldenrod1", "orange"),
                         breaks = breaks_log,
                         limits = limits_log) +
  geom_point(data = keywords_heatmap_long |> filter(diag == TRUE),
             aes(x = from,
                 y = to,
                 size = log2),
             color = "grey50",
             show.legend = FALSE) +
  scale_x_discrete(position = "top",
                   limits = rev(levels(keywords_heatmap_long$from))) +
  scale_y_discrete(limits = rev(levels(keywords_heatmap_long$to))) +
  guides(colour = g_log, size = g_log) +
  ggtitle("Co-occurrence of keywords in abstract (log2)") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 19),
        legend.position = c( .83, .35),
        plot.margin = unit(c(1, 2, 1, 2), "cm"))

ggsave(paste0("results/",
              user_prefix,
              "_log_cooccurences.svg"),
       plot = pubmed_keyword_coocurrence_heatmap,
       device = "svg",
       dpi = 150)


# ------------------------------------------------------------------------------
# Jaccard Index

# Jaccard index is the intersection over the union. So we join for each node -
# keyword the total occurrencies. The join is double because we have two columns
# of keywords and this way is easier for the calculations

keywords_heatmap_jaccard <- keywords_heatmap_long |>
    filter(diag == FALSE) |>
    left_join(keywords_n_papers, c("from" = "keyword")) |>
    left_join(keywords_n_papers, c("to" = "keyword")) |>
    mutate(jaccard_index = count / (n_papers.x + n_papers.y - count),
           random_expectation = (count / n_papers_pubmed) / (freq.x * freq.y))

colors_j <- c("steelblue1","yellow","goldenrod1","orange")
colors_j <- rev(viridis(4, alpha = 1, begin = 0.18, end = 1, direction = 1, option = "viridis"))
limits <- c(0,
           round(max(keywords_heatmap_jaccard$jaccard_index[keywords_heatmap_jaccard$jaccard_index < 1]) * 1.07,4))
# the limits are the maximum value multiplied by >1 so it is bigger than 
# the maximum value

pubmed_jaccard_heatmap <- ggplot() + 
  geom_tile(data = keywords_heatmap_jaccard,
            aes(x = from,
                y = to,
                fill = jaccard_index,
                width = .98, height = .98),
            alpha = 1,
            show.legend = TRUE,
            colour = "white") +
  scale_fill_gradientn(colours = colors_j, limits = limits) +
  scale_x_discrete(position = "top",
                  limits = rev(unique(keywords_heatmap_jaccard$from)),
                  expand = expansion(add = c(1.3, 1.3))) +
  scale_y_discrete(limits = rev(unique(keywords_heatmap_jaccard$to)),
                              expand = expansion(add = c(1.3, 1.3))) +
  ggtitle("Jaccard similarity between keywords") +
  xlab("") +
  ylab("") +
  guides(fill = guide_legend(title = "Jaccard similarity")) +
  theme_bw() +
    theme(plot.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(size = 15), 
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5),
      legend.position = c(.83, .20),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 19),
      plot.margin = unit(c(0.2, 5, 1, 1), "cm"))

ggsave(paste0("results/",
              user_prefix,
              "_jaccard.svg"),
       plot = pubmed_jaccard_heatmap,
       device = "svg",
       dpi = 150)
