Clam resazurin metabolic rate analysis
================
AS Huffmyer
2026

This script analyzes resazurin assays for clam samples using a 48-well
plate format (BioTek Synergy HTX reader). Column 1 wells contain clam
samples; column 2 wells contain blanks.

# Set up

Set up workspace, set options, and load required packages.

``` r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

Load libraries.

``` r
library(tidyverse)
library(ggplot2)
library(cowplot)
library(tools)
```

`knitr::knit()` does not populate `params`; `rmarkdown::render()` does.
When knitting outside `render()`, use the same default as `params:` in
the YAML header above.

# Load data

## Parse plate reader files

Define a function to parse BioTek Synergy HTX .txt files. Each file may
contain one or more plate readings separated by a “Results” header. For
files with multiple readings the fluorescence values are averaged across
readings.

``` r
parse_plate_file <- function(file) {
  lines <- readLines(file, warn = FALSE)
  
  # Find all lines that contain only "Results"
  results_idx <- which(trimws(lines) == "Results")
  
  if (length(results_idx) == 0) {
    stop(paste("No 'Results' section found in file:", file))
  }
  
  readings_list <- list()
  
  for (idx in results_idx) {
    # Line idx   = "Results"
    # Line idx+1 = column header (tab-separated well columns; 6 cols for 24-well, 8 for 48-well)
    # Data rows: A–D (24-well) or A–F (48-well), then blank line or end of section
    j <- idx + 2
    data_lines <- character()
    while (j <= length(lines)) {
      line <- lines[j]
      if (nchar(trimws(line)) == 0) break
      if (!grepl("^[A-Fa-f]\\t", line)) break
      data_lines <- c(data_lines, line)
      j <- j + 1
    }

    plate_data <- lapply(data_lines, function(row) {
      parts <- unlist(strsplit(row, "\t"))
      row_letter <- trimws(parts[1])
      # Last field is filter label (e.g. 528/20,590/20); preceding fields are well values
      n <- length(parts)
      value_parts <- if (n >= 3) parts[seq.int(2, n - 1)] else parts[-1]
      values <- suppressWarnings(as.numeric(value_parts))
      data.frame(
        Row            = row_letter,
        sample_fluor   = values[1],  # column 1 = clam sample well
        blank_fluor    = values[2],  # column 2 = blank well
        stringsAsFactors = FALSE
      )
    })
    
    readings_list[[length(readings_list) + 1]] <- bind_rows(plate_data)
  }
  
  # If multiple readings exist, average col1 and col2 per row
  if (length(readings_list) > 1) {
    bind_rows(readings_list, .id = "reading") %>%
      group_by(Row) %>%
      summarise(
        sample_fluor = mean(sample_fluor, na.rm = TRUE),
        blank_fluor  = mean(blank_fluor,  na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    readings_list[[1]]
  }
}
```

Read and combine all plate files from `data/clam/<experiment>/` (or
`data/clam/` if `experiment` is empty). Filenames may look like
`plate1_t0.txt` or `plate01-T3.25.txt`.

``` r
find_repo_root <- function() {
  candidates <- c(
    getwd(),
    normalizePath(file.path(getwd(), "..")),
    normalizePath(file.path(getwd(), "..", ".."))
  )
  for (d in candidates) {
    if (dir.exists(file.path(d, "data", "clam"))) return(d)
  }
  getwd()
}

repo_root   <- find_repo_root()
exp_subdir  <- params$experiment
folder_path <- if (nzchar(exp_subdir)) {
  file.path(repo_root, "data", "clam", exp_subdir)
} else {
  file.path(repo_root, "data", "clam")
}

file_list <- list.files(path = folder_path, pattern = "\\.txt$",
                        full.names = TRUE, ignore.case = TRUE)

if (length(file_list) == 0) {
  stop("No .txt files found in ", folder_path, ". Check folder path and params$experiment.")
}

fig_dir <- file.path(repo_root, "figures", "clam", if (nzchar(exp_subdir)) exp_subdir else "default")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

data_list <- list()

for (file in file_list) {
  file_name <- file_path_sans_ext(basename(file))
  fn_lower  <- tolower(file_name)

  # plate01 / plate1 -> plate1 (for consistent plot colours)
  plate_num <- str_match(fn_lower, "plate0*(\\d+)")[, 2]
  plate     <- paste0("plate", as.integer(plate_num))

  # -T0, -T0.5, -T3.25 (hours)
  tp_str <- str_extract(fn_lower, "(?<=-t)[0-9.]+")
  timepoint <- as.numeric(tp_str)
  
  parsed <- parse_plate_file(file) %>%
    mutate(
      plate     = plate,
      timepoint = timepoint,
      well      = paste0(Row, "1")
    )
  
  data_list[[file_name]] <- parsed
}

combined_data <- bind_rows(data_list)

head(combined_data)
```

    ##   Row sample_fluor blank_fluor  plate timepoint well
    ## 1   A         3025        2239 plate1       0.5   A1
    ## 2   B         1315        2484 plate1       0.5   B1
    ## 3   C         3642        2485 plate1       0.5   C1
    ## 4   D         1695        1465 plate1       0.5   D1
    ## 5   A          592         537 plate1       0.0   A1
    ## 6   B          419         425 plate1       0.0   B1

## Normalize to T0

Divide each well’s fluorescence by the T0 (baseline) value to express
fluorescence as fold-change relative to the start of the assay.

``` r
combined_data <- combined_data %>%
  group_by(plate, well) %>%
  arrange(timepoint, .by_group = TRUE) %>%
  mutate(
    sample_norm = sample_fluor / first(sample_fluor),
    blank_norm  = blank_fluor  / first(blank_fluor)
  ) %>%
  ungroup()

head(combined_data)
```

    ## # A tibble: 6 × 8
    ##   Row   sample_fluor blank_fluor plate  timepoint well  sample_norm blank_norm
    ##   <chr>        <dbl>       <dbl> <chr>      <dbl> <chr>       <dbl>      <dbl>
    ## 1 A              592         537 plate1      0    A1           1          1   
    ## 2 A             3025        2239 plate1      0.5  A1           5.11       4.17
    ## 3 A             5323        2833 plate1      1    A1           8.99       5.28
    ## 4 A             7447        3729 plate1      2    A1          12.6        6.94
    ## 5 A             7805        4480 plate1      2.5  A1          13.2        8.34
    ## 6 A             8300        5168 plate1      3.25 A1          14.0        9.62

## Blank correction

Calculate the mean normalized blank fluorescence per plate and
timepoint, then subtract it from the normalized sample fluorescence
(adding 1 to keep values centred at 1 at T0).

``` r
blanks <- combined_data %>%
  group_by(plate, timepoint) %>%
  summarise(mean_blank_norm = mean(blank_norm, na.rm = TRUE), .groups = "drop")

combined_data <- left_join(combined_data, blanks, by = c("plate", "timepoint")) %>%
  mutate(fluorescence_corr = sample_norm - mean_blank_norm + 1)

head(combined_data)
```

    ## # A tibble: 6 × 10
    ##   Row   sample_fluor blank_fluor plate  timepoint well  sample_norm blank_norm
    ##   <chr>        <dbl>       <dbl> <chr>      <dbl> <chr>       <dbl>      <dbl>
    ## 1 A              592         537 plate1      0    A1           1          1   
    ## 2 A             3025        2239 plate1      0.5  A1           5.11       4.17
    ## 3 A             5323        2833 plate1      1    A1           8.99       5.28
    ## 4 A             7447        3729 plate1      2    A1          12.6        6.94
    ## 5 A             7805        4480 plate1      2.5  A1          13.2        8.34
    ## 6 A             8300        5168 plate1      3.25 A1          14.0        9.62
    ## # ℹ 2 more variables: mean_blank_norm <dbl>, fluorescence_corr <dbl>

# Plotting

## Raw fluorescence over time

``` r
plot_raw <- combined_data %>%
  ggplot(aes(x = timepoint, y = sample_fluor,
             colour = plate, group = interaction(plate, well))) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = c("plate1" = "steelblue", "plate2" = "darkorange")) +
  labs(
    x      = "Timepoint (h)",
    y      = "Raw Fluorescence",
    colour = "Plate"
  ) +
  theme_classic()

plot_raw
```

<figure>
<img src="figure/unnamed-chunk-8-1.png"
alt="plot of chunk unnamed-chunk-8" />
<figcaption aria-hidden="true">plot of chunk
unnamed-chunk-8</figcaption>
</figure>

``` r
ggsave(plot_raw, filename = file.path(fig_dir, "raw_fluorescence.png"), width = 6, height = 4)
```

## Normalized fluorescence over time

``` r
plot_norm <- combined_data %>%
  ggplot(aes(x = timepoint, y = sample_norm,
             colour = plate, group = interaction(plate, well))) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = c("plate1" = "steelblue", "plate2" = "darkorange")) +
  labs(
    x      = "Timepoint (h)",
    y      = "Fluorescence (fold change from T0)",
    colour = "Plate"
  ) +
  theme_classic()

plot_norm
```

<figure>
<img src="figure/unnamed-chunk-9-1.png"
alt="plot of chunk unnamed-chunk-9" />
<figcaption aria-hidden="true">plot of chunk
unnamed-chunk-9</figcaption>
</figure>

``` r
ggsave(plot_norm, filename = file.path(fig_dir, "normalized_fluorescence.png"), width = 6, height = 4)
```

## Blank-corrected fluorescence over time

``` r
plot_corr <- combined_data %>%
  ggplot(aes(x = timepoint, y = fluorescence_corr,
             colour = plate, group = interaction(plate, well))) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = c("plate1" = "steelblue", "plate2" = "darkorange")) +
  labs(
    x      = "Timepoint (h)",
    y      = "Blank-corrected Fluorescence",
    colour = "Plate"
  ) +
  theme_classic()

plot_corr
```

<figure>
<img src="figure/unnamed-chunk-10-1.png"
alt="plot of chunk unnamed-chunk-10" />
<figcaption aria-hidden="true">plot of chunk
unnamed-chunk-10</figcaption>
</figure>

``` r
ggsave(plot_corr, filename = file.path(fig_dir, "corrected_fluorescence.png"), width = 6, height = 4)
```

## Mean metabolic activity over time

Plot group means ± SE.

``` r
summary_data <- combined_data %>%
  group_by(plate, timepoint) %>%
  summarise(
    mean_fluor = mean(fluorescence_corr, na.rm = TRUE),
    se_fluor   = sd(fluorescence_corr,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

plot_mean <- summary_data %>%
  ggplot(aes(x = timepoint, y = mean_fluor,
             colour = plate, fill = plate)) +
  geom_ribbon(aes(ymin = mean_fluor - se_fluor, ymax = mean_fluor + se_fluor),
              alpha = 0.2, colour = NA) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = c("plate1" = "steelblue", "plate2" = "darkorange")) +
  scale_fill_manual(values   = c("plate1" = "steelblue", "plate2" = "darkorange")) +
  labs(
    x      = "Timepoint (h)",
    y      = "Blank-corrected Fluorescence (mean ± SE)",
    colour = "Plate",
    fill   = "Plate"
  ) +
  theme_classic()

plot_mean
```

<figure>
<img src="figure/unnamed-chunk-11-1.png"
alt="plot of chunk unnamed-chunk-11" />
<figcaption aria-hidden="true">plot of chunk
unnamed-chunk-11</figcaption>
</figure>

``` r
ggsave(plot_mean, filename = file.path(fig_dir, "mean_fluorescence.png"), width = 6, height = 4)
```

## Individual well trajectories faceted by plate

``` r
plot_facet <- combined_data %>%
  ggplot(aes(x = timepoint, y = fluorescence_corr,
             colour = well, group = well)) +
  facet_wrap(~plate) +
  geom_point() +
  geom_line() +
  labs(
    x      = "Timepoint (h)",
    y      = "Blank-corrected Fluorescence",
    colour = "Well"
  ) +
  theme_classic()

plot_facet
```

<figure>
<img src="figure/unnamed-chunk-12-1.png"
alt="plot of chunk unnamed-chunk-12" />
<figcaption aria-hidden="true">plot of chunk
unnamed-chunk-12</figcaption>
</figure>

``` r
ggsave(plot_facet, filename = file.path(fig_dir, "individual_wells.png"), width = 8, height = 4)
```
