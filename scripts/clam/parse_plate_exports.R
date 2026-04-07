suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(lubridate)
})

# Parse assay time in hours from filenames like plate-A-T2.5.txt.
parse_time_hr <- function(path) {
  file_name <- basename(path)
  hit <- str_match(file_name, "(?i)-T([0-9]+(?:\\.[0-9]+)?)\\.txt$")
  as.numeric(hit[, 2])
}

# Parse plate id from filenames like plate-A-T2.5.txt.
parse_plate_id <- function(path) {
  file_name <- basename(path)
  hit <- str_match(file_name, "(?i)^plate-([A-Za-z0-9-]+)-T[0-9]+(?:\\.[0-9]+)?\\.txt$")
  id <- hit[, 2]
  ifelse(is.na(id), "unknown", id)
}

parse_read_datetime <- function(lines) {
  date_line <- lines[str_detect(lines, "^Date\\t")][1]
  time_line <- lines[str_detect(lines, "^Time\\t")][1]

  if (is.na(date_line) || is.na(time_line)) {
    return(as.POSIXct(NA))
  }

  date_val <- trimws(str_split_fixed(date_line, "\\t", 2)[, 2])
  time_val <- trimws(str_split_fixed(time_line, "\\t", 2)[, 2])

  parsed <- suppressWarnings(mdy_hms(paste(date_val, time_val), tz = "UTC"))
  if (is.na(parsed)) {
    parsed <- suppressWarnings(mdy_hm(paste(date_val, time_val), tz = "UTC"))
  }
  parsed
}

parse_actual_temperature <- function(lines) {
  temp_line <- lines[str_detect(lines, "^Actual Temperature:")][1]
  if (is.na(temp_line)) {
    return(NA_real_)
  }
  as.numeric(str_extract(temp_line, "-?[0-9]+(?:\\.[0-9]+)?"))
}

parse_filter_label <- function(lines) {
  ex_line <- lines[str_detect(lines, "Excitation:")][1]
  em_line <- lines[str_detect(lines, "Emission:")][1]

  if (!is.na(ex_line) && !is.na(em_line)) {
    ex <- str_extract(ex_line, "[0-9]+/[0-9]+")
    em <- str_extract(em_line, "[0-9]+/[0-9]+")
    if (!is.na(ex) && !is.na(em)) {
      return(paste0(ex, ",", em))
    }
  }

  NA_character_
}

extract_results_block <- function(lines) {
  results_idx <- which(trimws(lines) == "Results")
  if (length(results_idx) == 0) {
    stop("No Results section found")
  }

  idx <- results_idx[1]
  if (idx + 1 > length(lines)) {
    stop("Results section missing header")
  }

  header_line <- lines[idx + 1]
  header_tokens <- strsplit(header_line, "\\t", fixed = FALSE)[[1]]
  header_tokens <- trimws(header_tokens)
  col_ids <- header_tokens[header_tokens != "" & str_detect(header_tokens, "^[0-9]+$")]

  if (length(col_ids) == 0) {
    stop("Could not identify numeric well columns in Results header")
  }

  j <- idx + 2
  data_lines <- character()
  while (j <= length(lines)) {
    line <- lines[j]
    trimmed <- trimws(line)
    if (trimmed == "") {
      break
    }
    if (!str_detect(line, "^[A-Za-z]\\t")) {
      break
    }
    data_lines <- c(data_lines, line)
    j <- j + 1
  }

  if (length(data_lines) == 0) {
    stop("No row data found under Results section")
  }

  list(col_ids = col_ids, data_lines = data_lines)
}

parse_results_rows <- function(col_ids, data_lines) {
  map_dfr(data_lines, function(line) {
    tokens <- strsplit(line, "\\t", fixed = FALSE)[[1]]
    tokens <- trimws(tokens)
    tokens <- tokens[tokens != ""]

    row_letter <- tokens[1]
    row_rest <- tokens[-1]

    nums <- suppressWarnings(as.numeric(row_rest))
    valid_num_idx <- which(!is.na(nums))
    if (length(valid_num_idx) == 0) {
      return(tibble())
    }

    vals <- nums[valid_num_idx]
    n <- min(length(vals), length(col_ids))
    vals <- vals[seq_len(n)]
    cols <- col_ids[seq_len(n)]

    tibble(
      row_id = toupper(row_letter),
      col_id = as.integer(cols),
      well_id = paste0(toupper(row_letter), cols),
      value = vals
    )
  })
}

parse_plate_export <- function(path) {
  lines <- read_lines(path, lazy = FALSE, progress = FALSE)

  results <- extract_results_block(lines)
  wells <- parse_results_rows(results$col_ids, results$data_lines)

  read_datetime <- parse_read_datetime(lines)
  actual_temperature <- parse_actual_temperature(lines)
  filter_label <- parse_filter_label(lines)

  wells %>%
    mutate(
      experiment_dir = basename(dirname(path)),
      source_file = path,
      source_name = basename(path),
      plate_id = str_to_lower(parse_plate_id(path)),
      time_hr = parse_time_hr(path),
      excitation_emission = if_else(is.na(filter_label), NA_character_, filter_label),
      read_datetime = read_datetime,
      actual_temperature = actual_temperature
    ) %>%
    select(any_of(c(
      "experiment_dir",
      "source_file",
      "source_name",
      "plate_id",
      "time_hr",
      "read_datetime",
      "actual_temperature",
      "excitation_emission",
      "row_id",
      "col_id",
      "well_id",
      "value"
    )))
}

parse_layout_txt <- function(layout_path) {
  # Preferred mode: tabular layout with explicit columns.
  tabular <- tryCatch(
    read_tsv(layout_path, show_col_types = FALSE, col_types = cols(.default = "c")),
    error = function(e) NULL
  )

  if (!is.null(tabular) && nrow(tabular) > 0) {
    names(tabular) <- names(tabular) %>%
      str_to_lower() %>%
      str_replace_all("[^a-z0-9]+", "_") %>%
      str_replace_all("_+", "_") %>%
      str_replace("_$", "")

    if (all(c("plate_well") %in% names(tabular))) {
      plate_col <- if ("plate_id" %in% names(tabular)) "plate_id" else NA_character_
      sample_col <- if ("sample_id" %in% names(tabular)) "sample_id" else NA_character_
      treatment_col <- if ("treatment_group" %in% names(tabular)) "treatment_group" else NA_character_
      shell_col <- if ("shell_color" %in% names(tabular)) "shell_color" else NA_character_
      size_col <- if ("size_mm" %in% names(tabular)) "size_mm" else NA_character_
      weight_col <- if ("weight_mg" %in% names(tabular)) "weight_mg" else NA_character_
      blank_col <- if ("is_blank" %in% names(tabular)) "is_blank" else NA_character_

      out <- tibble(
        plate_id = if (!is.na(plate_col)) tabular[[plate_col]] else NA_character_,
        well_id = toupper(tabular[["plate_well"]]),
        sample_label = if (!is.na(sample_col)) tabular[[sample_col]] else NA_character_,
        treatment = if (!is.na(treatment_col)) tabular[[treatment_col]] else NA_character_,
        shell_color = if (!is.na(shell_col)) tabular[[shell_col]] else NA_character_,
        treatment_group = if (!is.na(treatment_col)) tabular[[treatment_col]] else NA_character_,
        size_mm = suppressWarnings(as.numeric(if (!is.na(size_col)) tabular[[size_col]] else NA_character_)),
        weight_mg = suppressWarnings(as.numeric(if (!is.na(weight_col)) tabular[[weight_col]] else NA_character_)),
        is_blank = if (!is.na(blank_col)) {
          toupper(tabular[[blank_col]]) %in% c("TRUE", "T", "1", "YES", "Y")
        } else {
          NA
        },
        layout_status = "parsed_table",
        layout_raw = NA_character_
      ) %>%
        mutate(
          plate_id = str_remove(str_to_lower(plate_id), "^plate-"),
          plate_id = if_else(is.na(plate_id) | plate_id == "", NA_character_, plate_id),
          well_id = if_else(str_detect(well_id, "^[A-Z]+[0-9]+$"), well_id, NA_character_)
        )

      return(out)
    }
  }

  raw <- read_lines(layout_path, lazy = FALSE, progress = FALSE)
  lines <- raw[trimws(raw) != ""]

  if (length(lines) == 0) {
    return(tibble(
      plate_id = character(),
      well_id = character(),
      sample_label = character(),
      treatment = character(),
      shell_color = character(),
      treatment_group = character(),
      size_mm = numeric(),
      weight_mg = numeric(),
      is_blank = logical(),
      layout_status = character(),
      layout_raw = character()
    ))
  }

  # Mode 1: explicit well mapping, e.g. A1<TAB>control<TAB>treated
  split_candidates <- strsplit(lines, "[\\t,;]+")

  explicit <- map_dfr(split_candidates, function(parts) {
    p <- trimws(parts)
    if (length(p) >= 2 && str_detect(p[1], "^[A-Za-z]+[0-9]+$")) {
      tibble(
        plate_id = NA_character_,
        well_id = toupper(p[1]),
        sample_label = p[2],
        treatment = ifelse(length(p) >= 3, p[3], NA_character_),
        shell_color = NA_character_,
        treatment_group = ifelse(length(p) >= 3, p[3], NA_character_),
        size_mm = NA_real_,
        weight_mg = NA_real_,
        is_blank = NA,
        layout_status = "parsed_explicit",
        layout_raw = paste(parts, collapse = "\t")
      )
    } else {
      tibble()
    }
  })

  if (nrow(explicit) > 0) {
    return(explicit)
  }

  # Mode 2: matrix with row letters and numbered columns.
  first_parts <- trimws(split_candidates[[1]])
  numeric_header <- first_parts[first_parts != "" & str_detect(first_parts, "^[0-9]+$")]

  if (length(numeric_header) > 0) {
    matrix_rows <- map_dfr(split_candidates[-1], function(parts) {
      p <- trimws(parts)
      p <- p[p != ""]
      if (length(p) < 2 || !str_detect(p[1], "^[A-Za-z]+$")) {
        return(tibble())
      }

      row_letter <- toupper(p[1])
      cells <- p[-1]
      n <- min(length(cells), length(numeric_header))
      tibble(
        plate_id = NA_character_,
        well_id = paste0(row_letter, numeric_header[seq_len(n)]),
        sample_label = cells[seq_len(n)],
        treatment = NA_character_,
        shell_color = NA_character_,
        treatment_group = NA_character_,
        size_mm = NA_real_,
        weight_mg = NA_real_,
        is_blank = NA,
        layout_status = "parsed_matrix",
        layout_raw = paste(parts, collapse = "\t")
      )
    })

    if (nrow(matrix_rows) > 0) {
      return(matrix_rows)
    }
  }

  tibble(
    plate_id = NA_character_,
    well_id = NA_character_,
    sample_label = NA_character_,
    treatment = NA_character_,
    shell_color = NA_character_,
    treatment_group = NA_character_,
    size_mm = NA_real_,
    weight_mg = NA_real_,
    is_blank = NA,
    layout_status = "unparsed",
    layout_raw = paste(lines, collapse = " || ")
  )
}

parse_layout_for_experiment <- function(experiment_dir_path) {
  layout_path <- file.path(experiment_dir_path, "layout.txt")
  if (!file.exists(layout_path)) {
    return(tibble(
      plate_id = character(),
      well_id = character(),
      sample_label = character(),
      treatment = character(),
      shell_color = character(),
      treatment_group = character(),
      size_mm = numeric(),
      weight_mg = numeric(),
      is_blank = logical(),
      layout_status = character(),
      layout_raw = character(),
      layout_file = character()
    ))
  }

  parsed <- parse_layout_txt(layout_path)
  parsed %>% mutate(layout_file = layout_path)
}
