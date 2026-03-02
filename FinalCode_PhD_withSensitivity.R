# ============================================================
# Title: MLS–TLS–Field tree matching, comparison metrics,
#        plots, statistical tests, regression, and MC uncertainty
# Author: Larissa Maria Granja
# ============================================================

# ============================================================
# Chunk 0 — Packages + global settings
# Description:
# Loads required packages, installs missing ones, defines file paths,
# and sets stand groupings used throughout the workflow.
# ============================================================

pkgs <- c(
  "dplyr", "readr", "stringr", "purrr", "tidyr", "ggplot2", "RANN",
  "broom", "patchwork", "clue" 
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install)

invisible(lapply(pkgs, library, character.only = TRUE))

# Optional package for Lin's Concordance Correlation Coefficient (CCC)
if (!("DescTools" %in% installed.packages()[, "Package"])) {
  install.packages("DescTools")
}
library(DescTools)

# ---------- User paths ----------
paths <- list(
  field = "C:/Users/larac/Desktop/3rd Year Report_PhD/AllData_ForR/field_clean.csv",
  mls   = "C:/Users/larac/Desktop/3rd Year Report_PhD/AllData_ForR/mls_3dfin_all.csv",
  tls   = "C:/Users/larac/Desktop/3rd Year Report_PhD/AllData_ForR/tls_3dfin_all.csv",
  stripes_dir = "C:/Users/larac/Desktop/3rd Year Report_PhD/AllData_ForR/TLS_StemStripes_Tables",
  out_dir = "C:/Users/larac/Desktop/3rd Year Report_PhD/R_outputs_27Feb_V3"
)

dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)

# Stands present in 3DFin outputs
stands_keep <- c("Oregon", "Lawson1", "Lawson2", "Alerce", "Roble")

# Stands where Field comparisons are expected
stands_field_ok <- c("Oregon", "Lawson1", "Lawson2", "Alerce", "Roble")

# Stands where a volume equation is implemented
stands_with_volume <- c("Oregon", "Alerce", "Roble")


# ============================================================
# Chunk 1 — Robust reader + helper functions
# Description:
# Defines generic file reading functions (delimiter detection),
# stand-name cleaning, plausibility filters, and basic error metrics.
# ============================================================

guess_delim_from_line <- function(line) {
  if (grepl(";", line)) return(";")
  if (grepl("\t", line)) return("\t")
  if (grepl(",", line)) return(",")
  return("whitespace")
}

read_any_delim <- function(file, comment_prefix = NULL) {
  raw <- readLines(file, warn = FALSE)
  raw <- raw[nzchar(raw)]
  
  if (!is.null(comment_prefix)) {
    raw <- raw[!grepl(paste0("^\\s*", comment_prefix), raw)]
  }
  
  stopifnot(length(raw) > 1)
  delim <- guess_delim_from_line(raw[1])
  
  if (delim == "whitespace") {
    readr::read_table(
      file = I(paste(raw, collapse = "\n")),
      col_types = cols(.default = col_guess()),
      show_col_types = FALSE
    )
  } else {
    readr::read_delim(
      file = I(paste(raw, collapse = "\n")),
      delim = delim,
      col_types = cols(.default = col_guess()),
      show_col_types = FALSE,
      trim_ws = TRUE
    )
  }
}

# Standardize stand labels: e.g., "Lawson 1" -> "Lawson1"
clean_stand <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("\\s+", "")
}

# Plausibility filter for 3DFin outputs:
# DBH in meters: (0, 3]; Height in meters: (1, 80]
is_plausible_3dfin <- function(dbh_m, h_m) {
  !is.na(dbh_m) & !is.na(h_m) &
    dbh_m > 0 & dbh_m <= 3 &
    h_m > 1 & h_m <= 80
}

# Error summary functions
rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))
bias <- function(x) mean(x, na.rm = TRUE)
mae  <- function(x) mean(abs(x), na.rm = TRUE)

# ------------------------------------------------------------
# One-to-one XY matching (Hungarian algorithm)
#Description: 
#Returns globally optimal 1:1 matches under a max distance threshold
## ------------------------------------------------------------
match_xy_one_to_one <- function(query_df, ref_df,
                                qx, qy, rx, ry,
                                q_id, r_id,
                                max_dist_m,
                                stand_col = "stand_id") {
  
  if (nrow(query_df) == 0 || nrow(ref_df) == 0) return(tibble())
  
  # Keep only rows with finite coordinates
  q_ok <- is.finite(query_df[[qx]]) & is.finite(query_df[[qy]])
  r_ok <- is.finite(ref_df[[rx]])   & is.finite(ref_df[[ry]])
  
  query_df <- query_df[q_ok, , drop = FALSE]
  ref_df   <- ref_df[r_ok, , drop = FALSE]
  
  if (nrow(query_df) == 0 || nrow(ref_df) == 0) return(tibble())
  
  # Distance matrix (query rows x ref rows)
  dx <- outer(query_df[[qx]], ref_df[[rx]], "-")
  dy <- outer(query_df[[qy]], ref_df[[ry]], "-")
  D  <- sqrt(dx^2 + dy^2)
  
  BIG <- 1e9
  D[!is.finite(D)] <- BIG
  
  nr <- nrow(D)
  nc <- ncol(D)
  n  <- max(nr, nc)
  
  # Pad to square matrix for Hungarian solver
  cost <- matrix(BIG, nrow = n, ncol = n)
  cost[1:nr, 1:nc] <- D
  
  # Hungarian assignment
  assign <- clue::solve_LSAP(cost)
  j <- as.integer(assign)[1:nr]  # assigned ref column for each query row
  
  # Safe extraction (avoid index bug)
  dist_vec <- rep(NA_real_, nr)
  ref_ids  <- rep(NA_character_, nr)
  
  valid <- which(j >= 1 & j <= nc)
  if (length(valid) > 0) {
    dist_vec[valid] <- D[cbind(valid, j[valid])]
    ref_ids[valid]  <- as.character(ref_df[[r_id]][j[valid]])
  }
  
  out <- tibble(
    !!stand_col := query_df[[stand_col]],
    query_id = as.character(query_df[[q_id]]),
    ref_id   = ref_ids,
    dist_m   = dist_vec
  ) %>%
    filter(is.finite(dist_m), dist_m <= max_dist_m)
  
  out
}

# ============================================================
# Chunk 2 — Read + standardize Field / MLS / TLS tables
# Description:
# Imports the datasets and harmonizes variable names/types.
# Converts 3DFin DBH from meters to centimeters.
# ============================================================

# ---- Field (DBH already in cm) ----
field <- read_any_delim(paths$field) %>%
  rename_with(tolower) %>%
  transmute(
    tree_id_field = as.character(tree_id_field),
    species = as.character(species),
    dbh_field_cm = as.numeric(dbh_field),
    h_field_m = as.numeric(h_field),
    x_field = as.numeric(x_field),
    y_field = as.numeric(y_field)
  )

# ---- MLS raw ----
mls_raw <- read_any_delim(paths$mls) %>% rename_with(tolower)

mls <- mls_raw %>%
  mutate(
    stand_id = clean_stand(stand_id),
    dbh_m = as.numeric(dbh_m),
    h = as.numeric(h)
  ) %>%
  filter(stand_id %in% stands_keep) %>%
  filter(is_plausible_3dfin(dbh_m, h)) %>%
  transmute(
    stand_id,
    stem_id = as.character(stem_id),
    species = as.character(species),
    x_mls = as.numeric(x_mls),
    y_mls = as.numeric(y_mls),
    dbh_mls_cm = dbh_m * 100,
    h_mls_m = h
  )

# ---- TLS raw ----
tls_raw <- read_any_delim(paths$tls) %>% rename_with(tolower)

tls <- tls_raw %>%
  mutate(
    stand_id = clean_stand(stand_id),
    dbh_h = as.numeric(dbh_h),
    h_m = as.numeric(h_m)
  ) %>%
  filter(stand_id %in% stands_keep) %>%
  filter(is_plausible_3dfin(dbh_h, h_m)) %>%
  transmute(
    stand_id,
    stem_id = as.character(stem_id),
    species = as.character(species),
    x_tls = as.numeric(x_tls),
    y_tls = as.numeric(y_tls),
    dbh_tls_cm = dbh_h * 100,
    h_tls_m = h_m
  )


# ============================================================
# Chunk 2.5 — Stand-specific QC filters
# Description:
# Applies stand-specific filtering thresholds for DBH and height
# to remove outliers and implausible values before matching.
# ============================================================

mls <- mls %>%
  filter(!(stand_id == "Alerce"  & (dbh_mls_cm < 5   | dbh_mls_cm > 34 | h_mls_m < 4    | h_mls_m > 16.6))) %>%
  filter(!(stand_id == "Lawson1" & (dbh_mls_cm < 13.5| dbh_mls_cm > 60 | h_mls_m < 10.8 | h_mls_m > 35))) %>%
  filter(!(stand_id == "Lawson2" & (dbh_mls_cm < 13.5| dbh_mls_cm > 60 | h_mls_m < 10.8 | h_mls_m > 35))) %>%
  filter(!(stand_id == "Oregon"  & (dbh_mls_cm < 5.3 | dbh_mls_cm > 70 | h_mls_m < 5.3  | h_mls_m > 28.5))) %>%
  filter(!(stand_id == "Roble"   & (dbh_mls_cm < 12  | dbh_mls_cm > 66 | h_mls_m < 13.3 | h_mls_m > 39)))

tls <- tls %>%
  filter(!(stand_id == "Alerce"  & (dbh_tls_cm < 5   | dbh_tls_cm > 34 | h_tls_m < 4    | h_tls_m > 16.6))) %>%
  filter(!(stand_id == "Lawson1" & (dbh_tls_cm < 13.5| dbh_tls_cm > 60 | h_tls_m < 10.8 | h_tls_m > 35))) %>%
  filter(!(stand_id == "Lawson2" & (dbh_tls_cm < 13.5| dbh_tls_cm > 60 | h_tls_m < 10.8 | h_tls_m > 35))) %>%
  filter(!(stand_id == "Oregon"  & (dbh_tls_cm < 5.3 | dbh_tls_cm > 70 | h_tls_m < 5.3  | h_tls_m > 28.5))) %>%
  filter(!(stand_id == "Roble"   & (dbh_tls_cm < 12  | dbh_tls_cm > 66 | h_tls_m < 13.3 | h_tls_m > 39)))


# ============================================================
# Chunk 3 — MLS → TLS matching (strict 1:1 using Hungarian assignment)
# ============================================================

# Optional: check duplicate IDs in source tables
check_source_duplicates <- function() {
  cat("\n--- Source duplicate ID checks ---\n")
  
  dup_mls <- mls %>% count(stand_id, stem_id) %>% filter(n > 1)
  dup_tls <- tls %>% count(stand_id, stem_id) %>% filter(n > 1)
  dup_fld <- field %>% count(tree_id_field) %>% filter(n > 1)
  
  cat("MLS duplicate stem_id within stand:", nrow(dup_mls), "rows\n")
  cat("TLS duplicate stem_id within stand:", nrow(dup_tls), "rows\n")
  cat("Field duplicate tree_id_field:", nrow(dup_fld), "rows\n")
  
  if (nrow(dup_mls) > 0) print(dup_mls)
  if (nrow(dup_tls) > 0) print(dup_tls)
  if (nrow(dup_fld) > 0) print(dup_fld)
}

# Run once after field/mls/tls are loaded
check_source_duplicates()

map_mls_to_tls_one <- function(stand, max_dist_m = 2.0) {
  m <- mls %>%
    filter(stand_id == stand) %>%
    filter(is.finite(x_mls), is.finite(y_mls))
  
  t <- tls %>%
    filter(stand_id == stand) %>%
    filter(is.finite(x_tls), is.finite(y_tls))
  
  out <- match_xy_one_to_one(
    query_df = m,
    ref_df   = t,
    qx = "x_mls", qy = "y_mls",
    rx = "x_tls", ry = "y_tls",
    q_id = "stem_id", r_id = "stem_id",
    max_dist_m = max_dist_m
  )
  
  out %>%
    transmute(
      stand_id,
      stem_id_mls = query_id,
      stem_id_tls = ref_id,
      dist_m
    )
}

map_mls_tls <- purrr::map_dfr(stands_keep, ~map_mls_to_tls_one(.x, max_dist_m = 2.0))

# QC for MLS↔TLS 1:1
qc_map_mls_tls <- map_mls_tls %>%
  summarise(
    n = n(),
    u_mls = n_distinct(paste(stand_id, stem_id_mls, sep = "__")),
    u_tls = n_distinct(paste(stand_id, stem_id_tls, sep = "__"))
  )
print(qc_map_mls_tls)

if (qc_map_mls_tls$n != qc_map_mls_tls$u_mls || qc_map_mls_tls$n != qc_map_mls_tls$u_tls) {
  warning("MLS<->TLS mapping is not strictly 1:1 within stand. Check source duplicates or matching thresholds.")
}

readr::write_csv(map_mls_tls, file.path(paths$out_dir, "map_mls_to_tls.csv"))


# ============================================================
# Chunk 4 — TLS → Field matching (strict 1:1 using Hungarian assignment)
# ============================================================

# Species/stand filter for Field data
field_for_stand <- function(stand, field_df) {
  sp <- stringr::str_to_lower(field_df$species)
  
  if (stand == "Alerce")  return(field_df[grepl("alerce|fitzroya", sp), ])
  if (stand == "Oregon")  return(field_df[grepl("oregon|pseudotsuga|douglas", sp), ])
  if (stand == "Roble")   return(field_df[grepl("roble|nothofagus", sp), ])
  if (stand %in% c("Lawson1", "Lawson2")) {
    return(field_df[grepl("lawson|chamaecyparis|cipres|ciprés|falso", sp), ])
  }
  
  field_df
}

map_tls_to_field_one <- function(stand, max_dist_m = 5.0) {
  t <- tls %>%
    filter(stand_id == stand) %>%
    filter(is.finite(x_tls), is.finite(y_tls))
  
  f <- field %>%
    filter(is.finite(x_field), is.finite(y_field)) %>%
    field_for_stand(stand, .) %>%
    mutate(stand_id = stand)  # for consistent output
  
  out <- match_xy_one_to_one(
    query_df = t,
    ref_df   = f,
    qx = "x_tls", qy = "y_tls",
    rx = "x_field", ry = "y_field",
    q_id = "stem_id", r_id = "tree_id_field",
    max_dist_m = max_dist_m
  )
  
  out %>%
    transmute(
      stand_id,
      stem_id_tls   = query_id,
      tree_id_field = ref_id,
      dist_m_field  = dist_m
    )
}

map_tls_field <- purrr::map_dfr(stands_field_ok, ~map_tls_to_field_one(.x, max_dist_m = 5.0))

# QC for TLS↔Field 1:1
qc_map_tls_field <- map_tls_field %>%
  summarise(
    n = n(),
    u_tls   = n_distinct(paste(stand_id, stem_id_tls, sep = "__")),
    u_field = n_distinct(paste(stand_id, tree_id_field, sep = "__"))
  )
print(qc_map_tls_field)

if (qc_map_tls_field$n != qc_map_tls_field$u_tls || qc_map_tls_field$n != qc_map_tls_field$u_field) {
  warning("TLS<->Field mapping is not strictly 1:1 within stand. Check source duplicates or matching thresholds.")
}

readr::write_csv(map_tls_field, file.path(paths$out_dir, "map_tls_to_field.csv"))

# ============================================================
# Chunk 5 — Build matched MLS + TLS + Field table (triplets)
# Description:
# Combines the MLS↔TLS and TLS↔Field pairings to produce a final
# three-way matched table, then joins DBH, height, and coordinates.
# ============================================================

if (nrow(map_mls_tls) == 0) warning("map_mls_tls is empty")
if (nrow(map_tls_field) == 0) warning("map_tls_field is empty")

stands_with_field <- unique(map_tls_field$stand_id)

match_all <- map_mls_tls %>%
  filter(stand_id %in% stands_with_field) %>%
  inner_join(map_tls_field, by = c("stand_id", "stem_id_tls"))

if (nrow(match_all) == 0) {
  warning("No overlapping MLS-TLS-Field triplets found. Field-based analyses/plots will be skipped.")
}

matched <- match_all %>%
  left_join(
    mls %>% select(stand_id, stem_id_mls = stem_id, species_mls = species, dbh_mls_cm, h_mls_m, x_mls, y_mls),
    by = c("stand_id", "stem_id_mls")
  ) %>%
  left_join(
    tls %>% select(stand_id, stem_id_tls = stem_id, species_tls = species, dbh_tls_cm, h_tls_m, x_tls, y_tls),
    by = c("stand_id", "stem_id_tls")
  ) %>%
  left_join(
    field %>% select(tree_id_field, species_field = species, dbh_field_cm, h_field_m, x_field, y_field),
    by = "tree_id_field"
  )

# Final QC: strict 1:1:1?
dup_check <- matched %>%
  summarise(
    n_rows = n(),
    n_unique_mls   = n_distinct(paste(stand_id, stem_id_mls, sep = "__")),
    n_unique_tls   = n_distinct(paste(stand_id, stem_id_tls, sep = "__")),
    n_unique_field = n_distinct(paste(stand_id, tree_id_field, sep = "__"))
  )
print(dup_check)

if (dup_check$n_rows != dup_check$n_unique_mls ||
    dup_check$n_rows != dup_check$n_unique_tls ||
    dup_check$n_rows != dup_check$n_unique_field) {
  warning("Final matched table is not strictly 1:1:1 within stand.")
}

readr::write_csv(matched, file.path(paths$out_dir, "matched_field_tls_mls.csv"))

# ============================================================
# Chunk 5.5 — Create final table + error columns
# Description:
# Computes pairwise errors (MLS-TLS, TLS-Field, MLS-Field) for DBH
# and height on the matched dataset.
# ============================================================

final <- matched %>%
  mutate(
    # MLS vs TLS
    err_dbh_mls_tls_cm = dbh_mls_cm - dbh_tls_cm,
    err_h_mls_tls_m    = h_mls_m   - h_tls_m,
    
    # TLS vs Field
    err_dbh_tls_field_cm = dbh_tls_cm - dbh_field_cm,
    err_h_tls_field_m    = h_tls_m    - h_field_m,
    
    # MLS vs Field
    err_dbh_mls_field_cm = dbh_mls_cm - dbh_field_cm,
    err_h_mls_field_m    = h_mls_m    - h_field_m
  )

# ============================================================
# Chunk 5.6 — Bark correction and species-specific volume equations
# Implements:
# 1) Roble bark correction (under-bark DBH estimate),
# 2) Stand-to-species keys,
# 3) Species-specific total volume equations,
# then computes volume and volume error columns.
# ============================================================

# Bark correction (Roble): estimate DBH under bark (DBH_ub) from DBH over bark + total height
# Source: Corvalán et al. (2019) model of (2Ec/DAP) vs z = h/H, evaluated at breast height (h = 1.3 m).

dbh_under_bark_roble_cm <- function(dbh_ob_cm, h_total_m, h_bh_m = 1.3) {
  b0 <- 0.62
  b1 <- -0.99
  b2 <- 3.18
  b3 <- -4.64
  b4 <- 1.86
  
  z <- h_bh_m / h_total_m
  R <- ifelse(is.finite(z), (b0 + b1*z + b2*z^2 + b3*z^3 + b4*z^4), NA_real_)
  
  # double bark thickness in mm, then convert to cm
  dbt_mm <- R * dbh_ob_cm
  dbt_cm <- dbt_mm / 10
  
  dbh_ub <- dbh_ob_cm - dbt_cm
  ifelse(is.finite(dbh_ub) & dbh_ub > 0, dbh_ub, NA_real_)
}

stand_to_species_key <- function(stand_id) {
  dplyr::case_when(
    stand_id == "Oregon"  ~ "PINO_OREGON_Pseudotsuga_menziesii",
    stand_id == "Alerce"  ~ "ALERCE_Fitzroya_cupressoides",
    stand_id == "Roble"   ~ "ROBLE_Nothofagus_obliqua",
    stand_id %in% c("Lawson1", "Lawson2") ~ "LAWSON_Chamaecyparis_lawsoniana",
    TRUE ~ NA_character_
  )
}

volume_total_species <- function(species_key, dbh_cm, h_m) {
  dbh_for_vol_cm <- dplyr::case_when(
    species_key == "ROBLE_Nothofagus_obliqua" ~ dbh_under_bark_roble_cm(dbh_cm, h_m),
    TRUE ~ dbh_cm
  )
  
  d2h <- (dbh_for_vol_cm^2) * h_m
  
  dplyr::case_when(
    species_key == "PINO_OREGON_Pseudotsuga_menziesii" &
      is.finite(d2h) ~ (0.030764 + 0.000028 * d2h),
    
    species_key == "ALERCE_Fitzroya_cupressoides" &
      is.finite(d2h) & d2h > 0 ~ exp(-10.291067 + 0.974113 * log(d2h)),
    
    species_key == "ROBLE_Nothofagus_obliqua" &
      is.finite(d2h) & d2h > 0 ~ exp(-10.51381 + 1.010275 * log(d2h)),
    
    TRUE ~ NA_real_
  )
}

final <- final %>%
  mutate(
    species_key = stand_to_species_key(stand_id),
    
    vol_field_m3 = volume_total_species(species_key, dbh_field_cm, h_field_m),
    vol_tls_m3   = volume_total_species(species_key, dbh_tls_cm,   h_tls_m),
    vol_mls_m3   = volume_total_species(species_key, dbh_mls_cm,   h_mls_m),
    
    err_vol_mls_tls_m3   = vol_mls_m3 - vol_tls_m3,
    err_vol_tls_field_m3 = vol_tls_m3 - vol_field_m3,
    err_vol_mls_field_m3 = vol_mls_m3 - vol_field_m3
  )

# MLS↔TLS-only dataset for all stands (includes Alerce and Lawson)
mls_tls_only <- map_mls_tls %>%
  left_join(
    mls %>% select(stand_id, stem_id_mls = stem_id, dbh_mls_cm, h_mls_m),
    by = c("stand_id", "stem_id_mls")
  ) %>%
  left_join(
    tls %>% select(stand_id, stem_id_tls = stem_id, dbh_tls_cm, h_tls_m),
    by = c("stand_id", "stem_id_tls")
  ) %>%
  mutate(
    species_key = stand_to_species_key(stand_id),
    
    err_dbh_mls_tls_cm = dbh_mls_cm - dbh_tls_cm,
    err_h_mls_tls_m    = h_mls_m    - h_tls_m,
    
    vol_tls_m3 = volume_total_species(species_key, dbh_tls_cm, h_tls_m),
    vol_mls_m3 = volume_total_species(species_key, dbh_mls_cm, h_mls_m),
    err_vol_mls_tls_m3 = vol_mls_m3 - vol_tls_m3
  )


# ============================================================
# Chunk 6 — Metrics tables
# Description:
# Produces stand-level summary metrics (n, mapping distances, bias,
# RMSE, MAE) for Field-based and MLS↔TLS-only comparisons.
# ============================================================

metrics_by_stand_field <- final %>%
  filter(stand_id %in% stands_field_ok) %>%
  group_by(stand_id) %>%
  summarise(
    n = n(),
    
    # Mapping distances
    map_dist_mls_tls_mean_m = mean(dist_m, na.rm = TRUE),
    map_dist_mls_tls_p95_m  = quantile(dist_m, 0.95, na.rm = TRUE, names = FALSE),
    
    map_dist_tls_field_mean_m = mean(dist_m_field, na.rm = TRUE),
    map_dist_tls_field_p95_m  = quantile(dist_m_field, 0.95, na.rm = TRUE, names = FALSE),
    
    # MLS vs TLS
    dbh_mls_tls_bias_cm = bias(err_dbh_mls_tls_cm),
    dbh_mls_tls_rmse_cm = rmse(err_dbh_mls_tls_cm),
    dbh_mls_tls_mae_cm  = mae(err_dbh_mls_tls_cm),
    h_mls_tls_bias_m    = bias(err_h_mls_tls_m),
    h_mls_tls_rmse_m    = rmse(err_h_mls_tls_m),
    h_mls_tls_mae_m     = mae(err_h_mls_tls_m),
    vol_mls_tls_bias_m3 = bias(err_vol_mls_tls_m3),
    vol_mls_tls_rmse_m3 = rmse(err_vol_mls_tls_m3),
    vol_mls_tls_mae_m3  = mae(err_vol_mls_tls_m3),
    
    # TLS vs Field
    dbh_tls_field_bias_cm = bias(err_dbh_tls_field_cm),
    dbh_tls_field_rmse_cm = rmse(err_dbh_tls_field_cm),
    dbh_tls_field_mae_cm  = mae(err_dbh_tls_field_cm),
    h_tls_field_bias_m    = bias(err_h_tls_field_m),
    h_tls_field_rmse_m    = rmse(err_h_tls_field_m),
    h_tls_field_mae_m     = mae(err_h_tls_field_m),
    vol_tls_field_bias_m3 = bias(err_vol_tls_field_m3),
    vol_tls_field_rmse_m3 = rmse(err_vol_tls_field_m3),
    vol_tls_field_mae_m3  = mae(err_vol_tls_field_m3),
    
    # MLS vs Field
    dbh_mls_field_bias_cm = bias(err_dbh_mls_field_cm),
    dbh_mls_field_rmse_cm = rmse(err_dbh_mls_field_cm),
    dbh_mls_field_mae_cm  = mae(err_dbh_mls_field_cm),
    h_mls_field_bias_m    = bias(err_h_mls_field_m),
    h_mls_field_rmse_m    = rmse(err_h_mls_field_m),
    h_mls_field_mae_m     = mae(err_h_mls_field_m),
    vol_mls_field_bias_m3 = bias(err_vol_mls_field_m3),
    vol_mls_field_rmse_m3 = rmse(err_vol_mls_field_m3),
    vol_mls_field_mae_m3  = mae(err_vol_mls_field_m3),
    
    .groups = "drop"
  )

readr::write_csv(metrics_by_stand_field, file.path(paths$out_dir, "metrics_by_stand_field.csv"))

metrics_by_stand_mls_tls <- mls_tls_only %>%
  group_by(stand_id) %>%
  summarise(
    n = n(),
    map_dist_mean_m = mean(dist_m, na.rm = TRUE),
    map_dist_p95_m  = quantile(dist_m, 0.95, na.rm = TRUE, names = FALSE),
    
    dbh_bias_cm = bias(err_dbh_mls_tls_cm),
    dbh_rmse_cm = rmse(err_dbh_mls_tls_cm),
    dbh_mae_cm  = mae(err_dbh_mls_tls_cm),
    
    h_bias_m = bias(err_h_mls_tls_m),
    h_rmse_m = rmse(err_h_mls_tls_m),
    h_mae_m  = mae(err_h_mls_tls_m),
    
    vol_bias_m3 = bias(err_vol_mls_tls_m3),
    vol_rmse_m3 = rmse(err_vol_mls_tls_m3),
    vol_mae_m3  = mae(err_vol_mls_tls_m3),
    
    .groups = "drop"
  )

readr::write_csv(metrics_by_stand_mls_tls, file.path(paths$out_dir, "metrics_by_stand_mls_tls_allstands.csv"))


# ============================================================
# Chunk 7 — Plots
# Defines plotting functions and exports figure files for:
# (1) DBH comparisons,
# (2) Height comparisons,
# (3) Volume comparisons,
# (4) Mapping distance vs absolute DBH error.
# ============================================================

facet_stats <- function(df, x, y, stand = "stand_id") {
  df %>%
    group_by(.data[[stand]]) %>%
    summarise(
      n = sum(is.finite(.data[[x]]) & is.finite(.data[[y]])),
      bias = bias(.data[[y]] - .data[[x]]),
      rmse = rmse(.data[[y]] - .data[[x]]),
      .groups = "drop"
    ) %>%
    mutate(
      facet_lab = paste0(.data[[stand]], "\n",
                         "n=", n,
                         " | bias=", round(bias, 2),
                         " | RMSE=", round(rmse, 2))
    )
}

library(patchwork)

# 5-stand panel with centered layout (3 top, 2 bottom centered)
plot_xy_5stands_centered <- function(df, x, y, xlab, ylab, title,
                                     stand_order = c("Alerce", "Lawson1", "Lawson2", "Oregon", "Roble"),
                                     free_scales = FALSE,
                                     add_lm = TRUE) {
  
  # keep only finite pairs for plotting/stats
  df <- df %>% dplyr::filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  if (nrow(df) == 0) return(NULL)
  
  st <- facet_stats(df, x, y) %>% mutate(stand_id = as.character(stand_id))
  
  lims <- NULL
  if (!free_scales) lims <- range(c(df[[x]], df[[y]]), na.rm = TRUE)
  
  make_one <- function(stand) {
    d <- df %>% dplyr::filter(stand_id == stand)
    lab <- st %>% dplyr::filter(stand_id == stand) %>% dplyr::pull(facet_lab)
    lab <- ifelse(length(lab) == 0, stand, lab)
    
    p <- ggplot(d, aes(x = .data[[x]], y = .data[[y]])) +
      geom_point(alpha = 0.7, size = 2) +
      geom_abline(slope = 1, intercept = 0) +
      labs(title = lab, x = xlab, y = ylab) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(size = 11, hjust = 0.5),
        plot.title.position = "plot"
      )
    
    if (!is.null(lims) && all(is.finite(lims))) {
      p <- p + coord_equal(xlim = lims, ylim = lims, expand = FALSE)
    } else {
      p <- p + coord_equal()
    }
    
    if (add_lm && nrow(d) >= 2) {
      p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8)
    }
    
    p
  }
  
  plots <- purrr::map(stand_order, make_one)
  
  (plots[[1]] | plots[[2]] | plots[[3]]) /
    (plot_spacer() | plots[[4]] | plots[[5]] | plot_spacer()) +
    plot_layout(widths = c(1, 3, 3, 1)) +
    plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold"))
}

# 5-stand panel for mapping distance relationships
plot_mapdist_5stands_centered <- function(df,
                                          x = "dist_m",
                                          y = "abs_err",
                                          xlab = "MLS→TLS mapping distance (m)",
                                          ylab = "|DBH error| (cm)",
                                          title = "Mapping distance vs |DBH error|",
                                          stand_order = c("Alerce", "Lawson1", "Lawson2", "Oregon", "Roble"),
                                          free_scales = TRUE) {
  
  # keep only finite pairs
  df <- df %>% dplyr::filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  if (nrow(df) == 0) return(NULL)
  
  make_one <- function(stand) {
    d <- df %>% dplyr::filter(stand_id == stand)
    
    ggplot(d, aes(x = .data[[x]], y = .data[[y]])) +
      geom_point(alpha = 0.6, size = 2) +
      labs(title = stand, x = xlab, y = ylab) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(size = 11, hjust = 0.5),
        plot.title.position = "plot"
      ) +
      (if (free_scales) NULL else coord_cartesian())
  }
  
  plots <- purrr::map(stand_order, make_one)
  
  (plots[[1]] | plots[[2]] | plots[[3]]) /
    (plot_spacer() | plots[[4]] | plots[[5]] | plot_spacer()) +
    plot_layout(widths = c(1, 3, 3, 1)) +
    plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold"))
}

# Generic facet plot function (for variable number of stands)
# ncol can be set per figure (e.g., 2 for 4 panels, 3 for 3 panels in one row)
plot_xy_facet <- function(df, x, y, xlab, ylab, title,
                          free_scales = FALSE,
                          add_lm = TRUE,
                          ncol = 3) {
  
  df <- df %>% filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
  if (nrow(df) == 0) return(NULL)
  
  st <- facet_stats(df, x, y)
  df2 <- df %>% left_join(st %>% select(stand_id, facet_lab), by = "stand_id")
  
  p <- ggplot(df2, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0) +
    facet_wrap(~facet_lab,
               scales = if (free_scales) "free" else "fixed",
               ncol = ncol) +
    labs(x = xlab, y = ylab, title = title) +
    coord_equal() +
    theme_bw(base_size = 13) +
    theme(
      strip.text = element_text(size = 11),
      plot.title = element_text(face = "bold")
    )
  
  if (add_lm) {
    p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.8)
  }
  
  p
}

# ------------------------------------------------------------
# Figure export section
# ------------------------------------------------------------

# ---- Field vs MLS (stands with field data) ----
df_dbh_field_mls <- final %>%
  filter(stand_id %in% stands_field_ok) %>%
  filter(is.finite(dbh_field_cm), is.finite(dbh_mls_cm))

p_dbh_field_mls <- plot_xy_facet(
  df_dbh_field_mls,
  x = "dbh_field_cm", y = "dbh_mls_cm",
  xlab = "DBH Field (cm)", ylab = "DBH MLS (cm)",
  title = "DBH: Field vs MLS",
  free_scales = FALSE,
  add_lm = TRUE,
  ncol = 2   # 4 stands -> 2x2 layout
)

df_h_field_mls <- final %>%
  filter(stand_id %in% stands_field_ok) %>%
  filter(is.finite(h_field_m), is.finite(h_mls_m))

p_h_field_mls <- plot_xy_facet(
  df_h_field_mls,
  x = "h_field_m", y = "h_mls_m",
  xlab = "Height Field (m)", ylab = "Height MLS (m)",
  title = "Height: Field vs MLS",
  free_scales = FALSE,
  add_lm = TRUE,
  ncol = 2   # 4 stands -> 2x2 layout
)

# Volume exists only for stands with implemented equation
df_v_field_mls <- final %>%
  filter(stand_id %in% stands_field_ok, stand_id %in% stands_with_volume) %>%
  filter(is.finite(vol_field_m3), is.finite(vol_mls_m3))

p_v_field_mls <- plot_xy_facet(
  df_v_field_mls,
  x = "vol_field_m3", y = "vol_mls_m3",
  xlab = "Volume Field (m³)", ylab = "Volume MLS (m³)",
  title = "Volume: Field vs MLS",
  free_scales = FALSE,
  add_lm = TRUE,
  ncol = 3   # Alerce + Oregon + Roble -> single row
)

if (!is.null(p_dbh_field_mls)) {
  ggsave(file.path(paths$out_dir, "FIG_dbh_field_vs_mls.png"),
         p_dbh_field_mls, width = 10, height = 6, dpi = 300)
} else {
  message("Skipping FIG_dbh_field_vs_mls: no valid matched DBH pairs.")
}

if (!is.null(p_h_field_mls)) {
  ggsave(file.path(paths$out_dir, "FIG_height_field_vs_mls.png"),
         p_h_field_mls, width = 10, height = 6, dpi = 300)
} else {
  message("Skipping FIG_height_field_vs_mls: no valid matched height pairs.")
}

if (!is.null(p_v_field_mls)) {
  ggsave(file.path(paths$out_dir, "FIG_volume_field_vs_mls.png"),
         p_v_field_mls, width = 14, height = 5, dpi = 300)  # wider for 3 panels in one row
} else {
  message("Skipping FIG_volume_field_vs_mls: no valid matched volume pairs.")
}

# ---- MLS vs TLS (all stands / volume stands) ----
p_dbh_mls_tls <- plot_xy_5stands_centered(
  mls_tls_only,
  x = "dbh_tls_cm", y = "dbh_mls_cm",
  xlab = "TLS DBH (cm)", ylab = "MLS DBH (cm)",
  title = "DBH: MLS vs TLS",
  free_scales = FALSE,
  add_lm = TRUE
)

p_h_mls_tls <- plot_xy_5stands_centered(
  mls_tls_only,
  x = "h_tls_m", y = "h_mls_m",
  xlab = "TLS Height (m)", ylab = "MLS Height (m)",
  title = "Height: MLS vs TLS",
  free_scales = FALSE,
  add_lm = TRUE
)

p_v_mls_tls <- plot_xy_facet(
  mls_tls_only %>%
    filter(stand_id %in% stands_with_volume) %>%
    filter(is.finite(vol_tls_m3), is.finite(vol_mls_m3)),
  x = "vol_tls_m3", y = "vol_mls_m3",
  xlab = "TLS Volume (m³)", ylab = "MLS Volume (m³)",
  title = "Volume: MLS vs TLS",
  free_scales = FALSE,
  add_lm = TRUE,
  ncol = 3   # Alerce + Oregon + Roble -> single row
)

if (!is.null(p_dbh_mls_tls)) {
  ggsave(file.path(paths$out_dir, "FIG_dbh_mls_vs_tls.png"),
         p_dbh_mls_tls, width = 12, height = 7, dpi = 300)
}

if (!is.null(p_h_mls_tls)) {
  ggsave(file.path(paths$out_dir, "FIG_height_mls_vs_tls.png"),
         p_h_mls_tls, width = 12, height = 7, dpi = 300)
}

if (!is.null(p_v_mls_tls)) {
  ggsave(file.path(paths$out_dir, "FIG_volume_mls_vs_tls.png"),
         p_v_mls_tls, width = 14, height = 5, dpi = 300)  # wider for 3 panels in one row
} else {
  message("Skipping FIG_volume_mls_vs_tls: no valid volume pairs.")
}

# ---- Mapping distance vs absolute DBH error ----
p_mapdist_absdbh <- mls_tls_only %>%
  mutate(abs_err = abs(err_dbh_mls_tls_cm)) %>%
  plot_mapdist_5stands_centered(
    x = "dist_m",
    y = "abs_err",
    xlab = "MLS→TLS mapping distance (m)",
    ylab = "|DBH error| (cm)",
    title = "Mapping distance vs |DBH error|",
    free_scales = TRUE
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

if (!is.null(p_mapdist_absdbh)) {
  ggsave(file.path(paths$out_dir, "FIG_mapdist_vs_absdbherror.png"),
         p_mapdist_absdbh, width = 12, height = 7, dpi = 300)
} else {
  message("Skipping FIG_mapdist_vs_absdbherror: no valid mapping-distance pairs.")
}


# ============================================================
# Chunk 8 — Shapiro tests, paired tests, correlations, and CCC
# Description:
# Performs normality checks on paired differences, selects paired t-test
# or Wilcoxon signed-rank test accordingly, and computes Pearson r and
# Lin's CCC (with CI when available) by stand for each comparison type.
# Exports test tables
# ============================================================

safe_shapiro <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || length(x) > 5000) return(NA_real_)
  tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
}

paired_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(list(test = NA_character_, p_value = NA_real_, shapiro_p = NA_real_))
  }
  
  d <- x - y
  p_shap <- safe_shapiro(d)
  
  if (!is.na(p_shap) && p_shap >= 0.05) {
    tt <- t.test(x, y, paired = TRUE)
    return(list(test = "paired_t", p_value = tt$p.value, shapiro_p = p_shap))
  } else {
    wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
    return(list(test = "wilcoxon", p_value = wt$p.value, shapiro_p = p_shap))
  }
}

pearson_r <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok], method = "pearson")
}

# Robust Lin's CCC extractor across DescTools versions
lin_ccc_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  
  if (length(x) < 3) {
    return(list(ccc = NA_real_, ccc_lwr = NA_real_, ccc_upr = NA_real_))
  }
  
  out <- tryCatch(DescTools::CCC(x, y), error = function(e) NULL)
  if (is.null(out)) {
    return(list(ccc = NA_real_, ccc_lwr = NA_real_, ccc_upr = NA_real_))
  }
  
  rho_obj <- out$rho.c
  
  ccc_est <- NA_real_
  ccc_lwr <- NA_real_
  ccc_upr <- NA_real_
  
  if (is.matrix(rho_obj) || is.data.frame(rho_obj)) {
    rn <- rownames(rho_obj)
    cn <- colnames(rho_obj)
    
    # Try rownames first
    if (!is.null(rn) && any(grepl("est", rn, ignore.case = TRUE))) {
      ccc_est <- as.numeric(rho_obj[grep("est", rn, ignore.case = TRUE)[1], 1])
    }
    if (!is.null(rn) && any(grepl("lower|lwr", rn, ignore.case = TRUE))) {
      ccc_lwr <- as.numeric(rho_obj[grep("lower|lwr", rn, ignore.case = TRUE)[1], 1])
    }
    if (!is.null(rn) && any(grepl("upper|upr", rn, ignore.case = TRUE))) {
      ccc_upr <- as.numeric(rho_obj[grep("upper|upr", rn, ignore.case = TRUE)[1], 1])
    }
    
    # Fallback by colname
    if (is.na(ccc_est) && !is.null(cn) && any(grepl("rho", cn, ignore.case = TRUE))) {
      ccc_est <- as.numeric(rho_obj[1, grep("rho", cn, ignore.case = TRUE)[1]])
    }
    
    # Final fallback: first numeric values
    vals <- suppressWarnings(as.numeric(rho_obj))
    vals <- vals[is.finite(vals)]
    if (is.na(ccc_est) && length(vals) >= 1) ccc_est <- vals[1]
    if (is.na(ccc_lwr) && length(vals) >= 2) ccc_lwr <- vals[2]
    if (is.na(ccc_upr) && length(vals) >= 3) ccc_upr <- vals[3]
    
  } else {
    vals <- suppressWarnings(as.numeric(rho_obj))
    vals <- vals[is.finite(vals)]
    if (length(vals) >= 1) ccc_est <- vals[1]
    if (length(vals) >= 2) ccc_lwr <- vals[2]
    if (length(vals) >= 3) ccc_upr <- vals[3]
  }
  
  list(ccc = ccc_est, ccc_lwr = ccc_lwr, ccc_upr = ccc_upr)
}

pair_stats <- function(df, stand_col, xcol, ycol, label) {
  out <- df %>%
    group_by(.data[[stand_col]]) %>%
    group_modify(~{
      out_test <- paired_test(.x[[xcol]], .x[[ycol]])
      ccc_out  <- lin_ccc_stats(.x[[xcol]], .x[[ycol]])
      
      tibble(
        comparison = label,
        n = sum(is.finite(.x[[xcol]]) & is.finite(.x[[ycol]])),
        test = as.character(out_test$test),
        p_value = as.numeric(out_test$p_value),
        shapiro_p = as.numeric(out_test$shapiro_p),
        pearson_r = as.numeric(pearson_r(.x[[xcol]], .x[[ycol]])),
        ccc     = as.numeric(ccc_out$ccc),
        ccc_lwr = as.numeric(ccc_out$ccc_lwr),
        ccc_upr = as.numeric(ccc_out$ccc_upr)
      )
    }) %>%
    ungroup()
  
  names(out)[names(out) == stand_col] <- "stand_id"
  out
}

# Publication-ready formatting for test tables
format_test_table <- function(df) {
  df %>%
    group_by(comparison) %>%
    mutate(
      p_value_adj = p.adjust(p_value, method = "BH"),
      sig = case_when(
        is.na(p_value_adj) ~ NA_character_,
        p_value_adj < 0.001 ~ "***",
        p_value_adj < 0.01  ~ "**",
        p_value_adj < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    ungroup() %>%
    mutate(
      across(c(p_value, p_value_adj, shapiro_p), ~ signif(.x, 3)),
      across(c(pearson_r, ccc, ccc_lwr, ccc_upr), ~ round(.x, 3))
    ) %>%
    relocate(stand_id, comparison, n, test, shapiro_p, p_value, p_value_adj, sig,
             pearson_r, ccc, ccc_lwr, ccc_upr)
}

tests_mls_tls <- bind_rows(
  pair_stats(mls_tls_only, "stand_id", "dbh_mls_cm", "dbh_tls_cm", "DBH MLS vs TLS"),
  pair_stats(mls_tls_only, "stand_id", "h_mls_m",    "h_tls_m",    "H MLS vs TLS"),
  pair_stats(mls_tls_only %>% filter(stand_id %in% stands_with_volume),
             "stand_id", "vol_mls_m3", "vol_tls_m3", "V MLS vs TLS")
) %>% format_test_table()

readr::write_csv(tests_mls_tls, file.path(paths$out_dir, "tests_mls_vs_tls_allstands.csv"))

tests_mls_field <- bind_rows(
  pair_stats(final %>% filter(stand_id %in% stands_field_ok),
             "stand_id", "dbh_mls_cm", "dbh_field_cm", "DBH MLS vs Field"),
  pair_stats(final %>% filter(stand_id %in% stands_field_ok),
             "stand_id", "h_mls_m",    "h_field_m",    "H MLS vs Field"),
  pair_stats(final %>% filter(stand_id %in% stands_field_ok, stand_id %in% stands_with_volume),
             "stand_id", "vol_mls_m3", "vol_field_m3", "V MLS vs Field")
) %>% format_test_table()

readr::write_csv(tests_mls_field, file.path(paths$out_dir, "tests_mls_vs_field.csv"))

tests_tls_field <- bind_rows(
  pair_stats(final %>% filter(stand_id %in% stands_field_ok),
             "stand_id", "dbh_tls_cm", "dbh_field_cm", "DBH TLS vs Field"),
  pair_stats(final %>% filter(stand_id %in% stands_field_ok),
             "stand_id", "h_tls_m",    "h_field_m",    "H TLS vs Field"),
  pair_stats(final %>% filter(stand_id %in% stands_field_ok, stand_id %in% stands_with_volume),
             "stand_id", "vol_tls_m3", "vol_field_m3", "V TLS vs Field")
) %>% format_test_table()

readr::write_csv(tests_tls_field, file.path(paths$out_dir, "tests_tls_vs_field.csv"))


# ============================================================
# Chunk 9 — Regression tables
# Description:
# Fits simple linear regressions by stand and comparison type,
# exporting one row per model with intercept/slope estimates,
# 95% CI, p-values, R², adjusted R², and model RMSE.
# ============================================================

lm_table <- function(df, y, x, label) {
  df %>%
    group_by(stand_id) %>%
    group_modify(~{
      dat <- .x %>% dplyr::filter(is.finite(.data[[x]]), is.finite(.data[[y]]))
      
      if (nrow(dat) < 3) {
        return(tibble(
          comparison = label,
          n = nrow(dat),
          intercept = NA_real_,
          intercept_se = NA_real_,
          intercept_lwr = NA_real_,
          intercept_upr = NA_real_,
          intercept_p = NA_real_,
          slope = NA_real_,
          slope_se = NA_real_,
          slope_lwr = NA_real_,
          slope_upr = NA_real_,
          slope_p = NA_real_,
          r2 = NA_real_,
          adj_r2 = NA_real_,
          rmse_model = NA_real_
        ))
      }
      
      fit <- lm(reformulate(x, y), data = dat)
      td <- broom::tidy(fit, conf.int = TRUE, conf.level = 0.95)
      gl <- broom::glance(fit)
      
      get_coef <- function(term_name, col_name) {
        v <- td %>%
          dplyr::filter(term == term_name) %>%
          dplyr::pull(.data[[col_name]])
        if (length(v) == 0) NA_real_ else as.numeric(v[1])
      }
      
      pred <- predict(fit, newdata = dat)
      rmse_fit <- sqrt(mean((dat[[y]] - pred)^2, na.rm = TRUE))
      
      tibble(
        comparison = label,
        n = nrow(dat),
        
        intercept     = get_coef("(Intercept)", "estimate"),
        intercept_se  = get_coef("(Intercept)", "std.error"),
        intercept_lwr = get_coef("(Intercept)", "conf.low"),
        intercept_upr = get_coef("(Intercept)", "conf.high"),
        intercept_p   = get_coef("(Intercept)", "p.value"),
        
        slope         = get_coef(x, "estimate"),
        slope_se      = get_coef(x, "std.error"),
        slope_lwr     = get_coef(x, "conf.low"),
        slope_upr     = get_coef(x, "conf.high"),
        slope_p       = get_coef(x, "p.value"),
        
        r2        = as.numeric(gl$r.squared),
        adj_r2    = as.numeric(gl$adj.r.squared),
        rmse_model = as.numeric(rmse_fit)
      )
    }) %>%
    ungroup() %>%
    relocate(stand_id, comparison, n)
}

format_reg_table <- function(df) {
  df %>%
    mutate(
      across(c(intercept, intercept_se, intercept_lwr, intercept_upr,
               slope, slope_se, slope_lwr, slope_upr,
               r2, adj_r2, rmse_model), ~ round(.x, 4)),
      across(c(intercept_p, slope_p), ~ signif(.x, 3))
    )
}

# --- MLS vs TLS regressions ---
reg_mls_tls <- bind_rows(
  lm_table(mls_tls_only, "dbh_mls_cm", "dbh_tls_cm", "DBH MLS~TLS"),
  lm_table(mls_tls_only, "h_mls_m",    "h_tls_m",    "H MLS~TLS"),
  lm_table(mls_tls_only %>% filter(stand_id %in% stands_with_volume),
           "vol_mls_m3", "vol_tls_m3", "V MLS~TLS")
) %>% format_reg_table()

readr::write_csv(reg_mls_tls, file.path(paths$out_dir, "reg_mls_vs_tls_allstands.csv"))

# --- MLS vs Field regressions ---
reg_mls_field <- bind_rows(
  lm_table(final %>% filter(stand_id %in% stands_field_ok),
           "dbh_mls_cm", "dbh_field_cm", "DBH MLS~Field"),
  lm_table(final %>% filter(stand_id %in% stands_field_ok),
           "h_mls_m",    "h_field_m",    "H MLS~Field"),
  lm_table(final %>% filter(stand_id %in% stands_field_ok, stand_id %in% stands_with_volume),
           "vol_mls_m3", "vol_field_m3", "V MLS~Field")
) %>% format_reg_table()

readr::write_csv(reg_mls_field, file.path(paths$out_dir, "reg_mls_vs_field.csv"))

# --- TLS vs Field regressions ---
reg_tls_field <- bind_rows(
  lm_table(final %>% filter(stand_id %in% stands_field_ok),
           "dbh_tls_cm", "dbh_field_cm", "DBH TLS~Field"),
  lm_table(final %>% filter(stand_id %in% stands_field_ok),
           "h_tls_m",    "h_field_m",    "H TLS~Field"),
  lm_table(final %>% filter(stand_id %in% stands_field_ok, stand_id %in% stands_with_volume),
           "vol_tls_m3", "vol_field_m3", "V TLS~Field")
) %>% format_reg_table()

readr::write_csv(reg_tls_field, file.path(paths$out_dir, "reg_tls_vs_field.csv"))

# ============================================================
# Chunk 10 — Monte Carlo volume uncertainty
# Description:
# Propagates uncertainty from MLS↔TLS DBH and height residuals into
# volume estimates using Monte Carlo simulation. Exports:
# (1) final matched table with MC outputs,
# (2) MLS↔TLS volume subset with stand-specific MC outputs,
# (3) stand-level SD table used for simulation.
# ============================================================

sd_dbh_cm <- sd(mls_tls_only$err_dbh_mls_tls_cm, na.rm = TRUE)
sd_h_m    <- sd(mls_tls_only$err_h_mls_tls_m,    na.rm = TRUE)
n_mc      <- 10000

mc_volume <- function(species_key, dbh_cm, h_m, sd_dbh_cm, sd_h_m, n_mc = 10000) {
  if (!is.finite(dbh_cm) || !is.finite(h_m)) {
    return(c(mean = NA_real_, lwr = NA_real_, upr = NA_real_))
  }

  dbh_s <- pmax(rnorm(n_mc, dbh_cm, sd_dbh_cm), 0.1)
  h_s   <- pmax(rnorm(n_mc, h_m, sd_h_m), 0.1)
  v_s   <- volume_total_species(species_key, dbh_s, h_s)

  c(
    mean = mean(v_s, na.rm = TRUE),
    lwr  = as.numeric(quantile(v_s, 0.025, na.rm = TRUE)),
    upr  = as.numeric(quantile(v_s, 0.975, na.rm = TRUE))
  )
}

# ------------------------------------------------------------
# Sensitivity check: add an optional "allometric model residual"
# When Drake (2003) does not report SEE/ECM for the specific equation,
# we bracket plausible residual sizes and re-run the MC.
#
# For log–log models (Alerce, Roble): multiplicative error via ln(V) residual:
#   V* = V_det * exp(eps), eps ~ N(0, sigma_lnV^2)
#
# For linear model (Oregon): additive error scaled to the predicted volume:
#   V* = V_det + eps, eps ~ N(0, (sigma_relV * |V_det|)^2)
# ------------------------------------------------------------
is_log_model_species <- function(species_key) {
  species_key %in% c("ALERCE_Fitzroya_cupressoides", "ROBLE_Nothofagus_obliqua")
}

mc_volume_sens <- function(species_key, dbh_cm, h_m,
                           sd_dbh_cm, sd_h_m,
                           sigma_lnV = 0,
                           sigma_relV_lin = 0,
                           n_mc = 10000) {
  if (!is.finite(dbh_cm) || !is.finite(h_m)) {
    return(c(mean = NA_real_, lwr = NA_real_, upr = NA_real_))
  }

  dbh_s <- pmax(rnorm(n_mc, dbh_cm, sd_dbh_cm), 0.1)
  h_s   <- pmax(rnorm(n_mc, h_m, sd_h_m), 0.1)

  v_det <- volume_total_species(species_key, dbh_s, h_s)

  if (is_log_model_species(species_key)) {
    # residual in ln-space -> multiplicative in original space
    v_s <- v_det * exp(rnorm(n_mc, mean = 0, sd = sigma_lnV))
  } else {
    # Oregon linear: relative additive error around prediction
    v_s <- v_det + rnorm(n_mc, mean = 0, sd = sigma_relV_lin * abs(v_det))
  }

  # keep volumes non-negative
  v_s <- ifelse(is.finite(v_s) & v_s >= 0, v_s, NA_real_)

  c(
    mean = mean(v_s, na.rm = TRUE),
    lwr  = as.numeric(quantile(v_s, 0.025, na.rm = TRUE)),
    upr  = as.numeric(quantile(v_s, 0.975, na.rm = TRUE))
  )
}

# Sensitivity scenarios (interpretation):
# sigma_lnV = 0.10 ~ ~10% multiplicative scatter (approx, per 1 SD in ln-space)
# sigma_lnV = 0.20, 0.30 widen uncertainty accordingly
sens_scenarios <- tibble::tibble(
  scenario = c("meas_only", "model_0.10", "model_0.20", "model_0.30"),
  sigma_lnV = c(0, 0.10, 0.20, 0.30),
  sigma_relV_lin = c(0, 0.10, 0.20, 0.30)
)

# Stand-specific SDs from MLS–TLS differences (used in stand-level MC below)
 #from MLS–TLS differences (used in stand-level MC below)
sd_by_stand <- mls_tls_only %>%
  group_by(stand_id) %>%
  summarise(
    sd_dbh_cm = sd(err_dbh_mls_tls_cm, na.rm = TRUE),
    sd_h_m    = sd(err_h_mls_tls_m,    na.rm = TRUE),
    n_pairs   = sum(is.finite(err_dbh_mls_tls_cm) & is.finite(err_h_mls_tls_m)),
    .groups = "drop"
  )

readr::write_csv(sd_by_stand, file.path(paths$out_dir, "sd_by_stand_mls_tls.csv"))

# MC for final matched table (global SDs)
mc_final <- final %>%
  mutate(
    mc_field = purrr::pmap(list(species_key, dbh_field_cm, h_field_m),
                           ~mc_volume(..1, ..2, ..3, sd_dbh_cm, sd_h_m, n_mc)),
    mc_tls   = purrr::pmap(list(species_key, dbh_tls_cm,   h_tls_m),
                           ~mc_volume(..1, ..2, ..3, sd_dbh_cm, sd_h_m, n_mc)),
    mc_mls   = purrr::pmap(list(species_key, dbh_mls_cm,   h_mls_m),
                           ~mc_volume(..1, ..2, ..3, sd_dbh_cm, sd_h_m, n_mc))
  ) %>%
  mutate(
    vol_field_mc_mean = purrr::map_dbl(mc_field, ~ .x[["mean"]]),
    vol_field_mc_lwr  = purrr::map_dbl(mc_field, ~ .x[["lwr"]]),
    vol_field_mc_upr  = purrr::map_dbl(mc_field, ~ .x[["upr"]]),
    
    vol_tls_mc_mean   = purrr::map_dbl(mc_tls, ~ .x[["mean"]]),
    vol_tls_mc_lwr    = purrr::map_dbl(mc_tls, ~ .x[["lwr"]]),
    vol_tls_mc_upr    = purrr::map_dbl(mc_tls, ~ .x[["upr"]]),
    
    vol_mls_mc_mean   = purrr::map_dbl(mc_mls, ~ .x[["mean"]]),
    vol_mls_mc_lwr    = purrr::map_dbl(mc_mls, ~ .x[["lwr"]]),
    vol_mls_mc_upr    = purrr::map_dbl(mc_mls, ~ .x[["upr"]])
  ) %>%
  select(-mc_field, -mc_tls, -mc_mls)

# MC for MLS↔TLS volume subset (stand-specific SDs)
mc_mls_tls <- mls_tls_only %>%
  filter(stand_id %in% stands_with_volume) %>%
  left_join(sd_by_stand, by = "stand_id") %>%
  mutate(
    mc_tls = purrr::pmap(list(species_key, dbh_tls_cm, h_tls_m, sd_dbh_cm, sd_h_m),
                         ~mc_volume(..1, ..2, ..3, ..4, ..5, n_mc)),
    mc_mls = purrr::pmap(list(species_key, dbh_mls_cm, h_mls_m, sd_dbh_cm, sd_h_m),
                         ~mc_volume(..1, ..2, ..3, ..4, ..5, n_mc))
  ) %>%
  mutate(
    vol_tls_mc_mean = purrr::map_dbl(mc_tls, ~ .x[["mean"]]),
    vol_tls_mc_lwr  = purrr::map_dbl(mc_tls, ~ .x[["lwr"]]),
    vol_tls_mc_upr  = purrr::map_dbl(mc_tls, ~ .x[["upr"]]),
    
    vol_mls_mc_mean = purrr::map_dbl(mc_mls, ~ .x[["mean"]]),
    vol_mls_mc_lwr  = purrr::map_dbl(mc_mls, ~ .x[["lwr"]]),
    vol_mls_mc_upr  = purrr::map_dbl(mc_mls, ~ .x[["upr"]])
  ) %>%
  select(-mc_tls, -mc_mls)

# Export MC outputs
readr::write_csv(mc_mls_tls, file.path(paths$out_dir, "mls_tls_with_volume_montecarlo_standSD3.csv"))
readr::write_csv(mc_final,   file.path(paths$out_dir, "final_with_volume_montecarlo.csv"))

# ============================================================
# Chunk 10.1 — Sensitivity check (optional allometric residual)
# Description:
# Re-runs MC under plausible "model residual" scenarios because
# Drake (2003) does not report SEE/ECM for these specific equations.
# Outputs LONG tables (one row per tree × method × scenario).
# ============================================================

run_mc_sensitivity_final <- function(df, sd_dbh_cm, sd_h_m, n_mc, scenarios_tbl) {
  # Returns long format: one row per tree × method × scenario
  purrr::map_dfr(seq_len(nrow(scenarios_tbl)), function(k) {
    sc <- scenarios_tbl[k, ]

    out <- df %>%
      mutate(
        mc_field = purrr::pmap(list(species_key, dbh_field_cm, h_field_m),
                               ~mc_volume_sens(..1, ..2, ..3, sd_dbh_cm, sd_h_m,
                                               sigma_lnV = sc$sigma_lnV,
                                               sigma_relV_lin = sc$sigma_relV_lin,
                                               n_mc = n_mc)),
        mc_tls   = purrr::pmap(list(species_key, dbh_tls_cm,   h_tls_m),
                               ~mc_volume_sens(..1, ..2, ..3, sd_dbh_cm, sd_h_m,
                                               sigma_lnV = sc$sigma_lnV,
                                               sigma_relV_lin = sc$sigma_relV_lin,
                                               n_mc = n_mc)),
        mc_mls   = purrr::pmap(list(species_key, dbh_mls_cm,   h_mls_m),
                               ~mc_volume_sens(..1, ..2, ..3, sd_dbh_cm, sd_h_m,
                                               sigma_lnV = sc$sigma_lnV,
                                               sigma_relV_lin = sc$sigma_relV_lin,
                                               n_mc = n_mc))
      ) %>%
      transmute(
        stand_id,
        species_key,
        scenario = sc$scenario,
        # keep IDs if present
        tree_id_field = if ("tree_id_field" %in% names(df)) as.character(tree_id_field) else NA_character_,
        stem_id_tls   = if ("stem_id_tls" %in% names(df)) as.character(stem_id_tls) else NA_character_,
        stem_id_mls   = if ("stem_id_mls" %in% names(df)) as.character(stem_id_mls) else NA_character_,
        # deterministic volumes (for reference)
        vol_field_m3 = vol_field_m3,
        vol_tls_m3   = vol_tls_m3,
        vol_mls_m3   = vol_mls_m3,
        # MC results
        field_mc_mean = purrr::map_dbl(mc_field, ~ .x[["mean"]]),
        field_mc_lwr  = purrr::map_dbl(mc_field, ~ .x[["lwr"]]),
        field_mc_upr  = purrr::map_dbl(mc_field, ~ .x[["upr"]]),

        tls_mc_mean   = purrr::map_dbl(mc_tls, ~ .x[["mean"]]),
        tls_mc_lwr    = purrr::map_dbl(mc_tls, ~ .x[["lwr"]]),
        tls_mc_upr    = purrr::map_dbl(mc_tls, ~ .x[["upr"]]),

        mls_mc_mean   = purrr::map_dbl(mc_mls, ~ .x[["mean"]]),
        mls_mc_lwr    = purrr::map_dbl(mc_mls, ~ .x[["lwr"]]),
        mls_mc_upr    = purrr::map_dbl(mc_mls, ~ .x[["upr"]])
      ) %>%
      tidyr::pivot_longer(
        cols = c(field_mc_mean, field_mc_lwr, field_mc_upr,
                 tls_mc_mean, tls_mc_lwr, tls_mc_upr,
                 mls_mc_mean, mls_mc_lwr, mls_mc_upr),
        names_to = c("method", ".value"),
        names_pattern = "^(field|tls|mls)_mc_(mean|lwr|upr)$"
      )

    out
  })
}

run_mc_sensitivity_mls_tls <- function(df_mls_tls, n_mc, scenarios_tbl) {
  purrr::map_dfr(seq_len(nrow(scenarios_tbl)), function(k) {
    sc <- scenarios_tbl[k, ]

    out <- df_mls_tls %>%
      mutate(
        mc_tls = purrr::pmap(list(species_key, dbh_tls_cm, h_tls_m, sd_dbh_cm, sd_h_m),
                             ~mc_volume_sens(..1, ..2, ..3, ..4, ..5,
                                             sigma_lnV = sc$sigma_lnV,
                                             sigma_relV_lin = sc$sigma_relV_lin,
                                             n_mc = n_mc)),
        mc_mls = purrr::pmap(list(species_key, dbh_mls_cm, h_mls_m, sd_dbh_cm, sd_h_m),
                             ~mc_volume_sens(..1, ..2, ..3, ..4, ..5,
                                             sigma_lnV = sc$sigma_lnV,
                                             sigma_relV_lin = sc$sigma_relV_lin,
                                             n_mc = n_mc))
      ) %>%
      transmute(
        stand_id,
        species_key,
        scenario = sc$scenario,
        stem_id_tls = as.character(stem_id_tls),
        stem_id_mls = as.character(stem_id_mls),
        vol_tls_m3  = vol_tls_m3,
        vol_mls_m3  = vol_mls_m3,
        tls_mc_mean = purrr::map_dbl(mc_tls, ~ .x[["mean"]]),
        tls_mc_lwr  = purrr::map_dbl(mc_tls, ~ .x[["lwr"]]),
        tls_mc_upr  = purrr::map_dbl(mc_tls, ~ .x[["upr"]]),
        mls_mc_mean = purrr::map_dbl(mc_mls, ~ .x[["mean"]]),
        mls_mc_lwr  = purrr::map_dbl(mc_mls, ~ .x[["lwr"]]),
        mls_mc_upr  = purrr::map_dbl(mc_mls, ~ .x[["upr"]])
      ) %>%
      tidyr::pivot_longer(
        cols = c(tls_mc_mean, tls_mc_lwr, tls_mc_upr,
                 mls_mc_mean, mls_mc_lwr, mls_mc_upr),
        names_to = c("method", ".value"),
        names_pattern = "^(tls|mls)_mc_(mean|lwr|upr)$"
      )

    out
  })
}

# --- Sensitivity MC outputs ---
mc_final_sens_long <- run_mc_sensitivity_final(
  final %>% dplyr::filter(stand_id %in% stands_with_volume),
  sd_dbh_cm, sd_h_m, n_mc, sens_scenarios
)
readr::write_csv(mc_final_sens_long, file.path(paths$out_dir, "final_with_volume_montecarlo_sensitivity_LONG.csv"))

mc_mls_tls_sens_long <- mls_tls_only %>%
  filter(stand_id %in% stands_with_volume) %>%
  left_join(sd_by_stand, by = "stand_id") %>%
  run_mc_sensitivity_mls_tls(n_mc, sens_scenarios)

readr::write_csv(mc_mls_tls_sens_long, file.path(paths$out_dir, "mls_tls_with_volume_montecarlo_sensitivity_LONG.csv"))

# Quick stand-level summary of sensitivity (CI width in m3)
sens_summary <- mc_mls_tls_sens_long %>%
  mutate(ci_width = upr - lwr) %>%
  group_by(stand_id, scenario, method) %>%
  summarise(
    n = n(),
    mean_ci_width_m3 = mean(ci_width, na.rm = TRUE),
    median_ci_width_m3 = median(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(sens_summary, file.path(paths$out_dir, "mc_sensitivity_summary_ciwidth.csv"))

#Plots

# ============================================================
# FIGURES — Sensitivity analysis
# ============================================================

sens_summary_path <- file.path(paths$out_dir, "mc_sensitivity_summary_ciwidth.csv")
sens_long_path    <- file.path(paths$out_dir, "mls_tls_with_volume_montecarlo_sensitivity_LONG.csv")

sens_summary <- readr::read_csv(sens_summary_path, show_col_types = FALSE)
sens_long    <- readr::read_csv(sens_long_path,    show_col_types = FALSE)

# --- Consistent ordering/labels ---
scenario_levels <- c("meas_only", "model_0.10", "model_0.20", "model_0.30")
scenario_labels <- c("0", "0.10", "0.20", "0.30")

stand_levels <- c("Alerce", "Oregon", "Roble")

sens_summary <- sens_summary %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels, labels = scenario_labels),
    method   = factor(toupper(method), levels = c("TLS", "MLS")),
    stand_id = factor(stand_id, levels = stand_levels)
  )

sens_long2 <- sens_long %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels, labels = scenario_labels),
    method   = factor(toupper(method), levels = c("TLS", "MLS")),
    stand_id = factor(stand_id, levels = stand_levels),
    ci_width = upr - lwr
  ) %>%
  filter(is.finite(ci_width))

# --- Color palette
method_cols <- c("TLS" = "black", "MLS" = "dodgerblue3")

# --- Thesis theme: 
theme_thesis <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = base_size),
      strip.background = element_rect(fill = "grey80", colour = "grey40", linewidth = 0.6),
      strip.text = element_text(size = base_size, face = "plain"),
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_line(colour = "grey95")
    )
}

# ============================================================
# 1) MAIN FIGURE: line plot of CI width vs σ (median)
# ============================================================
use_metric <- "median_ci_width_m3"
y_lab <- "Median 95% CI width (m³)"

p_sens_line <- ggplot(
  sens_summary,
  aes(x = scenario, y = .data[[use_metric]],
      group = method, colour = method)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = method_cols) +
  facet_wrap(~ stand_id, nrow = 1, scales = "free_y") +
  labs(
    title = "Sensitivity of volume uncertainty to added allometric residual error",
    x = expression("Allometric residual scenario ("*sigma*")"),
    y = y_lab,
    colour = NULL
  ) +
  theme_thesis()

ggsave(
  filename = file.path(paths$out_dir, "FIG_sensitivity_ciwidth_line_median.png"),
  plot = p_sens_line,
  width = 12, height = 4.2, dpi = 300
)

# Optional: mean CI width version, same style
p_sens_line_mean <- p_sens_line +
  aes(y = mean_ci_width_m3) +
  labs(y = "Mean 95% CI width (m³)")

ggsave(
  filename = file.path(paths$out_dir, "FIG_sensitivity_ciwidth_line_mean.png"),
  plot = p_sens_line_mean,
  width = 12, height = 4.2, dpi = 300
)

# ============================================================
# 2) Boxplots of per-tree CI width by scenario
# # ============================================================

p_sens_box <- ggplot(
  sens_long2,
  aes(x = scenario, y = ci_width, colour = method)
) +
  geom_boxplot(
    outlier.alpha = 0.25,
    fill = "grey90",
    linewidth = 0.8
  ) +
  scale_colour_manual(values = method_cols) +
  facet_grid(method ~ stand_id, scales = "free_y") +
  labs(
    title = "Per-tree uncertainty in volume under increasing allometric residual scenarios",
    x = expression("Allometric residual scenario ("*sigma*")"),
    y = "Per-tree 95% CI width (m³)",
    colour = NULL
  ) +
  theme_thesis()

ggsave(
  filename = file.path(paths$out_dir, "FIG_sensitivity_ciwidth_boxplots.png"),
  plot = p_sens_box,
  width = 12, height = 6.5, dpi = 300
)

# print to viewer (optional)
print(p_sens_line)
print(p_sens_box)

# ============================================================
# FINAL NOTES
# ============================================================
# - MLS↔TLS matching was performed spatially (XY) using a Hungarian assignment
#   (global 1:1 optimization) with a 2.0 m maximum distance threshold.
# - TLS↔Field matching was performed spatially (XY) using a Hungarian assignment
#   (global 1:1 optimization) with a 5.0 m maximum distance threshold,
#   after stand-specific species filtering of the Field data.
# - Coordinate QA/QC was performed before TLS↔Field matching; a Lawson field subset
#   had X/Y columns swapped and was corrected to projected XY (Easting/Northing).
# - Figures are exported at 300 dpi with consistent styling and facet labels (n, bias, RMSE).
# - 1:1 lines in comparison plots are visually correct because coord_equal() is used.
# - Volume comparisons involving Field data exclude Lawson stands because no volume equation
#   was implemented for Lawson in this workflow.
# - Statistical test tables include BH-adjusted p-values and Lin's CCC (with CI columns when available).
# - Regression tables are exported in publication-ready format (one row per model).
# - Monte Carlo volume uncertainty outputs include mean, lower (2.5%), and upper (97.5%) CI columns.
# Done.

