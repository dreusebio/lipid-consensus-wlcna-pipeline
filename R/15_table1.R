# R/15_table1.R
# -------------------------------------------------------
# Table 1 — Sample characteristics stratified by APO and PPWR
#
# Layout: Option C — three column sets:
#   Overall | APO vs No APO | Excessive PPWR vs No PPWR
#
# Statistics:
#   Continuous  → Mean ± SD, Mann-Whitney U p-value
#   Categorical → n (%), Fisher's exact p-value
#
# Variables controlled via config/config.yml:
#   table1.continuous_vars  — list of continuous variable names
#   table1.categorical_vars — list of categorical variable names
#   table1.var_labels       — display labels for each variable
#
# Outputs → results/15_table1/
#   Table1.docx  — Word document for manuscript
#   Table1.xlsx  — Excel for supplementary
# -------------------------------------------------------

source("R/00_load_packages.R")

library(dplyr)
library(openxlsx)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Config ────────────────────────────────────────────────
cfg <- yaml::read_yaml("config/config.yml")

out <- file.path("results", "15_table1")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────
demo <- readRDS(file.path(cfg$output$processed,
                           "demographic_data_bmi.rds"))

# race_eth_new and mother_s_education_base are already in demo
# confirmed from colnames(demo) — no merge needed
message("Demographic data loaded: ", nrow(demo), " participants, ",
        ncol(demo), " variables")

# ── Variable definitions from config ─────────────────────
# Pull from config; fall back to sensible defaults if not set

t1_cfg <- cfg$table1 %||% list()

continuous_vars <- t1_cfg$continuous_vars %||% c(
  "age", "bmi_base", "gwg", "ppwr"
)

categorical_vars <- t1_cfg$categorical_vars %||% c(
  "race_eth_new", "mother_s_education_base", "group",
  "apo", "apo_hdp", "apo_gdm", "preterm", "apo_other",
  "ppwr_e"
)

# Display labels — override in config or edit here
var_labels_cfg <- t1_cfg$var_labels %||% list()

var_labels_default <- c(
  age                      = "Age, years",
  bmi_base                 = "Pre-pregnancy BMI, kg/m\u00b2",
  gwg                      = "Gestational weight gain, kg",
  ppwr                     = "Postpartum weight retention, kg",
  race_eth_new             = "Race/ethnicity",
  mother_s_education_base  = "Highest education",
  group                    = "Study arm",
  apo                      = "Adverse pregnancy outcome (composite)",
  apo_hdp                  = "Hypertensive disorder of pregnancy",
  apo_gdm                  = "Gestational diabetes mellitus",
  preterm                  = "Preterm/early-term birth (\u226437 weeks)",
  apo_other                = "Other APO",
  ppwr_e                   = "Excessive PPWR (>105% preconception weight)"
)

# Merge config labels over defaults
var_labels <- var_labels_default
for (nm in names(var_labels_cfg)) {
  var_labels[nm] <- var_labels_cfg[[nm]]
}

get_label <- function(v) {
  unname(var_labels[v] %||% v)
}

# ── Race/ethnicity label map ──────────────────────────────
# Factor levels from traits_raw coding:
#   0 = Multi-racial/Other, 1 = Hispanic, 2 = Asian,
#   3 = Black, 4 = White
RACE_LABELS <- c(
  "0" = "Multi-racial/Other",
  "1" = "Hispanic",
  "2" = "Asian",
  "3" = "Black",
  "4" = "White"
)

# Education label map:
#   0 = High School, 1 = Post-Baccalaureate, 2 = Some College
EDU_LABELS <- c(
  "0" = "High school",
  "1" = "Post-baccalaureate",
  "2" = "Some college"
)

GROUP_LABELS <- c("0" = "Control", "1" = "Intervention")
BINARY_LABELS <- c("0" = "No", "1" = "Yes")

# ── Grouping vectors ──────────────────────────────────────
apo_vec  <- demo$apo   # 0/1
ppwr_vec <- demo$ppwr_e # 0/1

# ── Statistical helpers ───────────────────────────────────

# Continuous: Mean (SD), Mann-Whitney p
summarise_continuous <- function(x, group = NULL) {

  x <- as.numeric(x)

  # Overall
  n_obs    <- sum(!is.na(x))
  mean_val <- mean(x, na.rm = TRUE)
  sd_val   <- sd(x,   na.rm = TRUE)
  # Median (IQR)
    median_val <- median(x, na.rm = TRUE)
    q1_val     <- quantile(x, 0.25, na.rm = TRUE)
    q3_val     <- quantile(x, 0.75, na.rm = TRUE)
    overall    <- sprintf("%.2f (%.2f, %.2f)", median_val, q1_val, q3_val)

  if (is.null(group)) return(list(overall = overall, p = NA_character_))

  g     <- as.character(group)
  lvls  <- sort(unique(g[!is.na(g)]))

  cells <- sapply(lvls, function(lv) {
    xi     <- x[!is.na(g) & g == lv]
    med_i  <- median(xi, na.rm = TRUE)
    q1_i   <- quantile(xi, 0.25, na.rm = TRUE)
    q3_i   <- quantile(xi, 0.75, na.rm = TRUE)
    sprintf("%.2f (%.2f, %.2f)", med_i, q1_i, q3_i)
  })

  # Mann-Whitney U
  if (length(lvls) == 2) {
    x0 <- x[!is.na(g) & g == lvls[1]]
    x1 <- x[!is.na(g) & g == lvls[2]]
    p  <- tryCatch(
      wilcox.test(x0, x1, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    p_str <- format_p(p)
  } else {
    p_str <- NA_character_
  }

  list(overall = overall, cells = cells, p = p_str)
}

# Categorical: n (%), Fisher's exact p
summarise_categorical <- function(x, group = NULL,
                                   level_labels = NULL) {

  x <- as.character(x)
  x[x %in% c("NA", "NaN", "")] <- NA_character_

  lvls <- sort(unique(x[!is.na(x)]))
  N    <- sum(!is.na(x))

  if (!is.null(level_labels)) {
    display_lvls <- ifelse(lvls %in% names(level_labels),
                           level_labels[lvls], lvls)
  } else {
    display_lvls <- lvls
  }

  rows <- lapply(seq_along(lvls), function(i) {
    lv  <- lvls[i]
    dlv <- display_lvls[i]
    n_i <- sum(!is.na(x) & x == lv)
    pct <- if (N > 0) 100 * n_i / N else 0
    overall_cell <- sprintf("%d (%.1f%%)", n_i, pct)

    if (is.null(group)) {
      return(list(level = dlv, overall = overall_cell,
                  cells = NULL, p = NA_character_))
    }

    g    <- as.character(group)
    g[g %in% c("NA","NaN","")] <- NA_character_
    glvls <- sort(unique(g[!is.na(g)]))

    cells <- sapply(glvls, function(gl) {
      idx  <- !is.na(g) & g == gl
      n_g  <- sum(!is.na(x[idx]) & x[idx] == lv)
      N_g  <- sum(!is.na(x[idx]))
      pct_g <- if (N_g > 0) 100 * n_g / N_g else 0
      sprintf("%d (%.1f%%)", n_g, pct_g)
    })

    list(level = dlv, overall = overall_cell,
         cells = cells, p = NA_character_)
  })

  # Fisher's exact on full contingency table
  if (!is.null(group)) {
    g <- as.character(group)
    g[g %in% c("NA","NaN","")] <- NA_character_
    tbl <- table(x, g)
    p <- tryCatch(
      fisher.test(tbl, simulate.p.value = TRUE,
                  B = 10000)$p.value,
      error = function(e) NA_real_
    )
    # Assign p to first level row only (displayed once per variable)
    rows[[1]]$p <- format_p(p)
  }

  rows
}

# Format p-value
format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# ── Build table rows ──────────────────────────────────────
# Returns a data frame with columns:
#   Variable | Level | Overall |
#   APO_No | APO_Yes | p_APO |
#   PPWR_No | PPWR_Yes | p_PPWR

build_table <- function(data, cont_vars, cat_vars,
                         apo_group, ppwr_group) {

  rows <- list()

  # Header row: N per group
  n_total    <- nrow(data)
  n_apo_no   <- sum(!is.na(apo_group)  & apo_group  == 0)
  n_apo_yes  <- sum(!is.na(apo_group)  & apo_group  == 1)
  n_ppwr_no  <- sum(!is.na(ppwr_group) & ppwr_group == 0)
  n_ppwr_yes <- sum(!is.na(ppwr_group) & ppwr_group == 1)

  rows[["n"]] <- data.frame(
    Variable = "N",
    Level    = "",
    Overall  = as.character(n_total),
    APO_No   = as.character(n_apo_no),
    APO_Yes  = as.character(n_apo_yes),
    p_APO    = "",
    PPWR_No  = as.character(n_ppwr_no),
    PPWR_Yes = as.character(n_ppwr_yes),
    p_PPWR   = "",
    stringsAsFactors = FALSE
  )

  # Continuous variables
  for (v in cont_vars) {
    if (!v %in% colnames(data)) {
      message("  Skipping continuous var not found: ", v)
      next
    }
    x    <- data[[v]]
    lbl  <- get_label(v)

    s_apo  <- summarise_continuous(x, apo_group)
    s_ppwr <- summarise_continuous(x, ppwr_group)

    rows[[v]] <- data.frame(
      Variable = lbl,
      Level = "Median (Q1, Q3)",
      Overall  = s_apo$overall,
      APO_No   = s_apo$cells["0"],
      APO_Yes  = s_apo$cells["1"],
      p_APO    = s_apo$p,
      PPWR_No  = s_ppwr$cells["0"],
      PPWR_Yes = s_ppwr$cells["1"],
      p_PPWR   = s_ppwr$p,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  # Categorical variables
  for (v in cat_vars) {
    if (!v %in% colnames(data)) {
      message("  Skipping categorical var not found: ", v)
      next
    }
    x   <- data[[v]]
    lbl <- get_label(v)

    # Choose level labels
    lvl_map <- if (v == "race_eth_new") RACE_LABELS else
               if (v == "mother_s_education_base") EDU_LABELS else
               if (v == "group") GROUP_LABELS else
               BINARY_LABELS

    s_apo  <- summarise_categorical(x, apo_group,  lvl_map)
    s_ppwr <- summarise_categorical(x, ppwr_group, lvl_map)

    # Header row for variable (no values, just label + p-value)
    header_row <- data.frame(
      Variable = lbl,
      Level    = "",
      Overall  = "",
      APO_No   = "",
      APO_Yes  = "",
      p_APO    = s_apo[[1]]$p,
      PPWR_No  = "",
      PPWR_Yes = "",
      p_PPWR   = s_ppwr[[1]]$p,
      stringsAsFactors = FALSE
    )

    level_rows <- lapply(seq_along(s_apo), function(i) {
      data.frame(
        Variable = "",
        Level    = paste0("  ", s_apo[[i]]$level),
        Overall  = s_apo[[i]]$overall,
        APO_No   = s_apo[[i]]$cells["0"]   %||% "",
        APO_Yes  = s_apo[[i]]$cells["1"]   %||% "",
        p_APO    = "",
        PPWR_No  = s_ppwr[[i]]$cells["0"]  %||% "",
        PPWR_Yes = s_ppwr[[i]]$cells["1"]  %||% "",
        p_PPWR   = "",
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    })

    rows[[v]] <- dplyr::bind_rows(header_row, level_rows)
  }

  dplyr::bind_rows(rows)
}

# ── Run ───────────────────────────────────────────────────
message("Building Table 1...")

tbl <- build_table(
  data       = demo,
  cont_vars  = continuous_vars,
  cat_vars   = categorical_vars,
  apo_group  = apo_vec,
  ppwr_group = ppwr_vec
)

# Clean up NAs in display
tbl[is.na(tbl)] <- ""

# Rename columns for display
colnames(tbl) <- c(
  "Characteristic", "Level",
  sprintf("Overall\n(N=%d)", nrow(demo)),
  sprintf("No APO\n(n=%d)", sum(!is.na(apo_vec) & apo_vec == 0)),
  sprintf("APO\n(n=%d)",    sum(!is.na(apo_vec) & apo_vec == 1)),
  "p-value\n(APO)",
  sprintf("No Excessive PPWR\n(n=%d)", sum(!is.na(ppwr_vec) & ppwr_vec == 0)),
  sprintf("Excessive PPWR\n(n=%d)",    sum(!is.na(ppwr_vec) & ppwr_vec == 1)),
  "p-value\n(PPWR)"
)

message("  Table rows: ", nrow(tbl))
print(tbl)

# ── Save Excel ────────────────────────────────────────────
message("\nSaving Table1.xlsx...")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Table 1")

# Styles
header_style <- openxlsx::createStyle(
  fontName      = "Helvetica",
  fontSize      = 10,
  fontColour    = "#FFFFFF",
  fgFill        = "#2C3E50",
  halign        = "CENTER",
  valign        = "CENTER",
  textDecoration = "bold",
  wrapText      = TRUE,
  border        = "Bottom",
  borderColour  = "#FFFFFF"
)
body_style <- openxlsx::createStyle(
  fontName  = "Helvetica",
  fontSize  = 10,
  valign    = "CENTER",
  wrapText  = FALSE
)
bold_style <- openxlsx::createStyle(
  fontName       = "Helvetica",
  fontSize       = 10,
  textDecoration = "bold"
)
pval_style <- openxlsx::createStyle(
  fontName  = "Helvetica",
  fontSize  = 10,
  halign    = "CENTER"
)
shaded_style <- openxlsx::createStyle(
  fontName = "Helvetica",
  fontSize = 10,
  fgFill   = "#F7F9FB"
)

# Write data
openxlsx::writeData(wb, "Table 1", tbl,
                    startRow = 1, startCol = 1,
                    headerStyle = header_style,
                    borders = "all",
                    borderStyle = "thin")

# Apply body styles
n_rows <- nrow(tbl)
n_cols <- ncol(tbl)

openxlsx::addStyle(wb, "Table 1", body_style,
                   rows = 2:(n_rows + 1),
                   cols = 1:n_cols, gridExpand = TRUE)

# Bold variable name rows (where Level is empty, Variable is not)
for (i in seq_len(n_rows)) {
  if (tbl[i, 2] == "" && tbl[i, 1] != "") {
    openxlsx::addStyle(wb, "Table 1", bold_style,
                       rows = i + 1, cols = 1)
  }
}

# Shade alternate variable blocks
openxlsx::addStyle(wb, "Table 1", pval_style,
                   rows = 2:(n_rows + 1),
                   cols = c(6, 9), gridExpand = TRUE)

# Column widths
openxlsx::setColWidths(wb, "Table 1",
                        cols = 1:n_cols,
                        widths = c(32, 22, 12, 12, 12, 10, 18, 16, 10))

# Row heights
openxlsx::setRowHeights(wb, "Table 1",
                         rows = 1, heights = 40)
openxlsx::setRowHeights(wb, "Table 1",
                         rows = 2:(n_rows + 1), heights = 18)

# Footnote
openxlsx::writeData(wb, "Table 1",
                    paste0("Values are median (Q1, Q3) for continuous variables",
                           "and n (%) for categorical variables. ",
                           "p-values from Mann-Whitney U test (continuous) ",
                           "and Fisher's exact test (categorical). ",
                           "APO = adverse pregnancy outcome; ",
                           "PPWR = postpartum weight retention; ",
                           "HDP = hypertensive disorders of pregnancy; ",
                           "GDM = gestational diabetes mellitus. ",
                           "\u2020Single GDM case; interpret with extreme caution."),
                    startRow = n_rows + 3, startCol = 1)

openxlsx::saveWorkbook(wb,
  file.path(out, "Table1.xlsx"),
  overwrite = TRUE)
message("  Saved: Table1.xlsx")

# ── Save Word document ────────────────────────────────────
message("Saving Table1.docx...")

# Use officer package for Word output
if (!requireNamespace("officer",   quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")

library(officer)
library(flextable)

# Build flextable
ft <- flextable::flextable(tbl) %>%

  # Font throughout
  flextable::font(fontname = "Helvetica", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%

  # Header styling
  flextable::bold(part = "header") %>%
  flextable::bg(bg = "#2C3E50", part = "header") %>%
  flextable::color(color = "#FFFFFF", part = "header") %>%
  flextable::align(align = "center", part = "header") %>%

  # Bold variable name rows
  flextable::bold(
    i   = which(tbl[, 2] == "" & tbl[, 1] != ""),
    j   = 1,
    part = "body"
  ) %>%

  # Indent level rows
  flextable::padding(
    i   = which(tbl[, 2] != ""),
    j   = 2,
    padding.left = 15,
    part = "body"
  ) %>%

  # Centre p-value columns
  flextable::align(j = c(6, 9), align = "center", part = "body") %>%

  # Alternate row shading for readability
  flextable::bg(
    i   = seq(2, nrow(tbl), by = 2),
    bg  = "#F7F9FB",
    part = "body"
  ) %>%

  # Borders
  flextable::border_outer(
    border = officer::fp_border(color = "#CCCCCC", width = 1)) %>%
  flextable::border_inner_h(
    border = officer::fp_border(color = "#EEEEEE", width = 0.5)) %>%
  flextable::border_inner_v(
    border = officer::fp_border(color = "#EEEEEE", width = 0.5)) %>%

  # Column widths (in inches — 170mm total = ~6.7 inches)
  flextable::width(j = 1, width = 2.2) %>%
  flextable::width(j = 2, width = 1.5) %>%
  flextable::width(j = c(3,4,5,7,8), width = 0.75) %>%
  flextable::width(j = c(6,9), width = 0.6) %>%

  # Wrap header text
  flextable::set_header_labels(values = as.list(setNames(
    colnames(tbl), colnames(tbl)))) %>%

  flextable::autofit(add_w = 0, add_h = 0)

# Footnote
ft <- flextable::add_footer_lines(ft,
  values = paste0(
    "Values are mean (SD) for continuous variables and n (%) for categorical variables. ",
    "p-values from Mann-Whitney U test (continuous variables) and ",
    "Fisher\u2019s exact test (categorical variables). ",
    "APO = adverse pregnancy outcome; PPWR = postpartum weight retention; ",
    "HDP = hypertensive disorders of pregnancy; GDM = gestational diabetes mellitus. ",
    "\u2020Single GDM case (n=1); interpret with extreme caution."
  )) %>%
  flextable::font(fontname = "Helvetica", part = "footer") %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::color(color = "#666666", part = "footer")

# Save to docx
doc <- officer::read_docx() %>%
  officer::body_add_par("Table 1. Sample characteristics",
                         style = "heading 2") %>%
  flextable::body_add_flextable(ft) %>%
  officer::body_add_par("", style = "Normal")

print(doc, target = file.path(out, "Table1.docx"))
message("  Saved: Table1.docx")

message("\nScript 15 complete -> ", out)
message("  Table1.xlsx — Excel version")
message("  Table1.docx — Word version for manuscript")