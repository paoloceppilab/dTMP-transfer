#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx)) {
    return(default)
  }
  if (idx == length(args)) {
    stop(sprintf("Missing value after %s", flag), call. = FALSE)
  }
  args[[idx + 1]]
}

script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 0) {
    return(NA_character_)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

script_file <- script_path()
archive_dir <- if (is.na(script_file)) {
  normalizePath(getwd(), mustWork = TRUE)
} else {
  dirname(dirname(script_file))
}

data_path <- arg_value(
  "--data",
  file.path(archive_dir, "data", "alteration_summary.csv")
)
output_path <- arg_value(
  "--output",
  file.path(archive_dir, "results", "extended_3c_pie.svg")
)

required_columns <- c(
  "category",
  "count",
  "expected_percent",
  "display_label",
  "color"
)
expected_categories <- c(
  "Amplification in TYMS and TK1",
  "Amplification in TYMS",
  "Amplification in TK1",
  "Missense in TYMS",
  "Missense in TK1",
  "None"
)
expected_total <- 496L

summary_table <- read.csv(
  data_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns <- setdiff(required_columns, names(summary_table))
if (length(missing_columns) > 0) {
  stop(
    sprintf("Input table is missing required columns: %s", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}

if (!identical(summary_table$category, expected_categories)) {
  stop("Input categories are not in the expected Extended Data Fig. 3C order.", call. = FALSE)
}

if (any(is.na(summary_table$count)) || any(summary_table$count < 0)) {
  stop("Counts must be non-negative numeric values with no missing entries.", call. = FALSE)
}

if (!all(summary_table$count == as.integer(summary_table$count))) {
  stop("Counts must be integers.", call. = FALSE)
}

total_count <- sum(summary_table$count)
if (total_count != expected_total) {
  stop(
    sprintf("Expected total count %d, observed %d.", expected_total, total_count),
    call. = FALSE
  )
}

computed_percent <- round(summary_table$count / total_count * 100, 2)
if (any(abs(computed_percent - summary_table$expected_percent) > 0.005)) {
  details <- paste(
    sprintf(
      "%s: expected %.2f, computed %.2f",
      summary_table$category,
      summary_table$expected_percent,
      computed_percent
    ),
    collapse = "; "
  )
  stop(sprintf("Percentages do not match counts: %s", details), call. = FALSE)
}

expected_labels <- sprintf(
  "%s (%.2f%%)",
  summary_table$category,
  summary_table$expected_percent
)
if (!identical(summary_table$display_label, expected_labels)) {
  stop("Display labels do not match category names and expected percentages.", call. = FALSE)
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

escape_xml <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

point_on_circle <- function(cx, cy, radius, angle_radians) {
  c(
    x = cx + radius * cos(angle_radians),
    y = cy + radius * sin(angle_radians)
  )
}

width <- 800
height <- 480
cx <- 220
cy <- 240
radius <- 165
angle <- -pi / 2
fractions <- summary_table$count / total_count

paths <- character(nrow(summary_table))
for (i in seq_len(nrow(summary_table))) {
  next_angle <- angle + 2 * pi * fractions[[i]]
  start <- point_on_circle(cx, cy, radius, angle)
  end <- point_on_circle(cx, cy, radius, next_angle)
  large_arc <- ifelse((next_angle - angle) > pi, 1, 0)
  paths[[i]] <- sprintf(
    paste0(
      "<path d=\"M %.3f %.3f L %.3f %.3f ",
      "A %.3f %.3f 0 %d 1 %.3f %.3f Z\" ",
      "fill=\"%s\" stroke=\"#FFFFFF\" stroke-width=\"2\"/>"
    ),
    cx,
    cy,
    start[["x"]],
    start[["y"]],
    radius,
    radius,
    large_arc,
    end[["x"]],
    end[["y"]],
    summary_table$color[[i]]
  )
  angle <- next_angle
}

legend_x <- 440
legend_y <- 130
legend_gap <- 38
legend_lines <- unlist(lapply(seq_len(nrow(summary_table)), function(i) {
  y <- legend_y + (i - 1) * legend_gap
  c(
    sprintf(
      "<rect x=\"%d\" y=\"%d\" width=\"18\" height=\"18\" fill=\"%s\" stroke=\"none\"/>",
      legend_x,
      y - 14,
      summary_table$color[[i]]
    ),
    sprintf(
      "<text x=\"%d\" y=\"%d\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"16\" fill=\"#222222\">%s</text>",
      legend_x + 30,
      y,
      escape_xml(summary_table$display_label[[i]])
    )
  )
}))

svg_lines <- c(
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
  sprintf(
    "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\" role=\"img\" aria-label=\"TCGA-LUAD TYMS and TK1 alteration status pie chart\">",
    width,
    height,
    width,
    height
  ),
  "<rect width=\"100%\" height=\"100%\" fill=\"#FFFFFF\"/>",
  "<g id=\"pie\">",
  paths,
  "</g>",
  "<g id=\"legend\">",
  legend_lines,
  "</g>",
  "</svg>"
)

writeLines(svg_lines, con = output_path, useBytes = TRUE)

message(sprintf("Wrote %s", normalizePath(output_path, mustWork = TRUE)))
