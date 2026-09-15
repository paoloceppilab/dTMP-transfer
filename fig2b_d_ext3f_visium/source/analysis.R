options(future.globals.maxSize = 32000 * 1024^2)
library(Seurat)
library(tidyverse)
library(dbscan)
library(magrittr)
library(multcompView) # For significance markers
set.seed(123)

# Load data from the parent directory
lung <- Load10X_Spatial(
    data.dir = "..",
    filename = "CytAssist_FFPE_Human_Lung_Squamous_Cell_Carcinoma_filtered_feature_bc_matrix.h5"
)

genes_of_interest <- c("TYMS", "TK1")
raw_counts_filtered <- FetchData(lung, vars = genes_of_interest, layer = "counts")
raw_counts_filtered$TYMS_level <- ifelse(raw_counts_filtered$TYMS > 0, "Pos", "Neg")
raw_counts_filtered$TK1_level <- ifelse(raw_counts_filtered$TK1 > 0, "Pos", "Neg")

lung <- AddMetaData(lung, metadata = raw_counts_filtered[, c("TYMS_level", "TK1_level")])

tyms_pos_tk1_pos_filtered <- sum(raw_counts_filtered$TYMS > 0 & raw_counts_filtered$TK1 > 0)
tyms_pos_tk1_zero_filtered <- sum(raw_counts_filtered$TYMS > 0 & raw_counts_filtered$TK1 == 0)
tyms_zero_tk1_zero_filtered <- sum(raw_counts_filtered$TYMS == 0 & raw_counts_filtered$TK1 == 0)
tyms_zero_tk1_pos_filtered <- sum(raw_counts_filtered$TYMS == 0 & raw_counts_filtered$TK1 > 0)

coords <- Seurat::GetTissueCoordinates(lung, scale = "lowres") |>
    rownames_to_column(var = "cellid") |>
    as_tibble() |>
    rename(all_of(c(row_scaled = "x", col_scaled = "y")))

meta <- rownames_to_column(lung@meta.data, var = "cellid") %>% as_tibble()
coord <- list(meta, coords) |> purrr::reduce(inner_join, by = "cellid")

coord <- coord |> dplyr::mutate(pos_status = case_when(
    (TYMS_level == "Pos" & TK1_level == "Neg") ~ "TYMS_Pos",
    (TYMS_level == "Neg" & TK1_level == "Pos") ~ "TK1_Pos",
    (TYMS_level == "Pos" & TK1_level == "Pos") ~ "Double_Pos",
    (TYMS_level == "Neg" & TK1_level == "Neg") ~ "Double_Neg"
))

coord$pos_status <- as.factor(coord$pos_status)
coordmatscl <- coord |>
    dplyr::select(c("cellid", "row_scaled", "col_scaled")) |>
    column_to_rownames(var = "cellid") |>
    as.matrix()

# Neighborhood Analysis Functions
find_mode <- function(x) {
    u <- unique(x)
    tab <- tabulate(match(x, u))
    u[tab == max(tab)]
}

nn_df <- function(nn, coord, colname) {
    nnl <- lengths(nn$id) |>
        enframe() |>
        rename(all_of(c(
            cellid = "name", numneigh = "value"
        )))

    complete_neighbors <- find_mode(nnl$numneigh)

    nbdmat <- purrr::map(nn$id, ~ coord[[colname]][.x] %>% table())
    nn_matdf <- bind_rows(nbdmat, .id = "cellid") |>
        as.data.frame() %>%
        replace(is.na(.), 0)

    coord_nn <- list(coord, nn_matdf, nnl) |> purrr::reduce(inner_join, by = "cellid")

    col_vals <- unique(coord[[colname]]) |> as.character()

    outdata <- coord_nn |>
        dplyr::filter(numneigh == complete_neighbors) |>
        group_by(pos_status) |>
        summarise(across(all_of(col_vals), list(mean = mean, sd = sd), .names = "{.col}_{.fn}")) |>
        pivot_longer(
            cols = matches("_mean$"),
            names_to = "neighbor_status",
            values_to = "meanvalue",
            names_pattern = "(.*)_mean"
        ) |>
        pivot_longer(
            cols = matches("_sd$"),
            names_to = "neighbor_status_sd",
            values_to = "sd_value",
            names_pattern = "(.*)_sd"
        ) |>
        filter(neighbor_status == neighbor_status_sd) |>
        select(-neighbor_status_sd)

    outdata$numneigh <- complete_neighbors
    outdata$radius <- (nn$eps - 2) / 6
    outdata$perc <- (outdata$meanvalue / complete_neighbors) * 100
    outdata$sd_perc <- (outdata$sd_value / complete_neighbors) * 100

    return(outdata)
}

# Run Neighborhood Analysis
pos_status_list <- list()
for (ix in 1:10) {
    eps <- (ix * 6) + 2
    nn <- dbscan::frNN(x = coordmatscl, eps = eps)
    pos_status_list[[ix]] <- nn_df(nn, coord, "pos_status")
}
pos_status_nn_data <- bind_rows(pos_status_list)

# ANOVA and Significance
anova_results <- pos_status_nn_data |>
    group_by(radius) |>
    summarise(
        anova_pval = list(aov(perc ~ neighbor_status, data = cur_data())) |>
            purrr::map_dbl(~ summary(.x)[[1]][["Pr(>F)"]][1])
    ) |>
    mutate(
        significance = case_when(
            anova_pval < 0.0001 ~ "****",
            anova_pval < 0.001 ~ "***",
            anova_pval < 0.01 ~ "**",
            anova_pval < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    )

pos_status_nn_data <- left_join(pos_status_nn_data, anova_results, by = "radius")

# 1. Neighborhood Plot
pos_status_nn_plot <- pos_status_nn_data |>
    ggplot(aes(radius, perc, color = neighbor_status)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = perc - sd_perc, ymax = perc + sd_perc), width = 0.2) +
    facet_wrap(~pos_status) +
    scale_x_continuous(breaks = c(seq(0, 10, 2)), limits = c(0, 11)) +
    theme_bw() +
    labs(x = "radius", y = "percentage of cells", title = "Neighborhood Analysis by TYMS/TK1 Status") +
    scale_color_manual(values = c("Double_Neg" = "red", "Double_Pos" = "darkgray", "TK1_Pos" = "purple", "TYMS_Pos" = "orange")) +
    geom_text(data = anova_results, aes(x = radius, y = 100, label = significance), size = 5, color = "black")

ggsave("neighborhood_plot.svg", plot = pos_status_nn_plot, width = 12, height = 10)
ggsave("neighborhood_plot.png", plot = pos_status_nn_plot, width = 12, height = 10, dpi = 300)

# 2. Spatial Plot
lung$pos_status <- coord$pos_status
lung$pos_status <- factor(lung$pos_status, levels = c("Double_Neg", "Double_Pos", "TK1_Pos", "TYMS_Pos"))

spatial_plot <- SpatialPlot(
    lung,
    group.by = "pos_status",
    pt.size.factor = 3,
    images = NULL,
    stroke = 0,
    cols = c(
        "Double_Neg" = "red",
        "Double_Pos" = "darkgray",
        "TK1_Pos" = "purple",
        "TYMS_Pos" = "orange"
    )
) +
    theme_void() +
    labs(title = "Spatial Distribution of TYMS/TK1 Status")

ggsave("spatial_plot.svg", plot = spatial_plot, width = 10, height = 8)
ggsave("spatial_plot.png", plot = spatial_plot, width = 10, height = 8, dpi = 300)

# 3. Count Table
results_table <- data.frame(
    Category = c(
        "TYMS_pos_TK1_pos",
        "TYMS_pos_TK1_neg",
        "TYMS_neg_TK1_neg",
        "TYMS_neg_TK1_pos"
    ),
    Counts = c(
        tyms_pos_tk1_pos_filtered,
        tyms_pos_tk1_zero_filtered,
        tyms_zero_tk1_zero_filtered,
        tyms_zero_tk1_pos_filtered
    )
)

write.csv(results_table, "tyms_tk1_counts.csv", row.names = FALSE)

# 4. Neighborhood Analysis Data Points
write.csv(pos_status_nn_data, "neighborhood_data_points.csv", row.names = FALSE)

print("Analysis complete. Results saved in current directory.")
