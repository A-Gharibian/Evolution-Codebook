# NOTE: This script depends on `rho_filtered` and `har_filtered`
# produced in script 2.

# -------------------------------------------------------------------
# RHO: Stacked bar chart — top 10 organisms by SRA hit share
# -------------------------------------------------------------------

top_10_rho_organisms <- rho_filtered |>
  count(organism) |>
  slice_max(n, n = 10) |>
  pull(organism)

rho_plot_data <- rho_filtered |>
  filter(organism %in% top_10_rho_organisms)

ggplot(rho_plot_data, aes(x = Query_Species, fill = organism)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  labs(
    title = "NCBI Metadata vs. Actual RHO Sequence Matches",
    x     = "RHO Sequence Query",
    y     = "Percentage of SRA Hits",
    fill  = "NCBI Organism Label"
  ) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right"
  )

# -------------------------------------------------------------------
# HAR1A: Stacked bar chart — top 16 organisms by SRA hit share
# -------------------------------------------------------------------

top_16_har_organisms <- har_filtered |>
  count(organism, sort = TRUE) |>
  slice_max(n, n = 16) |>
  pull(organism)

har_plot_data <- har_filtered |>
  filter(organism %in% top_16_har_organisms) |>
  mutate(Query_Species = fct_reorder(
    Query_Species,
    as.numeric(gsub("Variation", "", Query_Species))
  ))

ggplot(har_plot_data, aes(x = Query_Species, fill = organism)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  labs(
    title = "NCBI Metadata vs. Actual HAR1A Sequence Matches",
    x     = "HAR1A Variation Query",
    y     = "Percentage of SRA Hits",
    fill  = "NCBI Organism Label"
  ) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right"
  )

