library(dplyr)
library(ggplot2)
library(ggnewscale)
library(viridisLite)

## 1. define your orders
method_order  <- c("Decipher","NicheNet","LIANA+","NATMI","Connectome")
dataset_order <- c(
  "5yr_pic","bcg","cord_pic","covid",
  "erp","lupus","sepsis","tnbc",
  "cz_influenza","cz_hpap_t1d_islets","cz_hnscc_hpv",
  "cz_human_kidney_v1.5","cz_cf_bronchial_biopsy",
  "SevCOVID_Azimuthl2","MilCOVID_Azimuthl2"
)

## 2. build the “full grid” of every (Method1,Method2) × (tile_row,tile_col)
grid_df <- expand.grid(
  Method1  = method_order,
  Method2  = method_order,
  tile_row = 1:4,
  tile_col = 1:4,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Method1 = factor(Method1, levels = method_order),
    Method2 = factor(Method2, levels = method_order),
    # figure out which dataset should live here (1–15 → real names; 16 → NA)
    dataset_index = (tile_row - 1)*4 + tile_col,
    Dataset       = ifelse(dataset_index <= length(dataset_order),
                           dataset_order[dataset_index],
                           NA_character_)
  )

## 3. left-join your real values on to that grid
heat_df <- grid_df %>%
  left_join(combined_df, by = c("Method1","Method2","Dataset")) %>%
  mutate(
    big_row = match(Method1, method_order) - 1,
    big_col = match(Method2, method_order) - 1,
    x       = big_col*4 + tile_col,
    y       = (4 - tile_row) + big_row*4
  )

# 5. build a small data.frame of the 5×5 block‐centres
border_df <- expand.grid(
  Method1 = method_order,
  Method2 = method_order,
  stringsAsFactors = FALSE
) %>% 
  mutate(
    big_row = match(Method1, method_order) - 1,
    big_col = match(Method2, method_order) - 1,
    # centre of each 4×4 block:
    x = big_col*4 + 2.5,
    y = big_row*4 + 1.5
  )

## 4. plot—with one fill scale for Spearman (upper triangle)
##    and a second for k_value (lower triangle)
p <- ggplot() +
  # Diagonal override
  geom_tile(
    data = heat_df %>% filter(big_row == big_col),
    aes(x = x, y = y),
    fill = "lightgray",
    inherit.aes = FALSE
  ) +
  # upper tri: Spearman
  geom_tile(
    data = heat_df %>% filter(big_row < big_col),
    aes(x, y, fill = Spearman)
  ) +
  scale_fill_gradient2(
    name    = "Spearman",
    low     = "#b2182b", mid = "white", high = "#008837",
    limits  = c(-1,1),
    na.value= "grey80"
  ) +

  # start a fresh fill mapping
  new_scale_fill() +

  # lower tri: k_value
  geom_tile(
    data = heat_df %>% filter(big_row > big_col),
    aes(x, y, fill = k_value)
  ) +
  scale_fill_viridis_c(
    name    = "Search-space\n(k value)",
    option  = "B", end = 0.95,
    limits  = c(100, max(heat_df$k_value, na.rm = TRUE)),
    na.value= "grey80"
  ) +
  # --- border layer ---
  geom_tile(
    data      = border_df,
    aes(x, y),
    width     = 4,        # span 4 tiles
    height    = 4,
    fill      = NA,       # transparent
    color     = "white",  # white outline
    size      = 0.8,      # line thickness
    inherit.aes = FALSE
  ) +
  # tidy up axes so each 4×4 block is labelled by method
  coord_fixed() +
  scale_x_continuous(
    expand = c(0,0),
    breaks = (0:4)*4 + 2.5,
    labels = method_order,
    position = "bottom"
  ) +
  scale_y_reverse(
    expand = c(0,0),
    breaks = (0:4)*4 + 1.5,
    labels = method_order
  ) +

  theme_minimal(base_size = 12) +
  theme(
    axis.title   = element_blank(),
    axis.text.x  = element_text(
      size   = 17,
      face   = "bold",
      angle  = 45,
      hjust  = 1              # right-justify along the diagonal
    ),
    axis.text.y  = element_text(face = "bold", size = 17,vjust=1),
    panel.grid   = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(figures_folder,"big_heatmap_v3.png"),p)