library(dplyr)
library(ggplot2)

# ================================================================
# load + merge datasets
# ================================================================

GHoct <- read.csv("greenhouse_octexp_analysis.csv") %>%
  mutate(
    dataset = "GH",
    week = as.integer(week),
    leaf = as.character(leaf),
    Predicted_N = as.numeric(Predicted_N),
    M_distance  = as.numeric(M_distance)
  ) %>%
  transmute(
    dataset,
    date,
    week,
    group,
    stage,
    pot,
    plant,
    plot = NA,              # GH doesn't have plot
    leaf,
    Predicted_N,
    M_distance
  )

B5oct <- read.csv("B5_octexp_analysis.csv") %>%
  mutate(
    dataset = "B5",
    week = as.integer(week),
    leaf = as.character(leaf),
    Predicted_N = as.numeric(Predicted_N),
    M_distance  = as.numeric(M_distance)
  ) %>%
  transmute(
    dataset,
    date,
    week,
    group,
    stage,
    pot = NA,               # B5 doesn't have pot/plant
    plant = NA,
    plot,
    leaf,
    Predicted_N,
    M_distance
  )

oct_all <- bind_rows(GHoct, B5oct)

# ================================================================
# average repeated measurements (within unit × leaf × week)
# ================================================================

oct_all_avg <- oct_all %>%
  mutate(
    unit_id = ifelse(dataset == "GH",
                     paste(pot, plant, sep = "_"),
                     as.character(plot))
  ) %>%
  group_by(dataset, date, week, group, stage, unit_id, leaf) %>%
  summarise(
    Predicted_N = ifelse(
      all(M_distance > 4),
      NA_real_,
      mean(Predicted_N[M_distance <= 4])
    ),
    .groups = "drop"
  )

# ================================================================
# average across units (week × leaf)
# ================================================================

oct_week_leaf <- oct_all_avg %>%
  group_by(dataset, week, leaf) %>%
  summarise(
    Predicted_N = ifelse(
      all(is.na(Predicted_N)),
      NA_real_,
      mean(Predicted_N, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  filter(
    (dataset == "GH" & leaf %in% c("l14", "l9")) |
      (dataset == "B5" & leaf %in% c("e1", "e3"))
  ) %>%
  arrange(dataset, leaf, week)

print(oct_week_leaf)

p_leaf <- ggplot(oct_week_leaf, aes(x = week, y = Predicted_N, color = leaf)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~dataset, scales = "free_y") +
  labs(x = "Weeks after flowering", y = "Average Predicted N", color = "Leaf") +
  theme_classic()

p_leaf

ggsave(
  filename = "output/OCTEXP_by_leaf.png",
  plot = p_leaf,
  width = 7,
  height = 5,
  dpi = 300
)

# ================================================================
# average across units (week × group, leaves combined)
# ================================================================

oct_week_group <- oct_all_avg %>%
  filter(
    (dataset == "GH" & leaf %in% c("l14", "l9")) |
      (dataset == "B5" & leaf %in% c("e1", "e3"))
  ) %>%
  group_by(dataset, week, group, unit_id) %>%
  summarise(
    Predicted_N = ifelse(
      all(is.na(Predicted_N)),
      NA_real_,
      mean(Predicted_N, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  group_by(dataset, week, group) %>%
  summarise(
    Predicted_N = ifelse(
      all(is.na(Predicted_N)),
      NA_real_,
      mean(Predicted_N, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  arrange(dataset, group, week)

print(oct_week_group)

p_group <- ggplot(
  oct_week_group,
  aes(x = week, y = Predicted_N, color = group, group = group)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~dataset, scales = "free_y") +
  scale_color_manual(
    values = c(
      hihi = "green3",
      lolo = "gold2",
      lohi = "purple3",
      hilo = "dodgerblue3",
      check = "grey60"
    )
  ) +
  labs(
    x = "Weeks after flowering",
    y = "Average Predicted N (leaves combined)",
    color = "Group"
  ) +
  theme_classic()

p_group

ggsave(
  filename = "output/OCTEXP_by_group.png",
  plot = p_group,
  width = 7,
  height = 5,
  dpi = 300
)

# ================================================================
# unit-level values (leaves combined)
# ================================================================

oct_unit_week <- oct_all_avg %>%
  filter(
    (dataset == "GH" & leaf %in% c("l14", "l9")) |
      (dataset == "B5" & leaf %in% c("e1", "e3"))
  ) %>%
  group_by(dataset, week, group, unit_id) %>%
  summarise(
    Predicted_N = ifelse(
      all(is.na(Predicted_N)),
      NA_real_,
      mean(Predicted_N, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(week)) %>%          # GH week 0 set to NA in the CSV; drop those rows
  filter(week >= 0 & week <= 4) %>% # enforce shared week range for GH and B5
  arrange(dataset, group, unit_id, week)

print(oct_unit_week)

# ================================================================
# deviation from group mean (unit − group, per week)
# ================================================================

oct_unit_deviation <- oct_unit_week %>%
  left_join(
    oct_week_group,
    by = c("dataset", "week", "group"),
    suffix = c("_unit", "_group")
  ) %>%
  mutate(
    deviation = Predicted_N_unit - Predicted_N_group
  )

print(oct_unit_deviation)

# ================================================================
# spaghetti plot (colored unit trajectories + emphasized group means)
# ================================================================

p_group_spaghetti_colored <- ggplot() +
  geom_line(
    data = oct_unit_week,
    aes(
      x = week,
      y = Predicted_N,
      group = unit_id,
      color = group
    ),
    linewidth = 0.4,
    alpha = 0.25
  ) +
  geom_line(
    data = oct_week_group,
    aes(
      x = week,
      y = Predicted_N,
      color = group,
      group = group
    ),
    linewidth = 1.6
  ) +
  geom_point(
    data = oct_week_group,
    aes(
      x = week,
      y = Predicted_N,
      color = group
    ),
    size = 3
  ) +
  facet_wrap(~dataset, scales = "free_y") +
  scale_color_manual(
    values = c(
      hihi  = "green3",
      lolo  = "gold2",
      lohi  = "purple3",
      hilo  = "dodgerblue3",
      check = "grey60"
    )
  ) +
  labs(
    x = "Weeks after flowering",
    y = "Predicted N (leaves combined)",
    color = "Group"
  ) +
  scale_y_continuous(limits = c(2.7, 4.7)) +
  theme_classic()

p_group_spaghetti_colored

ggsave(
  filename = "output/OCTEXP_group_spaghetti_colored.png",
  plot = p_group_spaghetti_colored,
  width = 8,
  height = 5,
  dpi = 300
)

ok# ================================================================
# b5 spaghetti by leaf (e1 left, e3 right)
# ================================================================

# unit-level values per leaf (each plot contributes one value per week per leaf)
b5_unit_week_leaf <- oct_all_avg %>%
  filter(dataset == "B5", leaf %in% c("e1", "e3")) %>%
  group_by(week, group, unit_id, leaf) %>%
  summarise(
    Predicted_N = ifelse(
      all(is.na(Predicted_N)),
      NA_real_,
      mean(Predicted_N, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(week)) %>%
  arrange(leaf, group, unit_id, week)

print(b5_unit_week_leaf)

# group means per leaf
b5_week_group_leaf <- b5_unit_week_leaf %>%
  group_by(week, group, leaf) %>%
  summarise(
    Predicted_N = ifelse(
      all(is.na(Predicted_N)),
      NA_real_,
      mean(Predicted_N, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  arrange(leaf, group, week)

print(b5_week_group_leaf)

p_b5_spaghetti_by_leaf <- ggplot() +
  # thin unit trajectories
  geom_line(
    data = b5_unit_week_leaf,
    aes(x = week, y = Predicted_N, group = unit_id, color = group),
    linewidth = 0.4,
    alpha = 0.25
  ) +
  # thick group means
  geom_line(
    data = b5_week_group_leaf,
    aes(x = week, y = Predicted_N, group = group, color = group),
    linewidth = 1.6
  ) +
  # big points on group means
  geom_point(
    data = b5_week_group_leaf,
    aes(x = week, y = Predicted_N, color = group),
    size = 3
  ) +
  facet_grid(
    . ~ leaf,
    labeller = as_labeller(c(
      e1 = "e1 (younger leaf)",
      e3 = "e3 (older leaf)"
    ))
  ) +
  scale_color_manual(
    values = c(
      hihi  = "green3",
      lolo  = "gold2",
      lohi  = "purple3",
      hilo  = "dodgerblue3",
      check = "grey60"
    )
  ) +
  labs(
    x = "Weeks after flowering",
    y = "Predicted N",
    color = "Group"
  ) +
  theme_classic()

p_b5_spaghetti_by_leaf

ggsave(
  filename = "output/B5_spaghetti_by_leaf.png",
  plot = p_b5_spaghetti_by_leaf,
  width = 8,
  height = 5,
  dpi = 300
)

# ================================================================
# b5 auc by leaf (e1 left, e3 right)
# ================================================================

# start from the unit-level weekly values you already made
# (b5_unit_week_leaf has: week, group, unit_id, leaf, Predicted_N)

b5_auc <- b5_unit_week_leaf %>%
  arrange(unit_id, leaf, week) %>%
  group_by(group, unit_id, leaf) %>%
  summarise(
    AUC = {
      w <- week
      y <- Predicted_N
      
      ok <- !is.na(w) & !is.na(y)
      w <- w[ok]
      y <- y[ok]
      
      if (length(y) < 2) {
        NA_real_
      } else {
        sum(((y[-1] + y[-length(y)]) / 2) * (w[-1] - w[-length(w)]))
      }
    },
    .groups = "drop"
  )

print(b5_auc)

b5_auc_means <- b5_auc %>%
  group_by(group, leaf) %>%
  summarise(
    mean_auc = mean(AUC, na.rm = TRUE),
    sd_auc   = sd(AUC, na.rm = TRUE),
    .groups  = "drop"
  )

print(b5_auc_means)

p_b5_auc_by_leaf <- ggplot(b5_auc, aes(x = group, y = AUC, fill = group)) +
  geom_jitter(
    shape = 21, size = 4, alpha = 0.8,
    position = position_jitter(width = 0.1),
    color = "grey30"
  ) +
  # mean +/- sd (error bar)
  geom_errorbar(
    data = b5_auc_means,
    aes(x = group, ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc),
    width = 0, linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  # mean point (black)
  geom_point(
    data = b5_auc_means,
    aes(x = group, y = mean_auc),
    shape = 16, size = 3,
    inherit.aes = FALSE,
    color = "black"
  ) +
  facet_grid(
    . ~ leaf,
    labeller = as_labeller(c(
      e1 = "e1 (younger leaf)",
      e3 = "e3 (older leaf)"
    ))
  ) +
  scale_fill_manual(
    values = c(
      hihi  = "green3",
      lolo  = "gold2",
      lohi  = "purple3",
      hilo  = "dodgerblue3",
      check = "grey60"
    )
  ) +
  labs(
    x = "Group",
    y = "AUC of predicted N (weeks after flowering)"
  ) +
  theme_classic() +
  guides(fill = "none")

p_b5_auc_by_leaf

ggsave(
  filename = "output/B5_AUC_by_leaf.png",
  plot = p_b5_auc_by_leaf,
  width = 8,
  height = 5,
  dpi = 300
)