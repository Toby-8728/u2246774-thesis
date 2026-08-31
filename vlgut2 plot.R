library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
.
data <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
data <- data %>%
  mutate(MODEL = factor(trimws(MODEL), levels = c("Cold Pain", "Cold Allodynia")),
    AGENT = factor(trimws(AGENT), levels = c("AMG333", "TRPM2-KO", "AMG333 + TRPM2-KO")),
    MARKER = trimws(MARKER),
    L1_PERCENT = ifelse(L1_CFOS > 0, L1_TOTAL / L1_CFOS * 100, NA),
    L2_PERCENT = ifelse(L2_CFOS > 0, L2_TOTAL / L2_CFOS * 100, NA),
    L3_PERCENT = ifelse(L34_CFOS > 0, L34_TOTAL / L34_CFOS * 100, NA))

agent_colours <- c("AMG333" = "#F0202D", "TRPM2-KO" = "#69A5F7",
  "AMG333 + TRPM2-KO" = "#6536D4")

make_total_plot <- function(marker, ymax, ybreak) {
  d <- subset(data, MARKER == marker)
  s <- d %>%
    group_by(MODEL, AGENT) %>%
    summarise(MEAN = mean(SLICE_TOTAL, na.rm = TRUE),
      SD = sd(SLICE_TOTAL, na.rm = TRUE), .groups = "drop")

  ggplot(s, aes(MODEL, MEAN, fill = AGENT)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.68,
      colour = "black", linewidth = 0.5) +
    geom_errorbar(aes(ymin = pmax(MEAN - SD, 0), ymax = MEAN + SD),
      position = position_dodge(width = 0.8), width = 0.15, linewidth = 0.8) +
    geom_point(data = d, aes(MODEL, SLICE_TOTAL, colour = AGENT), inherit.aes = FALSE,
      shape = 1, size = 2, stroke = 0.8,
      position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.8, seed = 123)) +
    scale_fill_manual(values = agent_colours, name = NULL) +
    scale_colour_manual(values = agent_colours, guide = "none") +
    scale_y_continuous(breaks = seq(0, ymax, ybreak), limits = c(0, ymax),
      expand = expansion(mult = c(0, 0.03))) +
    labs(x = NULL, y = paste0(marker, " count")) +
    theme_classic() +
    theme(axis.text = element_text(size = 11, colour = "black"),
      axis.text.x = element_text(face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "top", legend.text = element_text(size = 7))
}

make_lamina_plot <- function(marker, type, ymax, ybreak = NULL) {
  d <- subset(data, MARKER == marker)
  if (type == "count") {
    columns <- c("L1_TOTAL", "L2_TOTAL", "L34_TOTAL")
    y_label <- paste0(marker, " count")
    lower_limit <- -0.5
    panel_top <- ymax
    breaks <- if (is.null(ybreak)) pretty(c(0, ymax)) else seq(0, ymax, ybreak)
  } else {
    columns <- c("L1_PERCENT", "L2_PERCENT", "L3_PERCENT")
    y_label <- paste0("% ", marker, "+ cells")
    lower_limit <- -3
    panel_top <- 103
    breaks <- seq(0, 100, 20)
  }

  d <- d %>%
    select(MODEL, AGENT, all_of(columns)) %>%
    pivot_longer(all_of(columns), names_to = "LAMINA", values_to = "VALUE") %>%
    mutate(LAMINA_NUMBER = match(LAMINA, columns),
      X = LAMINA_NUMBER + if_else(MODEL == "Cold Pain", 0, 3) +
        c(-0.30, 0, 0.30)[as.integer(AGENT)])

  s <- d %>%
    filter(!is.na(VALUE)) %>%
    group_by(MODEL, AGENT, LAMINA, X) %>%
    summarise(MEAN = mean(VALUE), SD = sd(VALUE), .groups = "drop") %>%
    mutate(LOW = pmax(MEAN - SD, 0), HIGH = MEAN + SD)
  if (type == "percent") s$HIGH <- pmin(s$HIGH, 100)
  title_y <- panel_top + (panel_top - lower_limit) * 0.037

  ggplot() +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#C0C0C0") +
    annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = "#969696") +
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "#C0C0C0") +
    annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = "#969696") +
    annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "#C0C0C0") +
    annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "#969696") +
    geom_vline(xintercept = 3.5, linetype = "dashed", linewidth = 0.7) +
    geom_point(data = d, aes(X, VALUE, colour = AGENT), shape = 1,
      size = 1.9, stroke = 0.8, na.rm = TRUE,
      position = position_jitter(width = 0.025, height = 0, seed = 123)) +
    geom_errorbar(data = s, aes(X, ymin = LOW, ymax = HIGH),
      width = 0.10, linewidth = 0.8, colour = "#333333") +
    geom_point(data = s, aes(X, MEAN, fill = AGENT), shape = 21,
      size = 4, stroke = 0.7, colour = "black") +
    scale_x_continuous(breaks = 1:6, labels = rep(c("L1", "L2", "L3+"), 2),
      limits = c(0.5, 6.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = breaks, expand = c(0, 0)) +
    scale_colour_manual(values = agent_colours, guide = "none") +
    scale_fill_manual(values = agent_colours, name = NULL) +
    annotate("text", x = 2, y = title_y, label = "COLD PAIN", fontface = "bold", size = 3.7) +
    annotate("text", x = 5, y = title_y, label = "COLD ALLODYNIA", fontface = "bold", size = 3.7) +
    labs(x = NULL, y = y_label) +
    theme_classic() +
    theme(axis.text = element_text(size = 11, colour = "black"),
      axis.text.x = element_text(face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.ticks.x = element_blank(), legend.position = c(0.5, 1.18),
      legend.justification = c(0.5, 0.5), legend.direction = "vertical",
      legend.text = element_text(size = 7), legend.key.height = grid::unit(0.28, "cm"),
      legend.background = element_rect(fill = "white", colour = "black"),
      plot.margin = margin(t = 70, r = 8, b = 8, l = 8)) +
    coord_cartesian(ylim = c(lower_limit, panel_top), clip = "off")
}

VGLUT2_total <- make_total_plot("VGLUT2", 40, 5)
PAX2_total <- make_total_plot("PAX2", 10, 2)
VGLUT2_lamina_count <- make_lamina_plot("VGLUT2", "count", 17)
VGLUT2_lamina_percent <- make_lamina_plot("VGLUT2", "percent", 100, 20)
PAX2_lamina_count <- make_lamina_plot("PAX2", "count", 6, 1)
PAX2_lamina_percent <- make_lamina_plot("PAX2", "percent", 100, 20)

#fig output
plots <- list(VGLUT2_total = VGLUT2_total, PAX2_total = PAX2_total,
  VGLUT2_lamina_count = VGLUT2_lamina_count, VGLUT2_lamina_percent = VGLUT2_lamina_percent,
  PAX2_lamina_count = PAX2_lamina_count, PAX2_lamina_percent = PAX2_lamina_percent)
plot_width <- c(7, 7, 7, 8, 8, 8)
plot_height <- c(4.7, 4.7, 4.7, 6, 6, 6)
}
