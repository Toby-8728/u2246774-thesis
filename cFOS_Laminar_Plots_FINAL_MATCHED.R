library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

raw <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
raw$MODEL <- trimws(raw$MODEL)
raw$AGENT <- factor(trimws(raw$AGENT),
  levels = c("AMG333", "TRPM2-KO", "AMG333 + TRPM2-KO"))

raw$L1_PERCENT <- ifelse(raw$CFOS_TOTAL > 0, raw$L1_CFOS / raw$CFOS_TOTAL * 100, NA)
raw$L2_PERCENT <- ifelse(raw$CFOS_TOTAL > 0, raw$L2_CFOS / raw$CFOS_TOTAL * 100, NA)
raw$L3_PERCENT <- ifelse(raw$CFOS_TOTAL > 0, raw$L34_CFOS / raw$CFOS_TOTAL * 100, NA)
agent_colours <- c("AMG333" = "#F0202D", "TRPM2-KO" = "#69A5F7",
  "AMG333 + TRPM2-KO" = "#6536D4")

make_cfos_plot <- function(columns, ymax, ybreak, y_label, percentage = FALSE) {
  d <- raw %>%
    select(MODEL, AGENT, all_of(columns)) %>%
    pivot_longer(all_of(columns), names_to = "LAMINA", values_to = "VALUE") %>%
    mutate(LAMINA_NUMBER = match(LAMINA, columns),
      X = LAMINA_NUMBER + if_else(MODEL == "Cold Pain", 0, 3) +
        c(-0.30, 0, 0.30)[as.integer(AGENT)])

  s <- d %>%
    group_by(MODEL, AGENT, LAMINA, X) %>%
    summarise(MEAN = mean(VALUE, na.rm = TRUE), SD = sd(VALUE, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(LOW = pmax(MEAN - SD, 0), HIGH = MEAN + SD)
  if (percentage) s$HIGH <- pmin(s$HIGH, 100)

  ggplot() +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#C0C0C0") +
    annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = "#969696") +
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "#C0C0C0") +
    annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = "#969696") +
    annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "#C0C0C0") +
    annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf, fill = "#969696") +
    geom_vline(xintercept = 3.5, linetype = "dashed", linewidth = 0.7, colour = "black") +
    geom_point(data = d, aes(X, VALUE, colour = AGENT), shape = 1, size = 1.9, stroke = 0.8,
      position = position_jitter(width = 0.035, height = 0, seed = 123)) +
    geom_errorbar(data = s, aes(X, ymin = LOW, ymax = HIGH),
      width = 0.10, linewidth = 0.8, colour = "black") +
    geom_point(data = s, aes(X, MEAN, fill = AGENT), shape = 21,
      size = 4, stroke = 0.7, colour = "black") +
    scale_x_continuous(breaks = 1:6, labels = rep(c("L1", "L2", "L3+"), 2),
      limits = c(0.5, 6.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, ymax, ybreak), expand = c(0, 0)) +
    scale_colour_manual(values = agent_colours, guide = "none") +
    scale_fill_manual(values = agent_colours, guide = "none") +
    annotate("text", x = 2, y = ymax * 1.03, label = "COLD PAIN", fontface = "bold", size = 4.3) +
    annotate("text", x = 5, y = ymax * 1.03, label = "COLD ALLODYNIA", fontface = "bold", size = 4.3) +
    labs(x = NULL, y = y_label) +
    theme_classic() +
    theme(axis.text = element_text(size = 11, colour = "black"),
      axis.text.x = element_text(face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 8)),
      axis.line = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks.x = element_blank(), legend.position = "none",
      plot.margin = margin(t = 22, r = 8, b = 8, l = 8)) +
    coord_cartesian(ylim = c(0, ymax), clip = "off")
}

p_count <- make_cfos_plot(c("L1_CFOS", "L2_CFOS", "L34_CFOS"), 40, 5, "cFOS+ cells")
p_percent <- make_cfos_plot(c("L1_PERCENT", "L2_PERCENT", "L3_PERCENT"),
  100, 20, "% cFOS+ cells", percentage = TRUE)

if (capabilities("aqua")) grDevices::quartz(width = 10, height = 6.5)
print(p_count)
if (capabilities("aqua")) grDevices::quartz(width = 10, height = 6.5)
print(p_percent)


