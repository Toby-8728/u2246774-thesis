library(readxl)
library(dplyr)
library(ggplot2)

cfos <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
cfos <- cfos %>%
  mutate(MODEL = trimws(MODEL), AGENT = trimws(AGENT),
    AGENT = factor(AGENT, levels = c("AMG333", "TRPM2-KO", "AMG333 + TRPM2-KO")),
    MODEL_X = if_else(MODEL == "Cold Pain", 1, 2.2),
    AGENT_SHIFT = c(-0.28, 0, 0.28)[as.integer(AGENT)],
    X = MODEL_X + AGENT_SHIFT)

cfos_summary <- cfos %>%
  group_by(MODEL, AGENT, X) %>%
  summarise(MEAN = mean(CFOS_TOTAL), SD = sd(CFOS_TOTAL), .groups = "drop")

agent_colours <- c("AMG333" = "#F0202D", "TRPM2-KO" = "#69A5F7",
  "AMG333 + TRPM2-KO" = "#6536D4")

CFOS_total <- ggplot() +
  geom_col(data = cfos_summary, aes(X, MEAN, fill = AGENT),
    width = 0.22, colour = "black", linewidth = 0.7) +
  geom_errorbar(data = cfos_summary,
    aes(X, ymin = pmax(MEAN - SD, 0), ymax = MEAN + SD),
    width = 0.08, linewidth = 0.8, colour = "#333333") +
  geom_point(data = cfos, aes(X, CFOS_TOTAL, colour = AGENT),
    shape = 1, size = 1.9, stroke = 0.8,
    position = position_jitter(width = 0.035, height = 0, seed = 123)) +
  scale_fill_manual(values = agent_colours) +
  scale_colour_manual(values = agent_colours) +
  scale_x_continuous(breaks = NULL, limits = c(0.4, 2.8), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 70, 10), expand = c(0, 0)) +
  annotate("text", x = 1, y = 72.5, label = "COLD PAIN", fontface = "bold", size = 5.5) +
  annotate("text", x = 2.2, y = 72.5, label = "COLD ALLODYNIA", fontface = "bold", size = 5.5) +
  labs(y = "cFOS count", x = NULL, fill = NULL, colour = NULL) +
  coord_cartesian(ylim = c(0, 70), clip = "off") +
  guides(colour = "none", fill = guide_legend(ncol = 1, override.aes = list(size = 4))) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 11, colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line = element_line(linewidth = 0.8, colour = "black"),
    legend.position = c(0.5, 1.18),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.text = element_text(size = 8), legend.key.height = grid::unit(0.28, "cm"),
    plot.margin = margin(t = 70, r = 15, b = 10, l = 10))
if (capabilities("aqua")) grDevices::quartz(width = 8, height = 6)
print(CFOS_total)
dir.create("stats", showWarnings = FALSE)
ggsave("stats/cFOS_whole_slice.png", CFOS_total,
  width = 8, height = 6, units = "in", dpi = 600, bg = "white")
