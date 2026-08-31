library(readxl)
library(coin)
library(dplyr)

data <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
data$AGENT <- factor(trimws(data$AGENT))
data$MODEL <- factor(trimws(data$MODEL))
data$MARKER <- trimws(data$MARKER)
vg <- subset(data, MARKER == "VGLUT2")
p2 <- subset(data, MARKER == "PAX2")

VG_CP_TRPM2 <- wilcox_test(SLICE_TOTAL ~ AGENT,
  data = droplevels(subset(vg, MODEL == "Cold Pain" &
    AGENT %in% c("AMG333", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

VG_CP_OXA <- wilcox_test(SLICE_TOTAL ~ MODEL,
  data = droplevels(subset(vg, AGENT == "AMG333")),
  distribution = "exact")

VG_OXA_TRPM2 <- wilcox_test(SLICE_TOTAL ~ AGENT,
  data = droplevels(subset(vg, MODEL == "Cold Allodynia" &
    AGENT %in% c("AMG333", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

vg_p <- c(pvalue(VG_CP_TRPM2), pvalue(VG_CP_OXA), pvalue(VG_OXA_TRPM2))
VG_results <- data.frame(
  comparison = c("CPAMG vs CPM2AMG", "CPAMG vs OXAMG", "OXAMG vs OXM2AMG"),
  p_raw = vg_p,
  p_holm = p.adjust(vg_p, method = "holm"),

VG_CP_TRPM8 <- wilcox_test(SLICE_TOTAL ~ AGENT,
  data = droplevels(subset(vg, MODEL == "Cold Pain" &
    AGENT %in% c("TRPM2-KO", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

VG_M2_CP_OXA <- wilcox_test(SLICE_TOTAL ~ MODEL,
  data = droplevels(subset(vg, AGENT == "TRPM2-KO")),
  distribution = "exact")

VG_M2AMG_CP_OXA <- wilcox_test(SLICE_TOTAL ~ MODEL,
  data = droplevels(subset(vg, AGENT == "AMG333 + TRPM2-KO")),
  distribution = "exact")

vg_added <- data.frame(
  comparison = c("CPM2 vs CPM2AMG", "CPM2 vs OXM2", "CPM2AMG vs OXM2AMG"),
  p_raw = c(pvalue(VG_CP_TRPM8), pvalue(VG_M2_CP_OXA), pvalue(VG_M2AMG_CP_OXA)),
  p_holm = NA_real_,
VG_results <- rbind(VG_results, vg_added)

# PAX2 is tested separately, not against VGLUT2.
P2_CP_TRPM2 <- wilcox_test(SLICE_TOTAL ~ AGENT,
  data = droplevels(subset(p2, MODEL == "Cold Pain" &
    AGENT %in% c("AMG333", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

P2_CP_OXA <- wilcox_test(SLICE_TOTAL ~ MODEL,
  data = droplevels(subset(p2, AGENT == "AMG333")),
  distribution = "exact")

P2_OXA_TRPM2 <- wilcox_test(SLICE_TOTAL ~ AGENT,
  data = droplevels(subset(p2, MODEL == "Cold Allodynia" &
    AGENT %in% c("AMG333", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

p2_p <- c(pvalue(P2_CP_TRPM2), pvalue(P2_CP_OXA), pvalue(P2_OXA_TRPM2))
P2_results <- data.frame(
  comparison = c("CPAMG vs CPM2AMG", "CPAMG vs OXAMG", "OXAMG vs OXM2AMG"),
  p_raw = p2_p,
  p_holm = p.adjust(p2_p, method = "holm"),

.
P2_CP_TRPM8 <- wilcox_test(SLICE_TOTAL ~ AGENT,
  data = droplevels(subset(p2, MODEL == "Cold Pain" &
    AGENT %in% c("TRPM2-KO", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

P2_M2_CP_OXA <- wilcox_test(SLICE_TOTAL ~ MODEL,
  data = droplevels(subset(p2, AGENT == "TRPM2-KO")),
  distribution = "exact")

P2_M2AMG_CP_OXA <- wilcox_test(SLICE_TOTAL ~ MODEL,
  data = droplevels(subset(p2, AGENT == "AMG333 + TRPM2-KO")),
  distribution = "exact")

p2_added <- data.frame(
  comparison = c("CPM2 vs CPM2AMG", "CPM2 vs OXM2", "CPM2AMG vs OXM2AMG"),
  p_raw = c(pvalue(P2_CP_TRPM8), pvalue(P2_M2_CP_OXA), pvalue(P2_M2AMG_CP_OXA)),
  p_holm = NA_real_,
P2_results <- rbind(P2_results, p2_added)

group_means <- data %>%
  filter(MARKER %in% c("VGLUT2", "PAX2"), !is.na(SLICE_TOTAL)) %>%
  group_by(MARKER, MODEL, AGENT) %>%
  summarise(n = n(), mean = mean(SLICE_TOTAL), median = median(SLICE_TOTAL),
    .groups = "drop")

print(VG_results, row.names = FALSE)
print(P2_results, row.names = FALSE)
print(group_means)


dir.create("stats2", showWarnings = FALSE)
write.csv(VG_results, "stats2/VGLUT2_whole_slice_stats.csv", row.names = FALSE)
write.csv(P2_results, "stats2/PAX2_whole_slice_stats.csv", row.names = FALSE)
write.csv(group_means, "stats2/marker_whole_slice_means.csv", row.names = FALSE)
