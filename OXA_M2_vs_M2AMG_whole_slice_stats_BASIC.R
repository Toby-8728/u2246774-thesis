library(readxl)
library(dplyr)
library(coin)

raw <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
raw$MODEL <- trimws(raw$MODEL)
raw$AGENT <- trimws(raw$AGENT)
raw$MARKER <- trimws(raw$MARKER)
oxa <- subset(raw, MODEL == "Cold Allodynia" &
  AGENT %in% c("TRPM2-KO", "AMG333 + TRPM2-KO"))
oxa$AGENT <- factor(oxa$AGENT, levels = c("TRPM2-KO", "AMG333 + TRPM2-KO"))

# PAX2
pax2 <- subset(oxa, MARKER == "PAX2" & !is.na(SLICE_TOTAL))
pax2_test <- wilcox_test(SLICE_TOTAL ~ AGENT, data = pax2, distribution = "exact")
pax2_summary <- pax2 %>%
  group_by(AGENT) %>%
  summarise(n = n(), mean = mean(SLICE_TOTAL), median = median(SLICE_TOTAL),
    .groups = "drop")

# VGLUT2
vglut2 <- subset(oxa, MARKER == "VGLUT2" & !is.na(SLICE_TOTAL))
vglut2_test <- wilcox_test(SLICE_TOTAL ~ AGENT, data = vglut2, distribution = "exact")
vglut2_summary <- vglut2 %>%
  group_by(AGENT) %>%
  summarise(n = n(), mean = mean(SLICE_TOTAL), median = median(SLICE_TOTAL),
    .groups = "drop")

# cFOS uses the CFOS_TOTAL values from both staining groups.
cfos <- subset(oxa, !is.na(CFOS_TOTAL))
cfos_test <- wilcox_test(CFOS_TOTAL ~ AGENT, data = cfos, distribution = "exact")
cfos_summary <- cfos %>%
  group_by(AGENT) %>%
  summarise(n = n(), mean = mean(CFOS_TOTAL), median = median(CFOS_TOTAL),
    .groups = "drop")

results <- data.frame(
  comparison = "OXM2 vs OXM2AMG",
  marker = c("PAX2", "VGLUT2", "cFOS"),
  p_raw = c(pvalue(pax2_test), pvalue(vglut2_test), pvalue(cfos_test)))

results$p_holm <- results$p_raw
group_means <- bind_rows(PAX2 = pax2_summary, VGLUT2 = vglut2_summary,
  cFOS = cfos_summary, .id = "marker")

print(results, row.names = FALSE)
print(group_means)

