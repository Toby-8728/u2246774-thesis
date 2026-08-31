library(readxl)
library(coin)
library(dplyr)

data <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
data$AGENT <- factor(trimws(data$AGENT))
data$MODEL <- factor(trimws(data$MODEL))


AMGCPOXA <- wilcox_test(CFOS_TOTAL ~ MODEL,
  data = droplevels(subset(data, AGENT == "AMG333")),
  distribution = "exact")

TRPM2CP <- wilcox_test(CFOS_TOTAL ~ AGENT,
  data = droplevels(subset(data, MODEL == "Cold Pain" &
    AGENT %in% c("AMG333", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

TRPM2OXA <- wilcox_test(CFOS_TOTAL ~ AGENT,
  data = droplevels(subset(data, MODEL == "Cold Allodynia" &
    AGENT %in% c("AMG333", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

p <- c(pvalue(AMGCPOXA), pvalue(TRPM2CP), pvalue(TRPM2OXA))
results <- data.frame(
  comparison = c("CPAMG vs OXAMG", "CPAMG vs CPM2AMG", "OXAMG vs OXM2AMG"),
  p_raw = p,
  p_holm = p.adjust(p, method = "holm"),
  


TRPM8CP <- wilcox_test(CFOS_TOTAL ~ AGENT,
  data = droplevels(subset(data, MODEL == "Cold Pain" &
    AGENT %in% c("TRPM2-KO", "AMG333 + TRPM2-KO"))),
  distribution = "exact")

M2CPOXA <- wilcox_test(CFOS_TOTAL ~ MODEL,
  data = droplevels(subset(data, AGENT == "TRPM2-KO")),
  distribution = "exact")

M2AMGCPOXA <- wilcox_test(CFOS_TOTAL ~ MODEL,
  data = droplevels(subset(data, AGENT == "AMG333 + TRPM2-KO")),
  distribution = "exact")

added <- data.frame(
  comparison = c("CPM2 vs CPM2AMG", "CPM2 vs OXM2", "CPM2AMG vs OXM2AMG"),
  p_raw = c(pvalue(TRPM8CP), pvalue(M2CPOXA), pvalue(M2AMGCPOXA)),
  p_holm = NA_real_,
results <- rbind(results, added)

group_means <- data %>%
  filter(!is.na(CFOS_TOTAL)) %>%
  group_by(MODEL, AGENT) %>%
  summarise(n = n(), mean = mean(CFOS_TOTAL), median = median(CFOS_TOTAL),
    .groups = "drop")

print(results, row.names = FALSE)
print(group_means)

