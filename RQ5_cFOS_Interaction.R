library(readxl)

raw <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
raw$MODEL <- trimws(raw$MODEL)
raw$AGENT <- trimws(raw$AGENT)
cfos <- subset(raw, AGENT %in% c("AMG333", "AMG333 + TRPM2-KO"))
cfos$MODEL[cfos$MODEL == "Cold Pain"] <- "CP"
cfos$MODEL[cfos$MODEL == "Cold Allodynia"] <- "OXA"
cfos$TRPM2 <- ifelse(cfos$AGENT == "AMG333", "WT", "KO")

# These percentages are laminar cFOS divided by whole-slice cFOS.
cfos$L1_PERCENT <- ifelse(cfos$CFOS_TOTAL > 0, cfos$L1_CFOS / cfos$CFOS_TOTAL * 100, NA)
cfos$L2_PERCENT <- ifelse(cfos$CFOS_TOTAL > 0, cfos$L2_CFOS / cfos$CFOS_TOTAL * 100, NA)
cfos$L3_PERCENT <- ifelse(cfos$CFOS_TOTAL > 0, cfos$L34_CFOS / cfos$CFOS_TOTAL * 100, NA)

all_effects <- function(values, n_ko) {
  arrangements <- combn(seq_along(values), n_ko, simplify = FALSE)
  sapply(arrangements, function(ko_rows) {
    wt_rows <- setdiff(seq_along(values), ko_rows)
    mean(values[ko_rows]) - mean(values[wt_rows])
  })
}

interaction_test <- function(dat, value_column) {
  d <- data.frame(MODEL = dat$MODEL, TRPM2 = dat$TRPM2, value = dat[[value_column]])
  d <- subset(d, !is.na(value))
  cp_wt <- d$value[d$MODEL == "CP" & d$TRPM2 == "WT"]
  cp_ko <- d$value[d$MODEL == "CP" & d$TRPM2 == "KO"]
  oxa_wt <- d$value[d$MODEL == "OXA" & d$TRPM2 == "WT"]
  oxa_ko <- d$value[d$MODEL == "OXA" & d$TRPM2 == "KO"]

  cp_effect <- mean(cp_ko) - mean(cp_wt)
  oxa_effect <- mean(oxa_ko) - mean(oxa_wt)
  observed <- oxa_effect - cp_effect
  cp_effects <- all_effects(d$value[d$MODEL == "CP"], length(cp_ko))
  oxa_effects <- all_effects(d$value[d$MODEL == "OXA"], length(oxa_ko))

  extreme <- 0
  for (oxa_difference in oxa_effects) {
    extreme <- extreme + sum(abs(oxa_difference - cp_effects) >= abs(observed) - 1e-12)
  }
  p <- extreme / (length(cp_effects) * length(oxa_effects))

  data.frame(
    n_CP_WT = length(cp_wt), n_CP_KO = length(cp_ko),
    n_OXA_WT = length(oxa_wt), n_OXA_KO = length(oxa_ko),
    mean_CP_WT = mean(cp_wt), mean_CP_KO = mean(cp_ko),
    mean_OXA_WT = mean(oxa_wt), mean_OXA_KO = mean(oxa_ko),
    TRPM2_effect_CP = cp_effect, TRPM2_effect_OXA = oxa_effect,
    interaction = observed, p_raw = p)
}

count_results <- data.frame()
for (i in 1:3) {
  test <- interaction_test(cfos, c("L1_CFOS", "L2_CFOS", "L34_CFOS")[i])
  test$endpoint <- "cFOS laminar count"
  test$lamina <- c("L1", "L2", "L3+")[i]
  count_results <- rbind(count_results, test)
}
count_results$p_holm <- p.adjust(count_results$p_raw, method = "holm")

percent_results <- data.frame()
for (i in 1:3) {
  test <- interaction_test(cfos, c("L1_PERCENT", "L2_PERCENT", "L3_PERCENT")[i])
  test$endpoint <- "cFOS laminar % of total"
  test$lamina <- c("L1", "L2", "L3+")[i]
  percent_results <- rbind(percent_results, test)
}
percent_results$p_holm <- p.adjust(percent_results$p_raw, method = "holm")

whole_result <- interaction_test(cfos, "CFOS_TOTAL")
whole_result$endpoint <- "Whole-slice cFOS count"
whole_result$lamina <- "Whole slice"
whole_result$p_holm <- whole_result$p_raw

results <- rbind(count_results, percent_results, whole_result)
results <- results[, c("endpoint", "lamina", "n_CP_WT", "n_CP_KO", "n_OXA_WT", "n_OXA_KO",
  "mean_CP_WT", "mean_CP_KO", "mean_OXA_WT", "mean_OXA_KO", "TRPM2_effect_CP",
  "TRPM2_effect_OXA", "interaction", "p_raw", "p_holm")]

