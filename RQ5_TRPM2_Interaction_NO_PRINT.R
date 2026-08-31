library(readxl)

raw <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
raw$MODEL <- trimws(raw$MODEL)
raw$AGENT <- trimws(raw$AGENT)
raw$MARKER <- trimws(raw$MARKER)
rq5 <- subset(raw, AGENT %in% c("AMG333", "AMG333 + TRPM2-KO"))
rq5$MODEL[rq5$MODEL == "Cold Pain"] <- "CP"
rq5$MODEL[rq5$MODEL == "Cold Allodynia"] <- "OXA"
rq5$TRPM2 <- ifelse(rq5$AGENT == "AMG333", "WT", "KO")

rq5$L1_PERCENT <- ifelse(rq5$L1_CFOS > 0, rq5$L1_TOTAL / rq5$L1_CFOS * 100, NA)
rq5$L2_PERCENT <- ifelse(rq5$L2_CFOS > 0, rq5$L2_TOTAL / rq5$L2_CFOS * 100, NA)
rq5$L3_PERCENT <- ifelse(rq5$L34_CFOS > 0, rq5$L34_TOTAL / rq5$L34_CFOS * 100, NA)
rq5$WHOLE_PERCENT <- ifelse(rq5$CFOS_TOTAL > 0, rq5$SLICE_TOTAL / rq5$CFOS_TOTAL * 100, NA)

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

results <- data.frame()
for (marker_name in c("PAX2", "VGLUT2")) {
  marker_data <- subset(rq5, MARKER == marker_name)

  percent_results <- data.frame()
  for (i in 1:3) {
    test <- interaction_test(marker_data, c("L1_PERCENT", "L2_PERCENT", "L3_PERCENT")[i])
    test$marker <- marker_name
    test$endpoint <- "Marker/cFOS laminar %"
    test$lamina <- c("L1", "L2", "L3+")[i]
    percent_results <- rbind(percent_results, test)
  }
  percent_results$p_holm <- p.adjust(percent_results$p_raw, method = "holm")

  count_results <- data.frame()
  for (i in 1:3) {
    test <- interaction_test(marker_data, c("L1_TOTAL", "L2_TOTAL", "L34_TOTAL")[i])
    test$marker <- marker_name
    test$endpoint <- "Marker laminar count"
    test$lamina <- c("L1", "L2", "L3+")[i]
    count_results <- rbind(count_results, test)
  }
  count_results$p_holm <- p.adjust(count_results$p_raw, method = "holm")

  whole_count <- interaction_test(marker_data, "SLICE_TOTAL")
  whole_count$marker <- marker_name
  whole_count$endpoint <- "Whole-slice marker count"
  whole_count$lamina <- "Whole slice"
  whole_count$p_holm <- whole_count$p_raw

  whole_percent <- interaction_test(marker_data, "WHOLE_PERCENT")
  whole_percent$marker <- marker_name
  whole_percent$endpoint <- "Whole-slice marker/cFOS %"
  whole_percent$lamina <- "Whole slice"
  whole_percent$p_holm <- whole_percent$p_raw

  results <- rbind(results, percent_results, count_results, whole_count, whole_percent)
}

results <- results[, c("marker", "endpoint", "lamina", "n_CP_WT", "n_CP_KO",
  "n_OXA_WT", "n_OXA_KO", "mean_CP_WT", "mean_CP_KO", "mean_OXA_WT", "mean_OXA_KO",
  "TRPM2_effect_CP", "TRPM2_effect_OXA", "interaction", "p_raw", "p_holm")]
