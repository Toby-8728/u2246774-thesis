library(readxl)
library(coin)

data <- read_excel("THESIS COUNT 3.7.xlsx", sheet = "RAW")
data$MODEL <- trimws(data$MODEL)
data$AGENT <- trimws(data$AGENT)
data$MARKER <- trimws(data$MARKER)

data$L1_MARKER_PCT <- ifelse(data$L1_CFOS > 0, data$L1_TOTAL / data$L1_CFOS * 100, NA)
data$L2_MARKER_PCT <- ifelse(data$L2_CFOS > 0, data$L2_TOTAL / data$L2_CFOS * 100, NA)
data$L3_MARKER_PCT <- ifelse(data$L34_CFOS > 0, data$L34_TOTAL / data$L34_CFOS * 100, NA)
data$L1_CFOS_PCT <- ifelse(data$CFOS_TOTAL > 0, data$L1_CFOS / data$CFOS_TOTAL * 100, NA)
data$L2_CFOS_PCT <- ifelse(data$CFOS_TOTAL > 0, data$L2_CFOS / data$CFOS_TOTAL * 100, NA)
data$L3_CFOS_PCT <- ifelse(data$CFOS_TOTAL > 0, data$L34_CFOS / data$CFOS_TOTAL * 100, NA)

run_three_laminae <- function(df, columns, group_col, group1, group2,
  analysis_type, endpoint, marker, context, family) {
  out <- data.frame()
  for (i in 1:3) {
    d <- data.frame(value = df[[columns[i]]], group = df[[group_col]])
    d <- subset(d, !is.na(value) & group %in% c(group1, group2))
    d$group <- factor(d$group, levels = c(group1, group2))
    x <- d$value[d$group == group1]
    y <- d$value[d$group == group2]
    test <- wilcox_test(value ~ group, data = d, distribution = "exact")

    out <- rbind(out, data.frame(
      analysis_type = analysis_type, endpoint = endpoint, marker = marker,
      context = context, lamina = c("L1", "L2", "L3+")[i],
      group1 = group1, n1 = length(x), mean1 = mean(x), median1 = median(x),
      group2 = group2, n2 = length(y), mean2 = mean(y), median2 = median(y),
      p_raw = as.numeric(pvalue(test)), family = family))
  }
  # One Holm family contains L1, L2 and L3+ for this comparison and endpoint.
  out$p_holm <- p.adjust(out$p_raw, method = "holm")
  out
}

results <- data.frame()
for (marker_name in c("PAX2", "VGLUT2", "cFOS")) {
  if (marker_name == "cFOS") {
    marker_data <- data
    count_cols <- c("L1_CFOS", "L2_CFOS", "L34_CFOS")
    percent_cols <- c("L1_CFOS_PCT", "L2_CFOS_PCT", "L3_CFOS_PCT")
    count_endpoint <- "cFOS laminar count"
    percent_endpoint <- "cFOS laminar % of total"
  } else {
    marker_data <- subset(data, MARKER == marker_name)
    count_cols <- c("L1_TOTAL", "L2_TOTAL", "L34_TOTAL")
    percent_cols <- c("L1_MARKER_PCT", "L2_MARKER_PCT", "L3_MARKER_PCT")
    count_endpoint <- "Marker laminar count"
    percent_endpoint <- "Marker/cFOS laminar %"
  }

  for (measure in c("COUNT", "PERCENT")) {
    if (measure == "COUNT") {
      columns <- count_cols
      endpoint <- count_endpoint
    } else {
      columns <- percent_cols
      endpoint <- percent_endpoint
    }

    for (model_name in c("Cold Pain", "Cold Allodynia")) {
      d <- subset(marker_data, MODEL == model_name)
      results <- rbind(results, run_three_laminae(d, columns, "AGENT",
        "AMG333", "AMG333 + TRPM2-KO", "TRPM2 effect", endpoint,
        marker_name, model_name,
        paste("TRPM2", marker_name, measure, model_name, sep = " | ")))

      results <- rbind(results, run_three_laminae(d, columns, "AGENT",
        "TRPM2-KO", "AMG333 + TRPM2-KO", "TRPM8 blockade in TRPM2-KO", endpoint,
        marker_name, model_name,
        paste("M2_vs_M2AMG", marker_name, measure, model_name, sep = " | ")))
    }
    for (agent_name in c("AMG333", "TRPM2-KO", "AMG333 + TRPM2-KO")) {
      d <- subset(marker_data, AGENT == agent_name)
      results <- rbind(results, run_three_laminae(d, columns, "MODEL",
        "Cold Pain", "Cold Allodynia", "Model effect", endpoint,
        marker_name, agent_name,
        paste("CP_vs_OXA", marker_name, measure, agent_name, sep = " | ")))
    }
  }
}

results <- results[order(results$analysis_type, results$marker, results$endpoint,
  results$context, results$family, match(results$lamina, c("L1", "L2", "L3+"))), ]
results_export <- results
for (column in c("mean1", "mean2", "median1", "median2")) {
  results_export[[column]] <- round(results_export[[column]], 3)
}
results_export$p_raw <- signif(results_export$p_raw, 5)
results_export$p_holm <- signif(results_export$p_holm, 5)
notable <- subset(results_export, !is.na(p_holm) & p_holm < 0.10)

