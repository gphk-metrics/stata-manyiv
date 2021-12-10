#!/usr/bin/env Rscript

library(ggplot2)
library(stringr)
library(haven)
library(tidyr)
library(dplyr)
library(readr)

maindir  <- ".."
progmsg  <- "Graph Coefficient Comparison"

Main <- function () {
    print(version)
    cat(sep = "\n", paste(rep('-', 72), collapse = ''), date(), "", progmsg)

    coef_input   <- import_data()
    graph_input  <- reshape_data(coef_input)
    graph_output <- comparison_plot(graph_input)
    ggsave(file.path(maindir, "misc/graph_output.pdf"),
           graph_output, "pdf", width = 12, height = 14)

    cat(sep = "\n", paste(rep('-', 72), collapse = ''), date(), "", progmsg)
}

# ---------------------------------------------------------------------
# Aux functions

import_data <- function () {
    read_csv(file.path(maindir, "misc/graph_input.csv"))
}

reshape_data <- function(coef_input) {
    varlabel <- c("Financial Strain Index",
                  "Revolving Balance",
                  "Collection Balance",
                  "Have a Mortgage",
                  "Mortgage Balance",
                  "Have an Auto Loan",
                  "Auto Balance",
                  "Revolving Utilization",
                  "Non-Mortgage Inquiries")

    graph_input <- coef_input %>%
        bind_cols(tibble(var_label = varlabel)) %>%
        mutate(b_ols        = beta_ols,
               sig_ols      = se_ols,
               b_ujive_w    = beta_ujive_fe_w,
               b_ujive_no_w = beta_ujive_fe_no_w) %>%
        gather(variable, estimate, -varname, -var_label, -b_ols, -sig_ols, -sd, -b_ujive_w, -b_ujive_no_w, factor_key = TRUE) %>%
        separate(variable, into=c("coef", "type"), extra = "merge") %>%
        spread(coef, estimate) %>%
        mutate(interacted    = !str_detect(type, "no_w"),
               judge_fe      = str_detect(type, "fe_"),
               ujive         = str_detect(type, "ujive"),
               paper         = str_detect(type, "paper"),
               iv            = str_detect(type, "iv_"),
               ols           = str_detect(type, "ols"),
               jive          = !ujive & !iv  & ! ols,
               estimate_type = case_when(ujive ~ "UJIVE",
                                         paper ~ "UJIVE-Paper",
                                         jive  ~ "JIVE",
                                         iv    ~ "IV",
                                         ols   ~ "OLS"),
               shape_type = case_when(!ols &  judge_fe &  interacted ~ "Judge FE x W",
                                      !ols &  judge_fe & !interacted ~ "Judge FE",
                                      !ols & !judge_fe &  interacted ~ "Leniency using W",
                                      !ols & !judge_fe & !interacted ~ "Leniency",
                                       ols ~ "OLS"),
               shape_type2 = case_when(!ols &  judge_fe ~ "Judge FE",
                                       !ols & !judge_fe ~ "Leniency",
                                       ols ~ "OLS"))

    graph_input <- graph_input %>%
        bind_rows(graph_input %>%
                  filter(estimate_type == "OLS" | estimate_type == "UJIVE-Paper") %>%
                  mutate(interacted = FALSE)) %>%
        mutate(beta_diff = (beta - b_ols)/sig_ols) %>%
        mutate(t2 = beta/sig_ols) %>%
        mutate(scaled_beta = case_when(interacted == TRUE  ~ (beta - b_ujive_w)/sd,
                                       interacted == FALSE ~ (beta - b_ujive_w)/sd)) %>%
        mutate(scaled_se   = case_when(interacted == TRUE  ~ se/sd,
                                       interacted == FALSE ~ se/sd)) %>%
        mutate(var_label = factor(var_label, levels = varlabel),
        interacted = factor(interacted, labels = c("Uninteracted", "Interacted")))

    x <- c("UJIVE-Paper", "UJIVE", "JIVE", "IV", "OLS")
    graph_input$estimate_factor <-  factor(graph_input$estimate_type, x)
    return(graph_input)
}

comparison_plot <- function (graph_input) {
    theme <- theme_bw() +
        theme(panel.border     = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.caption     = element_text(hjust = 0, size = 10),
              axis.line        = element_line(colour = "black"),
              plot.title       = element_text(size = 24),
              strip.text.x     = element_text(size = 14),
              legend.title     = element_text(size = 16),
              legend.text      = element_text(size = 14),
              axis.title.x     = element_text(size = 16),
              axis.title.y     = element_text(size = 16),
              axis.text.x      = element_text(size = 14),
              axis.text.y      = element_text(size = 14))

    ggplot() +
        geom_pointrange(data = graph_input,
                        aes(y     = scaled_beta,
                            ymin  = scaled_beta - 1.96 * scaled_se,
                            ymax  = scaled_beta + 1.96 * scaled_se,
                            x     = var_label,
                            color = estimate_factor),
                        show.legend = FALSE,
                        alpha = 0.35,
                        size = 1.25,
                        fatten = 0,
                        position = position_dodge2(width=2)) +
        geom_point(data = graph_input,
                   aes(y     = scaled_beta,
                       x     = var_label,
                       color = estimate_factor,
                       shape = shape_type2),
                   size = 4.5,
                   position = position_dodge(width=2)) +
        scale_shape_manual(values = c(15, 16, 17)) +
        scale_x_discrete(breaks = graph_input$var_label,
                         limits = c(t(cbind(levels(graph_input$var_label), "skip", "skip")))) +
        coord_flip(ylim = c(-1.5, 1.5)) +
        labs(x     = "",
             y     = "(Estimated Coefficient - Judge FE x W UJIVE)/Control SD",
             shape = "Estimation Approach",
             color = "Leave-out Approach") +
        scale_color_manual(values = c("black", gg_color_hue(4)),
                           breaks = c("OLS", "IV", "JIVE", "UJIVE", "UJIVE-Paper")) +
        scale_shape_discrete(breaks = c("OLS", "Leniency", "Judge FE")) +
        facet_wrap(~interacted) +
        geom_hline(yintercept = 0) +
        theme
}

gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

# ---------------------------------------------------------------------
# Run the things

Main()
