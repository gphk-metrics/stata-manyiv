

dblue = rgb(0,114,178,maxColorValue = 255)
dred = rgb(213,94,0,maxColorValue = 255)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
library(readr)
graph_input <- read_csv("~/Dropbox/Papers/GPH_ExaminerDesign/Code/stata-manyiv/misc/graph_input.csv")
varlabel = c("Financial Strain Index", 
              "Revolving Balance", 
              "Collection Balance",
              "Have a Mortgage",
              "Mortgage Balance",
              "Have an Auto Loan",
              "Auto Balance",
              "Revolving Utilization",
              "Non-Mortgage Inquiries")


              

graph_input2 = graph_input %>% 
  bind_cols(tibble(var_label=varlabel)) %>% 
  rename(beta_ols = output1,
         se_ols   = output2,
         beta_paper_no_w = output3,
         se_paper_no_w = output4) %>%
  mutate(b_ols = beta_ols, sig_ols = se_ols,
         b_ujive_w = beta_ujive_fe_w,
         b_ujive_no_w = beta_ujive_fe_no_w) %>%
  gather(variable, estimate, -varname, -var_label, -b_ols, -sig_ols, -sd, -b_ujive_w, -b_ujive_no_w, factor_key = TRUE) %>% 
  separate(variable, into=c("coef", "type"), extra="merge") %>%
  spread(coef, estimate) %>%
  mutate(interacted = !str_detect(type, "no_w"),
         judge_fe = str_detect(type, "fe_"),
         ujive = str_detect(type, "ujive"),
         paper = str_detect(type, "paper"),
         iv = str_detect(type, "iv_"),
         ols = str_detect(type, "ols"),
         jive = !ujive & !iv  & ! ols,
         estimate_type = case_when(ujive ~ "UJIVE",
                                   paper ~ "UJIVE-Paper",
                                   jive ~ "JIVE",
                                   iv ~ "IV",
                                   ols ~ "OLS"),
         shape_type = case_when(!ols & judge_fe & interacted  ~ "Judge FE x W",
                                !ols & judge_fe & !interacted  ~ "Judge FE",
                                !ols & !judge_fe & interacted  ~ "Leniency using W",
                                !ols & !judge_fe & !interacted  ~ "Leniency",
                                ols ~ "OLS"
                                ),
         shape_type2 = case_when(!ols & judge_fe  ~ "Judge FE",
                                !ols & !judge_fe ~ "Leniency",
                                ols ~ "OLS"
         ))

graph_input2 = graph_input2 %>% bind_rows(
  graph_input2 %>% filter(estimate_type == "OLS" | estimate_type == "UJIVE-Paper") %>%
  mutate(interacted = FALSE)
  )  %>%
  mutate(beta_diff = (beta - b_ols)/sig_ols) %>%
  mutate(t2 = beta/sig_ols) %>%
  mutate(scaled_beta = case_when(interacted == TRUE ~ (beta - b_ujive_w)/sd,
                                 interacted == FALSE ~ (beta - b_ujive_w)/sd)) %>%
           mutate(var_label = factor(var_label, levels=varlabel),
                  interacted = factor(interacted, labels =c("Uninteracted", "Interacted")))



ggplot() +
  geom_point(data = graph_input2, 
             aes(y = scaled_beta, x = var_label,
                 color = estimate_type,
                 shape = shape_type2), size = 3) + coord_flip() +
  labs(x ="", y = "(Estimated Coefficient - Judge FE x W UJIVE)/Control SD",
       shape = "Estimation Approach",
       color = "Leave-out Approach") +
  scale_color_manual(values = c("black",gg_color_hue(4)), 
                     breaks=c("OLS", "IV", "JIVE", "UJIVE", "UJIVE-Paper")) +
  scale_shape_discrete(breaks=c("OLS", "Leniency", "Judge FE")) +
  facet_wrap(~interacted) +
  theme_classic() +
  geom_hline(yintercept = 0) 





ggplot() +
  geom_pointrange(data = graph_input2 %>% filter(var_label=="Financial Strain Index"), 
             aes(y = beta, ymin = beta-1.96*se, ymax = beta+1.96*se, x = var_label, 
                 color = estimate_type,
                 shape = shape_type), position=position_dodge2(width=0.4)) + coord_flip() 

ggplot() +
  geom_pointrange(data = graph_input2, 
                  aes(y = beta-b_ols, ymin = (beta-b_ols)-1.96*se, ymax = (beta--b_ols)+1.96*se, x = var_label, 
                      color = estimate_type,
                      shape = shape_type), position=position_dodge2(width=1)) + coord_flip() +
  ylim(-2,2)

