###setup
# required packages
library(dplyr)
library(ggplot2)
library(lmtest)
library(sandwich)
library(ggsci)

options(scipen = 999)

#################
# Generate list of all covariate combinations
covariates = c("DM1", "DM2", "P")
all_combinations = vector("list", length(covariates))
for (i in seq_along(covariates)) {
  all_combinations[[i]] <- combn(covariates, i, simplify = FALSE)
}
covariates = unlist(all_combinations, recursive = FALSE)


################
### Simulation for S-->EMM

# sample(1:10000,1)
seed_start = 2765
iterations = 1000
sample_sizes = c(2000,20000,100000)

#strength of collider association
soca = c("Strong" = 1, "Medium" = 2/3, "Weak" = 1/3)

#strength of interaction
soi = c("Strong" = 1, "Medium" = 2/3, "Weak" = 1/3)

### Monte Carlo simulation
output = list()
for (n in sample_sizes){
  for (i in soca){
    for (m in soi){
      for (j in 1:length(covariates)){
        for (t in 1:iterations){

          # data simulation
          seed = seed_start+t*3
          set.seed(seed)
          S = rbinom(n,1,.5)
          A = rbinom(n,1,.5)
          DM1 = rbinom(n,1,.5)
          P = rbinom(n,1,0.55 + DM1*0.3*i - S*0.5*i)
          DM2 = rbinom(n,1,0.2 + S*0.4)
          Y = rbinom(n,1,prob = 0.05 + A*0.1 + DM1*0.1 + DM2*0.1 +
                       A*DM1*0.3*m + A*DM2*0.3)

          # calculate TATE and SATE
          data = cbind(S,Y,A,DM1,DM2,P)
          source = data[data[,"S"] == 1,]
          target = data[data[,"S"] == 0,]

          SATE_mod = lm(Y ~ A, data = as.data.frame(source))
          SATE_mod_summary = coeftest(SATE_mod , vcov = vcovHC(SATE_mod , type = "HC2"))
          SATE = SATE_mod_summary["A","Estimate"]
          SATE_se = SATE_mod_summary["A","Std. Error"]
          SATE_lb = SATE - 1.96 * SATE_se
          SATE_ub = SATE + 1.96 * SATE_se

          TATE_mod = lm(Y ~ A, data = as.data.frame(target))
          TATE_mod_summary = coeftest(TATE_mod , vcov = vcovHC(TATE_mod , type = "HC2"))
          TATE = TATE_mod_summary["A","Estimate"]
          TATE_se = TATE_mod_summary["A","Std. Error"]
          TATE_lb = TATE - 1.96 * TATE_se
          TATE_ub = TATE + 1.96 * TATE_se

          data = select(as.data.frame(data), "S","Y","A",paste(covariates[[j]], sep = ","))
          covariates_plus = paste(covariates[[j]], collapse = "+")

          # effect estimations (outcome modelling based on code by Dahabreh et al. 2020 (https://doi.org/10.1002/sim.8426))
          # link to source code: https://github.com/serobertson/ExtendingInferences/blob/master/simulated_example_source_code4.R
          S1data_A1 = filter(data, S==1 & A==1)
          OM1mod = glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A1, family = "gaussian")
          p1 = predict(OM1mod,newdata=data, type="response")
          S1data_A0 = filter(data, S==1 & A==0)
          OM0mod = glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A0, family = "gaussian")
          p0 = predict(OM0mod,newdata=data, type="response")
          data$p1 = p1
          data$p0 = p0
          S0sub = subset(data, S==0)
          OM_1 = mean(S0sub$p1)
          OM_0 = mean(S0sub$p0)
          OM = mean(S0sub$p1)-mean(S0sub$p0)

          # estimate agreement of iteration
          ea_SATE = (OM >= SATE_lb & OM <= SATE_ub)
          ea_TATE = (OM >= TATE_lb & OM <= TATE_ub)

          # prepare output
          outtable=data.frame(covariates = paste(covariates[[j]], collapse = ", "), estimate = OM,
                              SATE, SATE_se, SATE_lb, SATE_ub, ea_SATE,
                              TATE, TATE_se, TATE_lb,TATE_ub, ea_TATE,
                              n, soca = i, soi = m, seed)
          output=append(output,list(outtable))
        }
      }
    }
  }
}

# generate results from Monte Carlo simulation
df_output = do.call(rbind, output)

df_output_means = df_output %>%
  group_by(covariates, n, soca , soi) %>%
  summarise(across(c(SATE, TATE, estimate), \(x) round(mean(x), 4), .names = "mean_{.col}"),
            across(c(SATE, TATE, estimate), \(x) round(sd(x), 4), .names = "{.col}_se"),
            across(c(ea_SATE, ea_TATE), \(x) 100* (sum(x)/iterations), .names = "{.col}_perc"),
            across(c(percent_bias), \(x) mean(x), .names = "{.col}"))

df_output_means = mutate(df_output_means,
                         percent_bias = round(abs(100*(mean_TATE - mean_estimate) /
                                              (mean_TATE - mean_SATE)),2))

df_output_means = mutate(df_output_means, sufficient  = case_when(covariates %in% c("DM2","DM2, P") ~ "S-dependent",
                                                                  covariates %in% c("DM1","DM1, P", "P") ~ "Insufficient",
                                                                  TRUE ~ "Sufficient"))

df_output_means$sufficient = factor(df_output_means$sufficient, levels = c("Sufficient",
                                                                           "S-dependent",
                                                                           "Insufficient"))
df_output_means$case = "S-->EMM"


# results for illustrating s-dependence for S-->EMM
arrange(data.frame(unique(select(filter(df_output_means, n == 20000, soi == 1, soca == 1),
                                        covariates, mean_estimate, estimate_se, mean_SATE, mean_TATE,
                                        ea_SATE_perc, ea_TATE_perc, percent_bias, sufficient))), sufficient)

# results for illustrating bias magnitude depending on its determinants
arrange(data.frame(unique(select(filter(df_output_means, covariates == "DM2, P"),covariates, mean_estimate, estimate_se,
                                mean_SATE, mean_TATE,ea_SATE_perc, ea_TATE_perc, percent_bias, sufficient))), sufficient)



###########################################
### Simulation for EMM-->S

# sample(1:10000,1)
seed_start = 1735
iterations = 1000
n = 20000

### Monte-carlo simulation
output2 = list()
for (j in 1:length(covariates)){
  for (t in 1:iterations){

    seed = seed_start+t*3
    set.seed(seed)
    A=rbinom(n,1,.5)
    DM2=rbinom(n,1,0.5)
    DM1=rbinom(n,1,0.4)
    P=rbinom(n,1,0.2+0.3*DM1)
    S=rbinom(n,1, 0.25 + DM2*0.25 +  P*0.4)
    Y = rbinom(n,1,prob = 0.05 + A*0.1 + DM1*0.1 + DM2*0.1 +
                 A*DM1*0.3 + A*DM2*0.3)

    # calculate TATE and SATE
    data = cbind(S,Y,A,DM1,DM2,P)
    source = data[data[,"S"] == 1,]
    target = data[data[,"S"] == 0,]

    SATE_mod = lm(Y ~ A, data = as.data.frame(source))
    SATE_mod_summary = coeftest(SATE_mod , vcov = vcovHC(SATE_mod , type = "HC2"))
    SATE = SATE_mod_summary["A","Estimate"]
    SATE_se = SATE_mod_summary["A","Std. Error"]
    SATE_lb = SATE - 1.96 * SATE_se
    SATE_ub = SATE + 1.96 * SATE_se

    TATE_mod = lm(Y ~ A, data = as.data.frame(target))
    TATE_mod_summary = coeftest(TATE_mod , vcov = vcovHC(TATE_mod , type = "HC2"))
    TATE = TATE_mod_summary["A","Estimate"]
    TATE_se = TATE_mod_summary["A","Std. Error"]
    TATE_lb = TATE - 1.96 * TATE_se
    TATE_ub = TATE + 1.96 * TATE_se

    data = select(as.data.frame(data), "S","Y","A",paste(covariates[[j]], sep = ","))
    covariates_plus = paste(covariates[[j]], collapse = "+")

    # effect estimations (outcome modelling based on code by Dahabreh et al. 2020 (https://doi.org/10.1002/sim.8426))
    # link to source code: https://github.com/serobertson/ExtendingInferences/blob/master/simulated_example_source_code4.R
    S1data_A1 = filter(data, S==1 & A==1)
    OM1mod = glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A1, family = "gaussian")
    p1 = predict(OM1mod,newdata=data, type="response")
    S1data_A0 = filter(data, S==1 & A==0)
    OM0mod = glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A0, family = "gaussian")
    p0 = predict(OM0mod,newdata=data, type="response")
    data$p1 = p1
    data$p0 = p0
    S0sub = subset(data, S==0)
    OM_1 = mean(S0sub$p1)
    OM_0 = mean(S0sub$p0)
    OM = mean(S0sub$p1)-mean(S0sub$p0)

    # estimate agreement of iteration
    ea_SATE = (OM >= SATE_lb & OM <= SATE_ub)
    ea_TATE = (OM >= TATE_lb & OM <= TATE_ub)

    # prepare output
    outtable=data.frame(covariates = paste(covariates[[j]], collapse = ", "), estimate = OM,
                        SATE, SATE_se, SATE_lb, SATE_ub, ea_SATE,
                        TATE, TATE_se, TATE_lb,TATE_ub, ea_TATE,
                        n, soca = 1, soi = 1, seed)
    output2=append(output2,list(outtable))
  }
}

df_output2 = do.call(rbind, output2)

df_output2_means = df_output2 %>%
  group_by(covariates) %>%
  summarise(across(c(SATE, TATE, estimate), \(x) round(mean(x), 4), .names = "mean_{.col}"),
            across(c(SATE, TATE, estimate), \(x) round(sd(x), 4), .names = "{.col}_se"),
            across(c(ea_SATE, ea_TATE), \(x) 100* (sum(x)/iterations), .names = "{.col}_perc"),
            across(c(percent_bias), \(x) mean(x), .names = "{.col}"))
data.frame(df_output2_means)

df_output2_means = mutate(df_output2_means,
                         percent_bias = round(abs(100*(mean_TATE - mean_estimate) /
                                                    (mean_TATE - mean_SATE)),2))

df_output2_means = mutate(df_output2_means, sufficient  = case_when(covariates %in% c("DM2","DM2, P") ~ "S-dependent",
                                                                  covariates %in% c("DM1","DM1, P", "P") ~ "Insufficient",
                                                                  TRUE ~ "Sufficient"))

df_output2_means$sufficient = factor(df_output2_means$sufficient, levels = c("Sufficient",
                                                                           "S-dependent",
                                                                           "Insufficient"))

df_output2_means$case = "EMM-->S"

# results for illustrating s-dependence for EMM-->S
arrange(data.frame(unique(select(df_output2_means,covariates, mean_estimate, estimate_se,
                                 mean_SATE, mean_TATE,
                                 ea_SATE_perc, ea_TATE_perc, percent_bias, sufficient))), sufficient)


#########################
### create density plot for supplement

df_output$case = "S-->EMM"
df_output2$case = "EMM-->S"

df_output_overall = rbind(df_output, df_output2)
df_output_overall = mutate(df_output_overall,
                           sufficient  = case_when(covariates %in% c("DM2","DM2, P") ~ "S-dependent",
                                                   covariates %in% c("DM1","DM1, P", "P") ~ "Insufficient",
                                                   TRUE ~ "Sufficient"))
covariate_order <- df_output_overall %>%
  distinct(covariates, sufficient) %>%
  mutate(sufficient = factor(sufficient, levels = c("Sufficient", "S-dependent", "Insufficient"))) %>%
  arrange(sufficient, covariates) %>%   # secondary sort within each group by name;
  pull(covariates)                      # swap 'covariates' for any other within-group order you prefer

df_output_overall <- df_output_overall %>%
  mutate(covariates = factor(covariates, levels = covariate_order))

# plot
ggplot(filter(df_output_overall, n == 20000, soi == 1, soca == 1)) +
  geom_density(aes(x = estimate, fill = "Estimate")) +
  geom_density(aes(x = SATE, fill = "SATE")) +
  geom_density(aes(x = TATE, fill = "TATE")) +
  scale_fill_manual(name = NULL,
                    values = c("Estimate" = pal_lancet("lanonc", alpha = 0.3)(9)[1],
                               "TATE"     = pal_lancet("lanonc", alpha = 0.3)(9)[2],
                               "SATE"     = pal_lancet("lanonc", alpha = 0.3)(9)[3])) +
  scale_x_continuous(name = NULL, breaks = c(0.3, 0.4, 0.5), expand = c(0, 0)) +
  scale_y_continuous("Density", expand = c(0, 0)) +

  facet_grid(case ~ covariates, axes = "all_x",
             labeller = labeller(case = c("EMM-->S" = "EMMs cause selection",
                                          "S-->EMM" = "Selection causing EMMs"))) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.justification = c(0.5,0),
        legend.margin    = margin(t = -5),
        plot.margin      = margin(b = 2, l = 5, t = 2, r = 2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        )
