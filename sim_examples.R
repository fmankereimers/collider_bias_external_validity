##############################
### required packages
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggh4x)
library(scales)

##############################
### Generate all covariate combinations
covariates = c("DM1", "DM2", "P")
all_combinations = vector("list", length(covariates))
for (i in seq_along(covariates)) {
  all_combinations[[i]] <- combn(covariates, i, simplify = FALSE)
}
covariates = unlist(all_combinations, recursive = FALSE)


##############################
### simulation for S --> EMM

bootstraps = 1000

## determine levels of sample size and bias parameters for multidimensional bias analysis
sample_sizes = c(2000, 20000, 100000)

#strength of collider association
soca = c("Strong" = 1, "Medium" = 2/3, "Weak" = 1/3)

#strength of interaction
soi = c("Strong" = 1, "Medium" = 2/3, "Weak" = 1/3)

## boostrap procedure
output = list()
for (n in sample_sizes){
  for (i in soca){
    for (m in soi){
      for (j in 1:length(covariates)){
        for (t in 1:bootstraps){
          
          # data simulation
          seed = 2765+t*3
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
          SATE = mean(source[source[,"A"]==1,"Y"]) - mean(source[source[,"A"]==0,"Y"])
          target = data[data[,"S"] == 0,]
          TATE = mean(target[target[,"A"]==1,"Y"]) - mean(target[target[,"A"]==0,"Y"])
          data = select(as.data.frame(data), "S","Y","A",paste(covariates[[j]], sep = ","))
          covariates_plus = paste(covariates[[j]], collapse = "+")
          
          # effect estimations (outcome modelling based on code by Dahabreh et al. 2020 (https://doi.org/10.1002/sim.8426))
          # link to source code: https://github.com/serobertson/ExtendingInferences/blob/master/simulated_example_source_code4.R
          S1data_A1 = filter(data, S==1 & A==1)
          OM1mod = glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A1, family = "binomial")
          p1 = predict(OM1mod,newdata=data, type="response")
          S1data_A0 = filter(data, S==1 & A==0)
          OM0mod = glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A0, family = "binomial")
          p0 = predict(OM0mod,newdata=data, type="response")
          data$p1 = p1
          data$p0 = p0
          S0sub = subset(data, S==0)
          OM_1 = mean(S0sub$p1)
          OM_0 = mean(S0sub$p0)
          OM = mean(S0sub$p1)-mean(S0sub$p0)
          
          # prepare output
          outtable=data.frame(covariates = paste(covariates[[j]], collapse = ", "), estimate = OM, SATE, TATE, n, soca = i, soi = m, seed)
          output=append(output,list(outtable))
        }
      }
    }
  }
}

# transform output list to dataframe
df_output = do.call(rbind, output)

# get means and standard errors across bootstraps within sample size, soca and soi levels
df_output_means = df_output %>%
  group_by(covariates, n, soca , soi) %>%
  summarise(across(c(SATE, TATE, estimate), \(x) mean(x, na.rm = TRUE), .names = "mean_{.col}"),
            across(c(SATE, TATE, estimate), \(x) sd(x, na.rm = TRUE), .names = "{.col}_se"))

# calculate 95% confidence intervalls based on bootstrap standard error
df_output_means = mutate(df_output_means,
                         estimate_ci_lb = mean_estimate - 1.96 * estimate_se,
                         estimate_ci_ub = mean_estimate + 1.96 * estimate_se,
                         SATE_ci_lb = mean_SATE - 1.96 * SATE_se,
                         SATE_ci_ub = mean_SATE + 1.96 * SATE_se,
                         TATE_ci_lb = mean_TATE - 1.96 * TATE_se,
                         TATE_ci_ub = mean_TATE + 1.96 * TATE_se)

# label covariate sets
df_output_means = mutate(df_output_means, sufficient  = case_when(covariates %in% c("DM2","DM2, P") ~ "S-dependent",
                                                                  covariates %in% c("DM1","DM1, P", "P") ~ "Insufficient",
                                                                  TRUE ~ "Sufficient"))

# calculate percent bias and estimate agreement
df_output_means = mutate(df_output_means,
                         percent_bias = abs(100*(mean_TATE - mean_estimate) /
                                              (mean_TATE - mean_SATE)),
                         estimate_agreement = case_when(mean_estimate >= TATE_ci_lb & mean_estimate <= TATE_ci_ub ~ "yes",
                                                        mean_estimate < TATE_ci_lb | mean_estimate > TATE_ci_ub ~ "no"))

# define case
df_output_means$case = "S-->EMM"



##############################
### simulation for EMM --> S

bootstraps = 1000

n = 20000

## boostrap procedure
output2 = list()
for (j in 1:length(covariates)){
  for (t in 1:bootstraps){
    
    # data simulation
    seed = 7120+t*3
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
    SATE = mean(source[source[,"A"]==1,"Y"]) - mean(source[source[,"A"]==0,"Y"])
    target = data[data[,"S"] == 0,]
    TATE = mean(target[target[,"A"]==1,"Y"]) - mean(target[target[,"A"]==0,"Y"])
    data = as.data.frame(data)
    data = select(data, "S","Y","A",paste(covariates[[j]], sep = ","))
    covariates_plus = paste(covariates[[j]], collapse = "+")
    
    # effect estimations (outcome modelling based on code by Dahabreh et al. 2020 (https://doi.org/10.1002/sim.8426))
    # link to source code: https://github.com/serobertson/ExtendingInferences/blob/master/simulated_example_source_code4.R
    S1data_A1<-filter(data, S==1 & A==1)
    OM1mod<-glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A1, family = "binomial")
    p1<- predict(OM1mod,newdata=data, type="response")
    S1data_A0<-filter(data, S==1 & A==0)
    OM0mod<-glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A0, family = "binomial")
    p0<- predict(OM0mod,newdata=data, type="response")
    data$p1<-p1
    data$p0<-p0
    S0sub<-subset(data, S==0)
    OM_1<-mean(S0sub$p1)
    OM_0<-mean(S0sub$p0)
    OM<-mean(S0sub$p1)-mean(S0sub$p0)
    
    # prepare output
    outtable=data.frame(covariates = paste(covariates[[j]], collapse = ", "), estimate = OM, SATE, TATE, seed, n)
    output2=append(output2,list(outtable))
  }
}

# transform output list to dataframe
df_output2 = do.call(rbind, output2)

# get means and standard errors
df_output_means2 = df_output2 %>%
  group_by(covariates) %>%
  summarise(across(c(SATE, TATE, estimate), \(x) mean(x, na.rm = TRUE), .names = "mean_{.col}"),
            across(c(SATE, TATE, estimate), \(x) sd(x, na.rm = TRUE), .names = "{.col}_se"))

# calculate 95% confidence intervalls based on bootstrap standard error
df_output_means2 = mutate(df_output_means2,
                          estimate_ci_lb = mean_estimate - 1.96 * estimate_se,
                          estimate_ci_ub = mean_estimate + 1.96 * estimate_se,
                          SATE_ci_lb = mean_SATE - 1.96 * SATE_se,
                          SATE_ci_ub = mean_SATE + 1.96 * SATE_se,
                          TATE_ci_lb = mean_TATE - 1.96 * TATE_se,
                          TATE_ci_ub = mean_TATE + 1.96 * TATE_se)

# label covariate sets
df_output_means2 = mutate(df_output_means2, sufficient  = case_when(covariates %in% c("DM2","DM2, P") ~ "S-dependent",
                                                                    covariates %in% c("DM1","DM1, P", "P") ~ "Insufficient",
                                                                    TRUE ~ "Sufficient"))
# calculate percent bias and estimate agreement
df_output_means2 = mutate(df_output_means2,
                          percent_bias = abs(100*(mean_TATE - mean_estimate) /
                                               (mean_TATE - mean_SATE)),
                          estimate_agreement = case_when(mean_estimate >= TATE_ci_lb & mean_estimate <= TATE_ci_ub ~ "yes",
                                                         mean_estimate < TATE_ci_lb | mean_estimate > TATE_ci_ub ~ "no"))

# assign sample size, soca, soi and case
df_output_means2$n = 20000
df_output_means2$soca = 1
df_output_means2$soi = 1
df_output_means2$case = "EMM-->S"


##############################
### create overview table
results_overview = bind_rows(df_output_means, df_output_means2)

# round numerical results
results_overview[] <- lapply(names(results_overview), function(col) {
  x <- results_overview[[col]]
  if (is.numeric(x) && col %in% names(results_overview)[!names(results_overview) == "percent_bias"]) {
    round(x, 2)
  }
  else if (is.numeric(x) && col == "percent_bias"){
    round(x, 0)
  }
  else {
    x
  }
})

# cosmetic changes
results_overview = mutate(results_overview, percent_bias = paste0(percent_bias,"%"))
results_overview = mutate(results_overview, estimate = paste0(mean_estimate," (",estimate_ci_lb,"-",estimate_ci_ub,")"),
                          TATE = paste0(mean_TATE," (",TATE_ci_lb,"-",TATE_ci_ub,")"),
                          SATE = paste0(mean_SATE," (",SATE_ci_lb,"-",SATE_ci_ub,")"))
results_overview = select(results_overview, covariates, case, n, soca, soi, estimate, TATE, SATE, sufficient, percent_bias, estimate_agreement)


##############
### tables for supplement

# illustrating s-dependence
etable1 = filter(results_overview, n == 20000, soi == 1, soca == 1)
etable1 = select(ungroup(etable1), case, covariates, sufficient, estimate, TATE, SATE, percent_bias, estimate_agreement)

# bias analysis
etable2 = filter(results_overview, case == "S-->EMM", covariates == "DM2, P")
etable2 = select(ungroup(etable2), n, soca, soi, estimate, TATE, SATE, percent_bias, estimate_agreement)


##############################
### plots

## Figure 4
# prepare plot data
figure4 = bind_rows(df_output_means, df_output_means2)
figure4 = mutate(figure4, case2 = case_when(case == "S-->EMM" ~ "Selection Causing EMMs",
                                            TRUE ~ "EMMs Cause Selection"))
figure4$case2 = factor(figure4$case2,
                       levels = c("Selection Causing EMMs", "EMMs Cause Selection"))
figure4$sufficient = factor(figure4$sufficient,
                            levels = c("Sufficient", "S-dependent", "Insufficient"))
figure4 = filter(figure4, n == 20000, soi == 1, soca == 1)
figure4 = arrange(figure4, sufficient)
figure4$covariates = factor(figure4$covariates,
                            levels = unique(figure4$covariates[order(figure4$sufficient)]))
# plot
figure4_plot = ggplot(figure4) +
  geom_rect(aes(xmin = SATE_ci_lb, xmax = SATE_ci_ub, ymin = 0,
                ymax = length(unique(figure4$covariates))+1),
            fill = pal_lancet("lanonc")(9)[8], color = NA, alpha = 0.6 /length(unique(figure4$covariates))) +
  geom_rect(aes(xmin = TATE_ci_lb, xmax = TATE_ci_ub, ymin = 0,
                ymax = length(unique(figure4$covariates))+1),
            fill = pal_lancet("lanonc")(9)[9], color = NA, alpha = 0.3 / length(unique(figure4$covariates))) +
  geom_vline(aes(xintercept = mean_TATE), color = pal_lancet("lanonc")(9)[9], size = 0.75) +
  geom_vline(aes(xintercept = mean_SATE), color = pal_lancet("lanonc")(9)[8], size = 0.75) +
  geom_hline(aes(yintercept = 0, color = "S = 1"), size = 0.01) +  # for legend only
  geom_hline(aes(yintercept = 0, color = "S = 0"), size = 0.01) +  # for legend only
  geom_point(aes(y = covariates, x = mean_estimate, fill = sufficient), size = 3, shape = 21) +
  geom_errorbar(aes(y=covariates, x = mean_estimate, xmax = estimate_ci_ub, xmin = estimate_ci_lb),
                size = 0.75, width = 0.7/length(unique(figure4$covariates)))+
  theme_minimal() + labs(y = "Covariate set", x = "Risk difference")+
  scale_x_continuous(breaks = c(0.3, 0.4, 0.5)) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(name = "Sample:",
                     values = c(pal_lancet("lanonc")(9)[9],
                                pal_lancet("lanonc")(9)[8])) +
  scale_fill_manual(name = "Sufficiency of covariate set:",
                    values = c(
                      "Sufficient"     = pal_lancet("lanonc")(9)[3],
                      "S-dependent"    = pal_lancet("lanonc")(9)[6],
                      "Insufficient" = pal_lancet("lanonc")(9)[7]
                    )) +
  guides(color = guide_legend(order = 2, override.aes = list(linewidth = 0.75)),
         fill = guide_legend(order = 1)
  ) +
  facet_grid2(cols = vars("",case2),
              space = "free_y") +
  theme(axis.ticks.x = element_line(), panel.border = element_rect(color = "black", fill = NA),
        legend.position = "bottom",legend.justification = c(0.5,0.5),legend.title.align = 0.5,
        legend.title.position = "left",
        legend.margin = margin(t = -10, unit = "pt"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()
  )
figure4_plot


## Figure 5
# prepare plot data

figure5 = bind_rows(df_output_means, df_output_means2)
figure5 = mutate(figure5, soca2 = case_when(round(soca,2) == 1 ~ "0.3*DM1-0.5*S",
                                            round(soca,2) == round(2/3,2) ~ "0.2*DM1-0.33*S",
                                            round(soca,2) == round(1/3,2) ~ "0.1*DM1-0.17*S"))

figure5 = mutate(figure5, soi2 = case_when(round(soi,2) == 1 ~ "0.3*Z*DM1",
                                           round(soi,2) == round(2/3,2) ~ "0.2*Z*DM1",
                                           round(soi,2) == round(1/3,2) ~ "0.1*Z*DM1"))

figure5$soca = factor(figure5$soca2, levels = c("0.1*DM1-0.17*S","0.2*DM1-0.33*S","0.3*DM1-0.5*S"))
figure5$soi = factor(figure5$soi2, levels = c("0.1*Z*DM1","0.2*Z*DM1","0.3*Z*DM1"))
figure5$n2 = comma(figure5$n)
figure5$n2 = factor(figure5$n2,
                    levels = c("2,000","20,000","100,000"))
figure5 = filter(figure5, case == "S-->EMM", covariates == "DM2, P")

# plot
figure5_plot = ggplot(figure5) +
  geom_rect(aes(xmin = SATE_ci_lb, xmax = SATE_ci_ub, ymin = 0,
                ymax = length(unique(figure5$soi))+1),
            fill = pal_lancet("lanonc")(9)[8], color = NA, alpha = 0.4/length(unique(figure5$soi))) +
  geom_rect(aes(xmin = TATE_ci_lb, xmax = TATE_ci_ub, ymin = 0,
                ymax = length(unique(figure5$soi))+1),
            fill = pal_lancet("lanonc")(9)[9], color = NA, alpha = 0.2/length(unique(figure5$soi))) +
  geom_vline(aes(xintercept = mean_TATE), color = pal_lancet("lanonc")(9)[9], linewidth = 0.75) +
  geom_vline(aes(xintercept = mean_SATE), color = pal_lancet("lanonc")(9)[8],linewidth = 0.75) +
  geom_hline(aes(yintercept = 0, color = "S = 1"), linewidth = 0.01) +  # for legend only
  geom_hline(aes(yintercept = 0, color = "S = 0"), linewidth = 0.01) +  # for legend only
  geom_point(aes(y = soca, x = mean_estimate, fill = "Transported"),
             shape = 21,
             size = 3) +
  geom_errorbar(aes(y=soca2,x = mean_estimate, xmax=estimate_ci_ub, xmin=estimate_ci_lb),
                size = 0.75, width = 0.15)+
  theme_minimal() + labs(y = "Strength of collider association (P = 0.5 + ...)", x = "Risk difference")+
  scale_y_discrete(limits = rev) +
  scale_color_manual(name = "ATE:",
                     values = c(pal_lancet("lanonc")(9)[9],
                                pal_lancet("lanonc")(9)[8])) +
  scale_fill_manual(name = NULL,
                    values = c(
                      pal_lancet("lanonc")(9)[6])) +
  guides(color = guide_legend(order = 1, override.aes = list(linewidth = 0.75)),
         fill = guide_legend(order = 2)) +
  facet_grid2(rows = vars("Strength of interaction",soi2),
              cols = vars("Sample size",n2),
              strip = strip_nested(),
              space = "free_y") +
  theme(axis.ticks.x = element_line(), panel.border = element_rect(color = "black", fill = NA),
        legend.position = "bottom",legend.justification = c(0.5,0.5),legend.title.align = 0.5,
        legend.title.position = "left",
        legend.margin = margin(t = -10, r = -5,
                               unit = "pt"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()
        
  )
figure5_plot
