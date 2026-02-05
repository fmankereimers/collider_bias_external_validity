##############################
### required packages
library(dplyr)

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
sample_sizes = c(2000,20000,100000)

#strength of collider association
soca = c("Strong" = 1, "Medium" = 2/3, "Weak" = 1/3)

#strength of interaction
soi = c("Strong" = 1, "Medium" = 2/3, "Weak" = 1/3)

output = list()
for (n in sample_sizes){
  for (i in soca){
    for (m in soi){
      for (j in 1:length(covariates)){
        for (t in 1:bootstraps){
          seed = 2765+t*3
          set.seed(seed)
          S = rbinom(n,1,.5)
          A = rbinom(n,1,.5)
          DM1 = rbinom(n,1,.5)
          P = rbinom(n,1,0.55 + DM1*0.3*i - S*0.5*i)
          DM2 = rbinom(n,1,0.2 + S*0.4)

          Y = rbinom(n,1,prob = 0.05 + A*0.1 + DM1*0.1 + DM2*0.1 +
                       A*DM1*0.3*m + A*DM2*0.3)
          sim_data = cbind(S,Y,A,DM1,DM2,P)
          source = sim_data[sim_data[,"S"] == 1,]
          SATE = mean(source[source[,"A"]==1,"Y"]) - mean(source[source[,"A"]==0,"Y"])
          target = sim_data[sim_data[,"S"] == 0,]
          TATE = mean(target[target[,"A"]==1,"Y"]) - mean(target[target[,"A"]==0,"Y"])
          sim_data = as.data.frame(sim_data)
          data = source[sample(nrow(source), size = nrow(source), replace = TRUE), ]
          data = data.frame(rbind(data, target))
          data = select(data, "S","Y","A",paste(covariates[[j]], sep = ","))
          S1data_A1<-filter(data, S==1 & A==1)
          #S1data_A1<-select(S1data_A1, -S,-A)
          covariates_plus = paste(covariates[[j]], collapse = "+")
          OM1mod<-glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A1, family = "binomial")
          p1<- predict(OM1mod,newdata=data, type="response")
          S1data_A0<-filter(data, S==1 & A==0)
          #S1data_A0<-select(S1data_A0, -S,-A)
          OM0mod<-glm(formula=paste("Y ~ ", covariates_plus ,sep = ""), data=S1data_A0, family = "binomial")
          p0<- predict(OM0mod,newdata=data, type="response")
          data$p1<-p1
          data$p0<-p0
          S0sub<-subset(data, S==0)
          OM_1<-mean(S0sub$p1)
          OM_0<-mean(S0sub$p0)
          OM<-mean(S0sub$p1)-mean(S0sub$p0)
          outtable=data.frame(covariates = paste(covariates[[j]], collapse = ", "), estimate = OM, SATE, TATE, n, soca = i, soi = m, seed)
          output=append(output,list(outtable))
        }
      }
    }
  }
}

df_output = do.call(rbind, output)
nrow(df_output)
tail(df_output)

df_output_means = df_output %>%
  group_by(covariates, n, soca , soi) %>%
  summarise(across(c(SATE, TATE, estimate), \(x) mean(x, na.rm = TRUE), .names = "mean_{.col}"),
            across(c(SATE, TATE, estimate), \(x) sd(x, na.rm = TRUE), .names = "{.col}_se"))

df_output_means = mutate(df_output_means,
                         estimate_ci_lb = mean_estimate - 1.96 * estimate_se,
                         estimate_ci_ub = mean_estimate + 1.96 * estimate_se,
                         SATE_ci_lb = mean_SATE - 1.96 * SATE_se,
                         SATE_ci_ub = mean_SATE + 1.96 * SATE_se,
                         TATE_ci_lb = mean_TATE - 1.96 * TATE_se,
                         TATE_ci_ub = mean_TATE + 1.96 * TATE_se)

df_output_means = mutate(df_output_means, sufficient  = case_when(covariates %in% c("DM2","DM2, P") ~ "S-dependent",
                                                                  covariates %in% c("DM1","DM1, P", "P") ~ "Non-sufficient",
                                                                  TRUE ~ "Sufficient"))
