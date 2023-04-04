library(tidyverse)
library(survey)
library(corrr)
library(readr)
library(ROCit)
library(broom)
library(writexl)

wvs_df <- read_csv("WVS_Cross-National_Wave_7_csv_v5_0.csv")

# Selecionar as variáveis a serem utilizadas
# Remover as observações com valores faltantes ou perguntas não respondidas
wvs <- wvs_df %>% 
  as_tibble() %>% 
  select(W_WEIGHT, country = B_COUNTRY_ALPHA, A_STUDY, A_WAVE, S018, S025,
         happy = Q46, gender = Q260, age = Q262, age6 = X003R, marital = Q273,
         children = Q274, educ = Q275, educ3 = Q275R, job = Q279, chief = Q285,
         mpol_leader = Q29, conf_wmov = Q80, member_wgroup = Q104, w_lcorrupt = Q119,
         religious = Q173, prostitution = Q183, abortion = Q184, divorce = Q185,
         m_beatw = Q189, interest_pol = Q199, discuss_pol = Q200, w_poleq = Q233,
         wsame_rights = Q249, incomel = Q288R) %>% 
  filter(if_all(where(is.numeric), ~ . >= 0),
         job != 8)

# Ajuste das variáveis e criação do índice identidade
wvs <- wvs %>% 
  mutate(across(c(wsame_rights, prostitution, abortion, divorce, m_beatw), ~ . / 10), # normalizando indices de 0 a 10
         mpol_leader = ifelse(mpol_leader <= 2, 0, 1),
         conf_wmov = ifelse(conf_wmov <= 2, 1, 0),
         w_lcorrupt = ifelse(w_lcorrupt <= 2, 1, 0),
         m_beatw = 1 - m_beatw) %>% 
  mutate(identidade = mpol_leader + conf_wmov + w_lcorrupt + wsame_rights +
                      prostitution + abortion + divorce + m_beatw) 

# Ajustes de tipos e calculo do multiplicador do peso amostral
wvs <- wvs %>% 
  mutate(across(c(country, gender, age6, marital, educ3, job, chief,
                  mpol_leader, conf_wmov, member_wgroup, w_lcorrupt,
                  religious, incomel), as.factor)) %>% 
  group_by(country) %>%
  mutate(S018R = 1000 / n()) %>%
  mutate(age2 = age^2,
         happy = ifelse(happy <=2, 1, 0),
         gender = fct_recode(gender,
                          "Male" = "1",
                          "Female" = "2"),
         # gender = fct_relevel(gender,
         #                      c("Female", "male")),
         marital = fct_collapse(marital, 
                                "Married" = c("1", "2"),
                                "Divorced" = "3",
                                "Separated" = "4",
                                "Widowed" = "5",
                                "Single" = "6"),
         job = fct_collapse(job,
                            "Full-time" = c("1", "3"),
                            "Not working" = c("2", "4", "6"),
                            "Housewife" = "5",
                            "Unemployed" = "7"),
         chief = fct_recode(chief,
                            "0" = "2",
                            "1" = "1"),
         member_wgroup = fct_collapse(member_wgroup,
                                      "0" = "0",
                                      "1" = c("1", "2")),
         religious = fct_collapse(religious,
                                  "Religious" = "1",
                                  "Not religious" = c("2", "3")))

# Primeiro estagio 
# Estimando a identidade
fit_idd <- lm(identidade ~ age + age2 + gender + educ + children + marital + 
                           chief + job + member_wgroup + country,
              data = wvs)
summary(fit_idd)

svy <- svydesign(id = ~1,
                 strata = ~S025,
                 weights = ~(W_WEIGHT * S018R), 
                 data = wvs)

svy$variables <- svy$variables %>% 
  mutate(identidadef = fit_idd$fitted.values) %>% 
  group_by(country) %>% 
  mutate(P = mean(identidade)) %>% 
  mutate(delta = (identidadef - P)^2)


fit_delta <- svyglm(delta ~ age + age2 + gender + children + marital + religious + job,
                    design = svy,
                    family = stats::gaussian(link = "identity"))
summary(fit_delta)

# MPL
fit_mpl <- svyglm(happy ~ age + age2 + gender + children + marital + religious + educ +
                    delta + delta:gender,
                  design = svy,
                  family = stats::gaussian(link = "identity"))
summary(fit_mpl)

# Probit
fit_probit <- svyglm(happy ~ age + age2 + gender + children + marital + religious + educ +
                       delta + delta:gender,
                     design = svy,
                     family = binomial(link = "probit"))
summary(fit_probit)

# Logit
fit_logit <- svyglm(happy ~ age + age2 + gender + children + marital + religious + educ +
                      delta + delta:gender,
                    design = svy,
                    family = binomial(link = "logit"))
summary(fit_logit)

exp(cbind(coef(fit_logit)))


# Calculando diagnósticos do modelo
roc <- rocit(score = fit_logit$fitted.values, class = fit_logit$y)
plot(roc)

# Estimando um segundo logit sem aplicar o quadrado no delta
svy2 <- svydesign(id = ~1,
                  strata = ~S025,
                  weights = ~(W_WEIGHT * S018R), 
                  data = wvs)

svy2$variables <- svy2$variables %>% 
  mutate(identidadef = idd_fit$fitted.values) %>% 
  group_by(country) %>% 
  mutate(P = mean(identidade)) %>% 
  mutate(delta = identidadef - P)

fit_logit2 <- svyglm(happy ~ age + age2 + gender + children + marital + religious + educ +
                      delta + delta:gender,
                    design = svy2,
                    family = binomial(link = "logit"))
summary(fit_logit2)

# Comparando os modelos
roc2 <- rocit(score = fit_logit2$fitted.values, class = fit_logit2$y)
plot(roc2)

roc2$AUC
roc$AUC

# Analise descritiva
stats_mean <- svymean(~ age + gender + educ + children + marital + 
                        chief + job + member_wgroup,
                      svy,
                      na.rm = TRUE,
                      parm=NA)

# Descritiva identidade
idd_mean <- data.frame(svymean(~ mpol_leader + conf_wmov + w_lcorrupt + 
                      wsame_rights + prostitution + abortion + 
                      divorce + m_beatw,
                    svy,
                    na.rm = TRUE,
                    parm=NA))


# Salvando as tabelas requisitadas como xlsx
write_xlsx(tidy(idd_fit), "fit_idd.xlsx")
write_xlsx(tidy(delta_fit), "fit_delta.xlsx")
write_xlsx(tidy(fit_mpl), "fit_mpl.xlsx")
write_xlsx(tidy(fit_probit), "fit_probit.xlsx")
write_xlsx(tidy(fit_logit), "fit_logit.xlsx")
write_xlsx(stats_mean, "stats.xlsx")
write_xlsx(idd_mean, "idd.xlsx")










