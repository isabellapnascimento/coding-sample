# Coding Sample 
# Author: Isabella Nascimento
# Purpose:

#   This script is adapted from an econometrics course where I served as a TA, 
#   rewritten in a pedagogical format for easy replication and learning of basic
#   econometric models. :)

#   It covers three standard tasks:
#     (Section 1) Binary outcome: Linear Probability Model (LPM) vs Probit
#     (Section 2) Difference-in-Differences (DiD)
#     (Section 3) Two-way Fixed Effects (firm & year)
#
# How to use?
#   1) Adjust the working directory below.
#   2) Ensure data files are present in your computer:
#        - smoking.dta  
#        - jtrain.dta  
#      The DiD section uses organ_donations from the 'causaldata' package.
#   3) Run in RStudio. Output is printed to the console.

set.seed(41400)           # reproducibility!
options(scipen = 999)     # cleaner numeric printing

# Load libraries
library(haven)      
library(mfx)        
library(tidyverse)  
library(fixest)     
library(ggplot2)    
library(causaldata) 
library(sandwich)   
library(broom)      
library(dplyr)      
library(lmtest)     
library(car)       

# Adjust this path to where your data files live.
setwd("/Users/isabellanascimento/Desktop/U Chicago/Reg Analysis")

# load datasets
smoking_raw <- read_dta("smoking.dta")
jtrain_raw  <- read_dta("jtrain.dta")

# =====================================================================
# Section 1 — Binary Dependent Variable: Smoking bans and smoking probability
# =====================================================================

# What we are doing in Section 1:
# - Goal: Estimate how smoking bans (smkban) relate to the probability of smoking (smoker).
# - Approach: Start with a simple Linear Probability Model (LPM), add controls, then compare
#   to a nonlinear Probit model. Report robust (HC1) standard errors and compute marginal
#   effects (Average Partial Effect, APE; and Partial Effect at the Average, PEA).

# (A) LPM: smoker ~ smkban
#   - Simple difference-in-means interpretation in an LPM framework.
mod1 <- lm(smoker ~ smkban, data = smoking_raw)
print(coeftest(mod1, vcov. = vcovHC(mod1, type = "HC1")))

# Interpretation: 
# All else equal, having a smoking ban decreases the probability of 
# someone smoking by 7.77 percentage points.

# (B) LPM with controls
#   - Adds demographics and education to adjust for confounding.
mod3 <- lm(
  smoker ~ smkban + female + age + I(age^2) +
    hsdrop + hsgrad + colsome + colgrad + black + hispanic,
  data = smoking_raw
)
print(coeftest(mod3, vcov. = vcovHC(mod3, type = "HC1")))

# Interpretation:
# After controlling for demographic and socioeconomic factors, the presence of a smoking ban
# decreases the probability of smoking by about 4.7 percentage points.
# Accounting for these covariates explains part of the variation in smoking behavior,
# reducing the unadjusted impact seen in model (A).
# In particular, gender and age are strongly correlated with smoking habits and likely
# overlap with some of the variation that was previously attributed to the ban.

# (C) Predicted probability range (LPM can go outside [0,1])
rng_lpm <- range(predict(mod3), na.rm = TRUE)
print(rng_lpm)

# Interpretation:
# The predicted values range outside the 0–1 interval (e.g., negative probabilities),
# which is a known limitation of the Linear Probability Model (LPM). :(
# This motivates using nonlinear probability models like Probit or Logit.

# (D) Joint F-test for the education dummies (HC1 robust)
ftest <- linearHypothesis(
  mod3,
  c("hsdrop = 0", "hsgrad = 0", "colsome = 0", "colgrad = 0"),
  vcov. = vcovHC(mod3, type = "HC1")
)
print(ftest)

# Interpretation:
# The joint F-test strongly rejects the null hypothesis that education has no effect
# on smoking probability. Lower levels of education (e.g., high school dropout)
# are associated with higher smoking rates, suggesting education plays
# a significant role in shaping health behaviors.


# (E) Probit with the same controls + robust SE; predicted probabilities always in [0,1]
probit <- glm(
  smoker ~ smkban + female + age + I(age^2) +
    hsdrop + hsgrad + colsome + colgrad + black + hispanic,
  family = binomial(link = "probit"),
  data = smoking_raw
)
print(coeftest(probit, vcov. = vcovHC(probit, type = "HC1")))
rng_probit <- range(predict(probit, type = "response"), na.rm = TRUE)
print(rng_probit)


# Interpretation:
# The Probit model confirms the same directional effects as the LPM.
# Its coefficients differ in scale because they describe the effect on
# the latent index of the probability rather than the probability itself.
# Predicted probabilities now fall cleanly between 0 and 1, resolving
# the logical limitation of the LPM.

# Note!
#   - LPM and Probit often agree on signs and significance!
#   - Probit coefficients are on the latent index scale; to compare magnitudes with LPM, we need to 
#     use marginal effects (APE/PEA), not the raw probit coefficient.

# (F) Average Partial Effect (APE) for smkban — via package and manual formula
ape_pkg <- probitmfx(
  smoker ~ smkban + female + age + I(age^2) +
    hsdrop + hsgrad + colsome + colgrad + black + hispanic,
  data = smoking_raw, atmean = FALSE
)
print(ape_pkg$mfxest["smkban", , drop = FALSE])

xb   <- predict(probit, type = "link")
phi  <- dnorm(xb)
beta <- coef(probit)[["smkban"]]
ape_manual <- mean(phi * beta, na.rm = TRUE)
print(ape_manual)

# Interpretation:
# The average partial effect (APE) of smkban in the Probit model is close to -0.047,
# nearly identical to the LPM coefficient. This consistency reinforces that both
# models capture the same underlying marginal effect: a small but statistically
# meaningful reduction in smoking probability due to smoking bans.


# (G) Partial Effect at the Average (PEA) for smkban — using covariate means
xbar <- smoking_raw %>%
  summarise(
    across(c(smkban, female, age), ~ mean(.x, na.rm = TRUE)),
    `I(age^2)` = mean(age^2, na.rm = TRUE),
    across(c(hsdrop, hsgrad, colsome, colgrad, black, hispanic), ~ mean(.x, na.rm = TRUE))
  ) %>%
  as.list()

b    <- coef(probit)
xvec <- c(
  1, xbar$smkban, xbar$female, xbar$age, xbar$`I(age^2)`,
  xbar$hsdrop, xbar$hsgrad, xbar$colsome, xbar$colgrad, xbar$black, xbar$hispanic
)
pea  <- dnorm(sum(b * xvec)) * b["smkban"]
print(pea)

# Interpretation:
# The Partial Effect at the Average (PEA) is nearly identical to the APE,
# again around -0.048. Both metrics confirm that smoking bans are associated
# with a roughly 5-percentage-point decline in smoking probability when other
# individual characteristics are held constant.

#                                 Wrap-up for Section 1:
#   - The LPM’s smkban effect is typically close to the Probit APE/PEA (comparable scales).
#   - Prefer Probit (or Logit) when probabilities outside [0,1] are a concern,
#     but LPM is a useful baseline in percentage-point units :)

# =====================================================================
# Section 2 — Difference-in-Differences: California organ-donation policy
# =====================================================================

# What we are doing in Section 2:
# - Context: California introduced a policy around 2011 affecting organ-donation registrations.
# - Goal: Estimate policy effect using Difference-in-Differences (DiD): compare California to
#   other states before vs after the change, and examine pre-treatment trends.

df <- organ_donations #our data!

# (A) Pretrend slice 
#   - We look at 2010Q4 to 2012Q1 to eyeball parallel trends pre-policy.
pre <- df %>%
  filter(Quarter_Num %in% 1:6) %>%
  mutate(Group = if_else(State == "CA", "California", "Other states")) %>%
  group_by(Group, Quarter_Num) %>%
  summarise(Avg_Rate = mean(Rate, na.rm = TRUE), .groups = "drop")
print(head(pre))

# Interpretation:
# Pre-treatment averages suggest broadly parallel trends between California and other states
# through early 2012. While not definitive, this descriptive evidence supports the DiD assumption
# that, absent the policy, California’s trend would have continued similarly to the control group.


# (B) DiD regression with clustered SEs by State
#   - post: indicator for post-policy period
#   - treated: California vs other states
#   - interact: the DiD term; the policy effect for CA relative to control
did_df <- df %>%
  mutate(
    post     = as.integer(Quarter_Num >= 4),
    treated  = as.integer(State == "CA"),
    interact = post * treated
  )

did_model <- feols(Rate ~ post + treated + interact, data = did_df, cluster = "State")
print(summary(did_model))

# Interpretation:
# - Intercept: Baseline organ-donation rate for control states in the pre-policy period.
# - post: Average change in organ-donation rates for control states after the policy; captures
#   nationwide or time-related shifts.
# - treated: Pre-policy difference between California and other states; measures baseline gap.
# - interact (DiD estimate): The policy’s estimated effect on California’s organ-donation rate,
#   relative to other states.

# The coefficient on `interact` is negative (around -0.022), indicating that, contrary to expectations,
# organ-donation rates in California fell by roughly 2 percentage points relative to other states after
# the policy. This unexpected result may reflect short-term implementation issues, data lags,
# or behavioral responses not captured by the model.

#                                 Wrap-up for Section 2:
#   - The Difference-in-Differences (DiD) model estimates the effect of California’s policy
#     by comparing changes in organ-donation rates in California to those in other states over time.
#   - While we cannot *prove* that the parallel trends assumption holds, examining pre-policy trends
#     helps assess whether treatment and control groups evolved similarly before the intervention.
#   - The negative and significant interaction term indicates that, after the policy, California’s
#     donation rate declined relative to other states—an unexpected finding that may reflect
#     short-term implementation effects or other unobserved factors.

# =====================================================================
# Section 3 — Two-way Fixed Effects: Do grants raise training hours in firms?
# =====================================================================

# What we are doing in Section 3:
# - Data: Firm-level panel with training grants and hours of training per employee (hrsemp).
# - Model: Two-way Fixed Effects (firm & year). This identifies within-firm changes over time,
#   controlling for common shocks each year and time-invariant firm heterogeneity.

fe_model <- feols(
  hrsemp ~ grant + grant_1 + lemploy | fcode + year,
  data = jtrain_raw, cluster = "fcode"
)
print(summary(fe_model))

# =====================================================================
# Section 3 — Two-way Fixed Effects: Do grants raise training hours in firms?
# =====================================================================

# This section uses firm-level panel data to estimate how receiving a training grant
# affects the average number of training hours per employee over time.
# A two-way Fixed Effects model (firm and year) controls for both firm-specific
# characteristics that do not change over time and year-specific shocks
# affecting all firms.

fe_model <- feols(
  hrsemp ~ grant + grant_1 + lemploy | fcode + year,
  data = jtrain_raw, cluster = "fcode"
)
print(summary(fe_model))

# Interpretation:
# - grant: Receiving a grant increases average training hours per employee by about 34 hours,
#   a statistically significant and economically meaningful effect.
# - grant_1: The lagged effect is small and statistically insignificant,
#   suggesting limited persistence once the grant period ends.
# - lemploy: The coefficient on firm size (log employees) is negative and insignificant,
#   which is expected given that within-firm changes in employment are small over time.

#                                 Wrap-up for Section 3:
#   - The two-way Fixed Effects model isolates within-firm variation by accounting for
#     both firm and year unobserved characteristics.
#   - Training grants have a strong contemporaneous effect on training intensity but
#     no lasting impact after the funding period.
#   - This approach demonstrates how Fixed Effects models capture causal relationships
#     in panel data by controlling for time-invariant differences across firms.
