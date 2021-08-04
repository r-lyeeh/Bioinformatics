# installation
# install.packages("devtools") #!important
devtools::install_github("PHP2560-Statistical-Programming-R/r-framingham") 
# load library
library(frisk)
sample_size <- 100
# simulate patients data
df <- data.frame(age=sample(30:70,sample_size,rep=TRUE),
                 gender=sample(c("M","F"),100,rep=TRUE),
                 bmi=sample(16:48, sample_size, rep = TRUE),
                 hdl=sample(10:100,sample_size,rep=TRUE),
                 chl=sample(100:400,sample_size,rep=TRUE),
                 sbp=sample(90:200,sample_size,rep=TRUE),
                 isSbpTreated=sample(c(TRUE,FALSE),sample_size,rep=TRUE),
                 smoking=sample(c(TRUE,FALSE),sample_size,rep=TRUE),
                 diabetes=sample(c(TRUE,FALSE),sample_size,rep=TRUE)
)

# call frisk function case no bmi
calc_card_10(df, age="age", gender="gender", cholesterol="chl", 
             hdl="hdl", sbp="sbp", is_sbp_under_treatment="isSbpTreated",
             smoking_status="smoking", diabetes_status="diabetes"
)

# Install from CRAN
install.packages("CVRisk")

# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("vcastro/CVRisk")
library(CVrisk)
ascvd_10y_accaha(race = "aa", gender = "male", age = 55, 
                 totchol = 213, hdl = 50, sbp = 140, 
                 bp_med = FALSE, smoker=0, diabetes=0)
library(dplyr)

sample_data %>% 
  compute_CVrisk(age = "age", race = "race", gender = "gender", bmi = "BMI", 
                 sbp = "sbp", hdl = "hdl", totchol = "totchol", bp_med = "bp_med", 
                 smoker = "smoker", diabetes = "diabetes")

library(framinghamRiskEquation)
packageVersion("framinghamRiskEquation")

framingham_riskequation(
  data.frame(InternalID = 1, BP = "135/80", Sex = "Female",
             Age = 55, SmokingStatus = "Smoker",
             CholHDLRatio = 230/48, Diabetes = TRUE, LVH = FALSE,
             CardiovascularDisease = FALSE, PersistentProteinuria = FALSE,
             eGFRValue = NA, eGFRUnits = NA,
             UrineAlbuminValue = NA, UrineAlbuminUnits = NA,
             FamilialHypercholesterolaemia = NA, Cholesterol = 5.96,
             Ethnicity = NA),
  outcome = "CHD", years = 10
)

framingham_riskequation(
  data.frame(InternalID = 1, BP = "135/80", Sex = "Female",
             Age = 55, SmokingStatus = "Smoker",
             CholHDLRatio = 230/48, Diabetes = TRUE, LVH = FALSE,
             CardiovascularDisease = FALSE, PersistentProteinuria = FALSE,
             eGFRValue = NA, eGFRUnits = NA,
             UrineAlbuminValue = NA, UrineAlbuminUnits = NA,
             FamilialHypercholesterolaemia = NA, Cholesterol = 5.96,
             Ethnicity = NA),
  outcome = "CVD", years = 5
)
