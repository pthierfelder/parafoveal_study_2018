# Load packages
library(broom.mixed)
library(dplyr)
library(huxtable)
library(lme4)
library(optimx)

Data <- Study_1_1_Data

# Pre-processing ----------------------------------------------------------
## Set high and low fixation thresholds for outlier filtration
Min_FFD <- 100
Max_GD <- 1500
Max_TRT <- 2000

## Remove missing observations 
Missing_Observations <- 1 - length(which(Data$FFD_ms > 0))/length(Data$FFD_ms)
Data <- filter(Data, FFD_ms > 0) # Remove missing observations

## Set unrelated/control condition as reference
Data$Preview <- factor(Data$Preview, levels = c("UNR", "COR", "SIGN"))
Data$Sign <- factor(Data$Sign, levels = c("Unrelated", "HLM", "LM", "HM", "HL", "Identical"))

## Filter and report percentages of data filtered
Filtered_Data_Hearing <- filter(Data, FFD_ms > Min_FFD & GD_ms < Max_GD & Pre_Target_Skip == 0)

Hearing_Pre_Tar_Skips_Filtered <- 1 - length(which(Data$FFD_ms > 0 & Data$Pre_Target_Skip == 0))/length(which(Data$FFD_ms > 0))

Hearing_Filtered <- 1 - length(which(Data$Pre_Target_Skip == 0 & Data$GD_ms < Max_GD & Data$FFD_ms > Min_FFD))/length(which(Data$FFD_ms > 0 & Data$Pre_Target_Skip == 0))

cat("\n", round(Missing_Observations * 100), "percent of FFD/GD trials are missed observations.\n",
    "\n", round(Hearing_Pre_Tar_Skips_Filtered * 100), "percent of the FFD/GD data had pre-target skips.\n",
    "\n", round(Hearing_Filtered * 100), "percent of the FFD/GD data was outside of min/max thresholds.\n")

Filtered_Data_Hearing_TRT <- filter(Filtered_Data_Hearing, TRT_ms < Max_TRT)

Hearing_Filtered_TRT <- 1 - length(which(Filtered_Data_Hearing$TRT_ms < Max_TRT))/length(which(Filtered_Data_Hearing$TRT_ms > 0))

cat("\n", round(Hearing_Filtered_TRT * 100), "percent of the TRT data was outside of threshold.\n")

# Models ------------------------------------------------------------------
## Build first set of models with all sign phonological previews as a single variable
## In some cases, we needed to use different optimizers to resolve convergence failures

Hearing_FFD_Preview_GLMM <- glmer(FFD_ms ~ Preview +  (1 | Participant) + (1 | Frame), data = Filtered_Data_Hearing, Gamma(link = "identity"),
                                  glmerControl(optimizer = "bobyqa"))

Hearing_GD_Preview_GLMM <- glmer(GD_ms ~ Preview + (1 | Participant) + (1 | Frame), data = Filtered_Data_Hearing, Gamma(link = "identity"),
                                 glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")))

Hearing_TRT_Preview_GLMM <- glmer(TRT_ms ~ Preview +  (1 | Participant) + (1 | Frame), data = Filtered_Data_Hearing_TRT, Gamma(link = "identity"),
                                  glmerControl(optimizer = "bobyqa"))

## Build second set of models with separate HLM/HL/HM/LM sign phonological conditions
Hearing_FFD_Sign_GLMM <- glmer(FFD_ms ~ Sign +  (1 | Participant) + (1 | Frame), data = Filtered_Data_Hearing, Gamma(link = "identity"),
                            glmerControl(optimizer = "bobyqa"))
summary(Hearing_FFD_Sign_GLMM)

Hearing_GD_Sign_GLMM <- glmer(GD_ms ~ Sign + (1 | Participant) + (1 | Frame), data = Filtered_Data_Hearing, Gamma(link = "identity"),
                                                               glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")))
                                                   
summary(Hearing_GD_Sign_GLMM)

Hearing_TRT_Sign_GLMM <- glmer(TRT_ms ~ Sign +  (1 | Participant) + (1 | Frame), data = Filtered_Data_Hearing_TRT, Gamma(link = "identity"),
                            glmerControl(optimizer = "Nelder_Mead"))

summary(Hearing_TRT_Sign_GLMM)

# Tables ------------------------------------------------------------------
## Table for first set of models
Study_1_1_Hearing_M1 <- huxreg("FFD" = Hearing_FFD_Preview_GLMM, "GD" = Hearing_GD_Preview_GLMM, "TRT" = Hearing_TRT_Preview_GLMM,
                            coefs= c("Unrelated" = "(Intercept)", 
                                     "Sign Phonological" = "PreviewSIGN", 
                                     "Identical" = "PreviewCOR"), 
                            statistics = c("# of Observations" = "nobs", "logLik", "AIC"),
                            number_format = 2,
                            bold_signif = 0.05,
                            error_format = "(SE={std.error}, t={statistic})",
                            stars = NULL)
print(Study_1_1_Hearing_M1)
huxtable::quick_docx(Study_1_1_Hearing_M1, file = 'Study_1_1_Hearing_M1.docx')

## Table for second set of models
Study_1_1_Hearing_M2 <- huxreg("FFD" = Hearing_FFD_Sign_GLMM, "GD" = Hearing_GD_Sign_GLMM, "TRT" = Hearing_TRT_Sign_GLMM,
                            coefs= c("Unrelated" = "(Intercept)", 
                                     "HLM" = "SignHLM", 
                                     "HL" = "SignHL", 
                                     "HM" = "SignHM", 
                                     "LM" = "SignLM", 
                                     "Identical" = "SignIdentical"), 
                            statistics = c("# of Observations" = "nobs", "logLik", "AIC"),
                            number_format = 2,
                            bold_signif = 0.05,
                            error_format = "(SE={std.error}, t={statistic})",
                            stars = NULL)

print(Study_1_1_Hearing_M2)
huxtable::quick_docx(Study_1_1_Hearing_M2, file = 'Study_1_1_Hearing_M2.docx')
