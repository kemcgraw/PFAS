setwd("D://KM//Columbia//Project VIVA//data")
library(tidyverse)
library(sas7bdat)
library(haven)
library(tableone)
library(broom)
library(labelled)
library(ggthemes)
library(ggcorrplot)
library(Hmisc)
library(sandwich)
library(lme4)
library(nlme)
library(gridExtra)
library(lattice)
library(gfoRmula)
require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
library(stats)
library(foreach)
library(doParallel)


##Read Data---------------------------------------------------------------------
pfas <- read_sas('katlyn_mcgraw_040124.sas7bdat') %>% 
  mutate(age1c = age_mom_enroll_d,                #age mom enrolled
         age2c = mom_age_days_blood_mbl17/365.25,        #age mom at midlife blood draw
         edu = education_mom_epi_epia_d,          #education mom enrolled
         race = case_when(                        #race mom enrolled
           race2_mom_epi_epia_d == 'white' ~ 0,             #white               
           race2_mom_epi_epia_d == 'black' ~ 1,             #black
           race2_mom_epi_epia_d == 'hispa' ~ 2,             #hispanic
           race2_mom_epi_epia_d == 'asian' ~ 3,             #asian
           race2_mom_epi_epia_d == 'amind' ~ 4,             #american indian
           race2_mom_epi_epia_d == 'more than 1 race' ~ 5,  #more than 1 race
           race2_mom_epi_epia_d == 'other' ~ 5,             #other   
           race2_mom_epi_epia_d == '' ~ 5),                 #no response???????
         maritalstat1c = case_when(               #marital status mom enrolled
           marital_epi_epia_d == 1 ~ 1,             #married               
           marital_epi_epia_d == 2 ~ 2,             #not married living with partner
           marital_epi_epia_d == 3 ~ 3,             #never married
           marital_epi_epia_d == 4 ~ 4,             #divorce
           marital_epi_epia_d == 5 ~ 5,             #separated
           marital_epi_epia_d == 6 ~ 6),            #widowed
         maritalstat2c = case_when(               #marital status mom midlife
           marital_status_qu17 == 1|          
           marital_status_qu17 == 2 ~ 1,            #married to viva father/#married to not viva father                       
           marital_status_qu17 == 3 ~ 2,            #not married but living with partner
           marital_status_qu17 == 4 ~ 3,            #never married
           marital_status_qu17 == 5 ~ 4,            #divorced
           marital_status_qu17 == 6 ~ 5,            #separated
           marital_status_qu17 == 7 ~ 6),           #widowed
         hhinc1c = case_when(                     #household income mom enrolled
           income_hh_epq_epqa_d == 1|                #<5k
           income_hh_epq_epqa_d == 2|                #5-10k
           income_hh_epq_epqa_d == 3|                #10-20k                       
           income_hh_epq_epqa_d == 4|                #<40k
           income_hh_epq_epqa_d == 5 ~ 0,            #40k-70k
           income_hh_epq_epqa_d == 6 ~ 1,            #>70k
           income_hh_epq_epqa_d == 9 ~ 0),           #don't know / not sure if I should code as 0 or NA
         hhinc2c = case_when(                     #household income mom midlife
           household_income_qu17 == 1|               #<40k          
           household_income_qu17 == 2 ~ 0,           #40-70k
           household_income_qu17 == 3|               #>70k   (70-100k)
           household_income_qu17 == 4|               #100k-150k
           household_income_qu17 == 5 ~ 1,           #>150k
           household_income_qu17 == 9 ~ 0),          #don't know / not sure if I should code as 0 or NA
         highbp = highbp_mom_epi_epia_d,          #history of high BP at early pregnancy interview
         t1d = t1diab_mom_epi_epia_d,             #history of T1DM at early pregnancy interview
         t2d = t2diab_mom_epi_epia_d,             #history of T2DM at early pregnancy interview
         htnpreg = preE_vitals_hosp_d,            #hypertensive disorders in pregnancy
         sercreat = tri1_creatinine_mg_dL,        #trimester 1 serum creatinine
         albumin =  tri1_albumin_g_dL,            #trimester 1 albumin
         egfr = CKD_EPIcr_R,                      #trimester 1 egfr - ckdepi without race equation
         menopause = mom_period_stopped_qu17,     #menopause at midlife
         gestage = ga_wks_BLD1,                   #gestational age at 1st trimester blood draw
         bfdur1c = bfdur_mos_oyq_d,               #breastfeeding duration for the index pregnancy at 1y postpartum
         bfdur2c = lifetime_bf_dur_mos_MT,        #breastfeeding lifetime duration months at midlife
         hei1c = AHEI_FFQ1_n2128,                 #AHEI
         hei2c = hei2015_total_score,             #HEI
         parity1c = case_when(                    #number of pregnancies at early pregnancy interview
           parity_d == 0 ~ 0,                           
           parity_d == 1 ~ 1,   
           parity_d == 2 ~ 2,
           parity_d == 3 ~ 3,
           parity_d == 4 ~ 4,
           parity_d >= 5 ~ 5),
         parity2c = case_when(                    #number of pregnancies at midlife interview
           preg_total_qu17 == 0 ~ 0,                           
           preg_total_qu17 == 1 ~ 1,   
           preg_total_qu17 == 2 ~ 2,
           preg_total_qu17 == 3 ~ 3,
           preg_total_qu17 == 4 ~ 4,
           preg_total_qu17 >= 5 ~ 5),
         #parity1c = case_when(                    #number of pregnancies at early pregnancy interview
           #parity_d == 0 ~ 0,                           
           #parity_d == 1 ~ 1,   
           #parity_d >= 2 ~ 2),
         #parity2c = case_when(                    #number of pregnancies at midlife interview
           #preg_total_qu17 == 0 ~ 0,                           
           #preg_total_qu17 == 1 ~ 1,   
           #preg_total_qu17 >= 2 ~ 2),
         alcohol1c = alc_d_f1,                      #alcohol intake at early pregnancy
         alcohol2c = a_drinks_MT,                   #alcohol intake at midlife
         dash = DASH_f1,                          #dash diet score at early pregnancy
         bmiprepreg = bmi_mom_prepreg_d,          #####not sure if i can use this as time-varying information or not since not technically time one
         bmi_midlf = mom_bmi_mt,
         cig1c = ifelse(smokpreg_final_d == 'xnever', 0, ifelse(smokpreg_final_d == 'former',1,2)),
         cig2c = case_when(
           cigarettes_ever_qu17 == 2 ~ 0,                               ##neversmoker
           cigarettes_ever_qu17 == 1 & cigarettes_12mon_qu17 == 2 ~ 1,  ##formersmoker 
           cigarettes_ever_qu17 == 1 & cigarettes_12mon_qu17 == 1 ~ 2),
         ecig2c = case_when(
           ecigarettes_ever_qu17 == 2 ~ 0,                              ##nevervaper
           ecigarettes_ever_qu17 == 1 & ecigarettes_freq_qu17 == 1 ~ 1, ##formervaper 
           ecigarettes_ever_qu17 == 1 & (ecigarettes_freq_qu17 == 2|ecigarettes_freq_qu17 == 3|ecigarettes_freq_qu17 == 4|
                                           ecigarettes_freq_qu17 == 5|ecigarettes_freq_qu17 == 6) ~ 2),
         mepfosa1c = Me_PFOSA_AcOH2,
         pfhxs1c = PFHxS,
         pfna1c = PFNA2,
         pfda1c = PFDeA,
         pfoa1c = PFOA,
         pfos1c = PFOS,
         etpfosa1c = Et_PFOSA_AcOH,
         pfosa1c = PFOSA,
         mepfosa2c = MeFOSAA_MT,
         pfhxs2c = PFHxS_MT,
         pfna2c = PFNA_MT,
         pfda2c = PFDA_MT,
         pfunda2c = PFUnDA_MT,
         sm_pfos2c = Sm_PFOS_MT,
         n_pfos2c = n_PFOS_MT,
         n_pfoa2c = n_PFOA_MT,
         sb_pfoa2c = Sb_PFOA_MT,
         pfoa2c = PFOA_MT,
         pfos2c = PFOS_MT,
         meansbp = sys_mean_mbp17,                          #meanSBP at midlife
         meandbp = dias_mean_mbp17,                         #meanDBP at midlife
         hdl = HDLC_mg_dL_mom_mt,                           #hdl at midlife
         ldl = Cholesterol_mg_dL_mom_mt-HDLC_mg_dL_mom_mt,  #ldl at midlife
         chol = Cholesterol_mg_dL_mom_mt,                   #chol at midlife
         trig = Triglyceride_mg_dL_mom_mt,                  #trig at midlife
         glucose = glucose_mgl17,                           #fasting glucose at midlife 
         glucosetol = glucose_tolerance_mgl17,              #glucose tolerance at midlife
         hba = pct_HbA1c_mom_mt,                            #hba1c at midlife
         insulin = Insulin_uU_mL_mom_mt,                    #insulin at midlife
         homair = HOMA_IR_mom_MT,                           #homair at midlife
         fasting = blood_fasting_mbl17,                     #fasting status
         diab_mell_med_in17=ifelse(diab_mell_med_in17=='NA'|is.na(diab_mell_med_in17),2,diab_mell_med_in17),
         dm_med = ifelse(diab_mell_med_in17==1,1,0)) %>%    #diabetes meds at midlife
  select(c(aid,age1c,age2c,edu,race,hhinc1c,hhinc2c,nullip, #maritalstat1c, maritalstat2c,
           highbp,t1d,t2d,htnpreg,egfr,bmiprepreg,cig1c,cig2c,bmi_midlf,menopause,
           dm_med, parity2c,parity1c, 
           #bfdur2c, #bfdur1c, 
           pfos1c, pfoa1c, pfhxs1c, pfna1c, pfda1c, mepfosa1c,
           pfos2c, pfoa2c, pfhxs2c, pfna2c, pfda2c, mepfosa2c,
           #sercreat,albumin,gestage,
           #dash,#bfdur,#alcohol1c,alcohol2c,
           #meansbp,meandbp,hdl,ldl,chol,trig,
           #glucose,homair,
           #glucosetol,
           hba,insulin)) %>% 
  mutate_at(vars(matches("pf")), .funs = list(log = ~log(.))) %>% 
  #mutate_at(c('t1d', 't2d', 'dm_med', 'menopause'), as.factor) %>% 
  drop_na() %>% 
  filter(t1d!=1) %>% #remove individuals with diabetes
  filter(t2d!=1) #remove individuals with diabetes


summary(pfas)
names(pfas) <- tolower(names(pfas))

##Labels------------------------------------------------------------------------
pfas$edu = factor(pfas$edu, levels = 1:5, 
                  labels = c("Less than 12th Grade", "High School Diploma/GED",
                             "Some College", "Bachelor Degree", "Graduate Degree"))
pfas$maritalstat= factor(pfas$maritalstat, levels = 1:6, 
                         labels = c("Married", "Living with partner",
                                    "Never Married", "Divorced", "Separated",
                                    "Widowed")) 
pfas$hhinc1c = factor(pfas$hhinc1c, levels = 0:1, 
                    labels = c("\u226470k", "\u226570k")) 
pfas$hhinc2c = factor(pfas$hhinc2c, levels = 0:1, 
                      labels = c("\u226470k", "\u226570k")) 
pfas$parity1c = factor(pfas$parity1c, levels = 0:5, 
                     labels = c('Nulliparous', '1 Birth', '2 Births', '3 Births', '4 Births', '\u22655 Births')) 
pfas$parity2c = factor(pfas$parity2c, levels = 0:5, 
                       labels = c('Nulliparous', '1 Birth', '2 Births', '3 Births', '4 Births', '\u22655 Births')) 
pfas$race= factor(pfas$race, levels = 0:5, 
                  labels = c("White", "Black", "Hispanic", "Asian", "American Indian","Other")) 
pfas$cig1c= factor(pfas$cig1c, levels = 0:2, 
                 labels = c("Never", "Former", "Current")) 
pfas$cig2c= factor(pfas$cig2c, levels = 0:2, 
                   labels = c("Never", "Former", "Current")) 


##Save RDA File-----------------------------------------------------------------
save(pfas, file="vivapfas.rda")
write.csv(pfas, file = "vivapfas.csv")


##Participant Flowchart---------------------------------------------------------
a <- pfas
b <- a %>% 
  select(c(aid,pfos1c, pfoa1c, pfhxs1c, pfna1c, pfda1c, mepfosa1c,
           pfos2c, pfoa2c, pfhxs2c, pfna2c, pfda2c, mepfosa2c,
           age1c,age2c,race,edu,hhinc1c,hhinc2c,bmiprepreg,bmi_midlf,
           cig1c,cig2c, parity1c, parity2c, dm_med, menopause,egfr,
           highbp,t1d,t2d,htnpreg,hba,insulin))%>% drop_na() %>% 
  filter(t1d!=1) %>% #remove individuals with diabetes
  filter(t2d!=1) #remove individuals with diabetes
           #glucose,homair,#glucosetol,
           #hba,insulin)) %>% 
  mutate_at(vars(matches("pf")), .funs = list(log = ~log(.))) %>% 
  drop_na() %>% 
  filter(t1d!=1) %>% #remove individuals with diabetes
  filter(t2d!=1) #remove individuals with diabetes


##Table 1. Participant Characteristics------------------------------------------
var_label(pfas) = list(age1c="Age", age2c="Age",
                            race="Race/Ethnicity", edu="Education", hhinc1c="Income",hhinc2c="Income",
                            bmiprepreg="Pre-Pregnancy BMI", bmi_midlf="BMI", #dash="DASH Diet Score",
                            cig1c="Smoking Status", cig2c="Smoking Status", egfr='eGFR',
                            #bfdur="Breastfeeding Duration", 
                            parity1c="Parity", parity2c="Parity",#maritalstat="Marital Status", 
                            menopause="Menopause",
                            #sercreat="Serum Creatinine", albumin="Albumin",
                            #meansbp="Systolic BP (mmHg)", meandbp="Diastolic BP (mmHg)", 
                            #hdl="HDL (mg/dL)", ldl="LDL (mg/dL)", trig="Triglycerides (mg/dL)", chol="Total Cholesterol (mg/dL)",
                            hba="HbA1c", insulin="Insulin",
                            #glucose="Fasting Glucose", homair="HOMA-IR", 
                            #glucosetol="Glucose Tolerance", 
                            dm_med='DM Medications',
                            #bp_med='BP Medications', lipid_med="Lipid Medications", 
                            highbp="History of High BP",
                            t1d="Type 1 Diabetes", t2d="Type 2 Diabetes", htnpreg="Hypertensive Disorders in Pregnancy",
                            pfos1c="PFOS", pfoa1c="PFOA", pfhxs1c="PFHxS", pfna1c="PFNA", pfda1c="PFDA", mepfosa1c="MEPFOSA",
                            pfos2c="PFOS", pfoa2c="PFOA", pfhxs2c="PFHxS", pfna2c="PFNA", pfda2c="PFDA", mepfosa2c="MEPFOSA")

vars <- c("age1c", "age2c", "race", "edu", "hhinc1c", "hhinc2c","bmiprepreg", "bmi_midlf", "cig1c","cig2c", #"dash",
          #"bfdur",
          "parity1c", "parity2c", "menopause",
          "htnpreg","highbp","t1d", "t2d",#"sercreat", "albumin",
          #"meansbp", "meandbp", "hdl", "ldl", "trig", "chol",
          "hba", "insulin",
          "glucose", "homair",
          "glucosetol", "dm_med",#"bp_med", "lipid_med", 
          "pfos1c", "pfoa1c", "pfhxs1c", "pfna1c", "pfda1c", "mepfosa1c",
          "pfos2c", "pfoa2c", "pfhxs2c", "pfna2c", "pfda2c", "mepfosa2c")
tab1 <- CreateTableOne(vars = vars,data = pfas)#, addOverall=T)


print(tab1, nonnormal = vars, varLabels=F, contDigits=2, formatOptions = list(big.mark = ","), showAllLevels=T)
write.csv(print(tab1, nonnormal = vars, varLabels=T, contDigits=2, formatOptions = list(big.mark = ","), showAllLevels=T), 'tab1.csv')


table(pfas$race)
pfas %>% group_by(exam, race) %>% 
  summarise(n = n())

##ICC---------------------------------------------------------------------------
icc <- pfas_long %>%
  select(c("aid", "exam","mepfosa", "pfda","pfhxs", "pfna","pfoa", "pfos"))
pfa <-c("mepfosa", "pfda","pfhxs", "pfna","pfoa", "pfos"); length(pfa) #11
pfas.icc <- subset(icc, complete.cases(icc[,pfa])); dim(pfas.icc) #753
f <- subset(pfas.icc, aid %in% names(which(table(aid) > 1))) %>% ungroup
f %>% group_by(aid)
str(f)

a <- f %>%
  group_by(aid) %>% 
  summarize_at(c("mepfosa", "pfda","pfhxs", "pfna","pfoa", "pfos"),
               .funs = list(avg = ~mean(., na.rm=T))) %>% 
  summarize_at(c("mepfosa_avg", "pfda_avg","pfhxs_avg", "pfna_avg","pfoa_avg", "pfos_avg"),
               .funs = list(btwvar = ~var(., na.rm=T))) %>% 
  gather(key=pfasbtwvar, value=btwvar) # %>% ungroup()

b <- f %>% 
  summarize_at(c("mepfosa", "pfda","pfhxs", "pfna","pfoa", "pfos"),
               .funs = list(totalvar = ~var(., na.rm=T)))%>% 
  gather(key=pfastotvar, value=totalvar) %>% 
  cbind(a) %>% 
  mutate(icc=btwvar/totalvar)



##Correlation Table-------------------------------------------------------------
plot <- pfas_long %>% select(c(aid, exam, pfos, pfoa, pfhxs, pfna, pfda, mepfosa))%>% 
  drop_na() %>% 
  subset(aid %in% names(which(table(aid) > 1)))  %>% #group_by(aid)
  pivot_wider(
    names_from = exam,
    names_sep = ".",
    values_from = c(pfos, pfoa, pfhxs, pfna, pfda, mepfosa)) %>% 
  mutate_at(vars(matches("pf")), .funs = list(log = ~log(.))) %>% 
  select(contains('log'))


summary(plot)

plot <- na.omit(plot)
corr <- round(cor(plot, method='spearman'), 2)
head(corr[, 1:12])
p.mat <- cor_pmat(plot)
head(p.mat[, 1:12])

tiff(file=paste0("D://KM//Columbia//Project VIVA//data", "/corr_tvpfaslog.tiff"), res=600, 
     width = 7, height = 7, units = 'in', compression = c("lzw"))
ggcorrplot(corr, outline.color='white', lab=T)+
  ggplot2::theme(
    axis.text.x = element_text(angle=0, color='black'),
    axis.text.y = element_text(angle=0, color='black'))
dev.off()

library(corrplot)

par(mar=c(4,5,1,1))
tiff(file="D://KM//Columbia//Project VIVA//data//corr_tvpfaslog.tiff", units="in", width=10, height=10, res=600) 
corr %>%  
  corrplot(method = "circle", type = "upper", diag=FALSE, is.corr = FALSE)
dev.off()

library(corrplot)
tiff(file="D://KM//Columbia//Project VIVA//data//corr_tvpfaslog_circlesvalues.tiff", units="in", width=10, height=10, res=600) 
par(mar=c(4,5,1,1))
cor(plot, method='spearman')%>%  
  corrplot(method = "circle", type = "upper", diag=FALSE,addCoef.col = 'black')
dev.off()




plot <- pfas_long %>%
  filter(exam==1) %>% 
  dplyr::select(c("pfos", "pfoa", "pfhxs", "pfna","pfda","mepfosa","et_pfosa_acoh"))#, 
#"pfunda", "sb_pfoa","sm_pfos", "n_pfos","n_pfoa"))
colnames(plot) <- c( "PFOS", "PFOA", "PFHXS", "PFNA", "PFDA", "MEPFOSAA", "ETPFOSA")

plot <- pfas_long %>%
  filter(exam==1) %>% 
  mutate_at(c("pfos", "pfoa", "pfhxs", "pfna","pfda","me[fosa","et_pfosa_acoh"),
            .funs = list(log = ~log(.))) %>% 
  select(c("pfos_log", "pfoa_log", "pfhxs_log", "pfna_log","pfda_log","mefosaa_log","et_pfosa_acoh_log"))
colnames(plot) <- c( "PFOS", "PFOA", "PFHXS", "PFNA", "PFDA", "MEPFOSAA", "ETPFOSA")
summary(plot)

plot <- na.omit(plot)
corr <- round(cor(plot, method='spearman'), 2)
head(corr[, 1:7])
p.mat <- cor_pmat(plot)
head(p.mat[, 1:7])


library(corrplot)

c <- cor.mtest(plot, conf.level = 0.95)

corrplot(corr, p.mat = p.mat, insig = "p-value", method = "circle", type = "upper", diag=FALSE, is.corr = FALSE)


par(mar=c(4,5,1,1))
tiff(file="D://KM//Columbia//Project VIVA//coorplot_logpfas.tiff", units="in", width=10, height=10, res=600) 
corrplot(corr, p.mat = p.mat, method = 'circle', type = 'lower', insig ='blank',
         addCoef.col ='black', number.cex = 0.8, diag = FALSE)
dev.off()


par(mar=c(4,5,1,1))
a <- cor(corr)
colnames(a) <- c("PFOS", "PFOA", "PFHXS", "PFNA", "PFDA", "MEFOSAA", "ETPFOSA")
rownames(a) <- c("PFOS", "PFOA", "PFHXS", "PFNA", "PFDA", "MEFOSAA", "ETPFOSA")
b <- a %>%  
  corrplot(a$r, p.mat = a$p, insig = "p-value", method = "circle", type = "upper", diag=FALSE, is.corr = FALSE)

tiff(file="D:/KM/Columbia/MESA/Data/coorplot_new_2.tiff", units="in", width=10, height=10, res=600) 
a %>%  
  corrplot(method = "circle", type = "upper", diag=FALSE, is.corr = FALSE)
dev.off()

getwd()
##Are there associations between time-varying covariates and PFAS0?-------------
rg <- map_dfr(pfas[,c("pfos1c","pfoa1c","pfhxs1c","pfna1c","pfda1c", "mepfosa1c")],
              function(x) quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=T)) %>% 
  mutate(p10=log(`10%`), p25=log(`25%`), p50=log(`50%`), p75=log(`75%`), p90=log(`90%`), iqr=p75-p25,
         term=c("log(pfos1c)","log(pfoa1c)","log(pfhxs1c)","log(pfna1c)","log(pfda1c)","log(mepfosa1c)")) %>% 
  select(c(term, iqr))

##simple linear regression models at baseline
dependent_variables <- c("pfos1c","pfoa1c","pfhxs1c","pfna1c","pfda1c", "mepfosa1c")
a <- map_df(dependent_variables, function(y) {
  map_df(names(pfas)[17], function(x) {
    formula <- as.formula(paste0(y, " ~ ", x, ""))
    model <- glm(formula, data = pfas)
    tidy(model, conf.int=T) %>% 
      mutate(outcome = y)})})%>% 
  filter(grepl('bmi', term))

dependent_variables <- c("bmi_midlf")
b <- map_df(dependent_variables, function(y) {
  map_df(names(pfas)[24:29], function(x) {
    formula <- as.formula(paste0(y, " ~ ", x, ""))
    model <- glm(formula, data = pfas)
    tidy(model, conf.int=T) %>% 
      mutate(outcome = y)})})%>% 
  filter(grepl('pf', term)) 

names(pfas)
dependent_variables <- c("pfos1c","pfoa1c","pfhxs1c","pfna1c","pfda1c", "mepfosa1c",
                         "pfos2c","pfoa2c","pfhxs2c","pfna2c","pfda2c", "mepfosa2c") 
c <- map_df(dependent_variables, function(y) {##maybe don't need to include dmmed as time varying covariate
  map_df(names(pfas)[18:22], function(x) {##parity and dm_meds not significant
    formula <- as.formula(paste0(y, " ~ ", x, ""))
    model <- glm(formula, data = pfas)
    tidy(model, conf.int=T) %>% 
      mutate(outcome = y)})})%>% 
  filter(grepl('menopause|lipid_med|bp_med|dm_med|parity', term)) 

dependent_variables <- c("menopause","bp_med","lipid_med","dm_med", "parity2c") ##can't code parity as 0/1
d <- map_df(dependent_variables, function(y) {
  map_df(names(pfas)[24:29], function(x) {
    formula <- as.formula(paste0(y, " ~ ", x, ""))
    model <- glm(formula, data = pfas, family=binomial(link=logit))
    tidy(model, conf.int=T) %>% 
      mutate(outcome = y)})})%>% 
  filter(grepl('pf', term)) 

##Model 1: Model similar to Causal model----------------------------------------
##Pregnancy PFAS and Midlife Glycemia
rg <- map_dfr(pfas[,c("pfos1c","pfoa1c","pfhxs1c","pfna1c","pfda1c", "mepfosa1c")],
              function(x) quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=T)) %>% 
  mutate(p10=log(`10%`), p25=log(`25%`), p50=log(`50%`), p75=log(`75%`), p90=log(`90%`), iqr=p75-p25,
         term=c("log(pfos1c)","log(pfoa1c)","log(pfhxs1c)","log(pfna1c)","log(pfda1c)","log(mepfosa1c)")) %>% 
  select(c(term, iqr))

names(pfas)
summary(pfas)

dependent_variables <- c("hba", "insulin") #n=399
dependent_variables <- c("homair", "glucose") #n=374
dependent_variables <- c("glucosetol") #n=176
a <- map_df(dependent_variables, function(y) {
  map_df(names(pfas)[24:29], function(x) {
    formula <- as.formula(paste0(y, " ~ log(", x, ") + age1c + race + edu + bmiprepreg + hhinc1c + cig1c + parity1c"))
    model <- glm(formula, data = pfas) #dash, bfdur, egfr not significant in model
    tidy(model, conf.int=T) %>% 
      mutate(outcome = y)})})%>% 
  filter(grepl('log', term)) %>% 
  left_join(rg) %>% 
  mutate(est=exp(iqr)*estimate,
         lci=exp(iqr)*conf.low,
         hci=exp(iqr)*conf.high) %>% 
  mutate(across(c('est', 'lci', 'hci'), round, 2)) %>% 
  unite("est.gmr", lci, hci, sep = ", ", remove = FALSE) %>%
  mutate("est.gmr" = paste0(est, " (", est.gmr, ")")) %>% 
  select(c(term, est.gmr, outcome))

table(pfas$bfdur)

##Midlife PFAS and Midlife Glycemia
rg <- map_dfr(pfas[,c("pfos2c","pfoa2c","pfhxs2c","pfna2c","pfda2c", "mepfosa2c")],
              function(x) quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=T)) %>% 
  mutate(p10=log(`10%`), p25=log(`25%`), p50=log(`50%`), p75=log(`75%`), p90=log(`90%`), iqr=p75-p25,
         term=c("log(pfos2c)","log(pfoa2c)","log(pfhxs2c)","log(pfna2c)","log(pfda2c)","log(mepfosa2c)")) %>% 
  select(c(term, iqr))

names(pfas)
summary(pfas)

dependent_variables <- c("hba", "insulin") #n=399
dependent_variables <- c("homair", "glucose") #n=374
dependent_variables <- c("glucosetol") #n=176
b <- map_df(dependent_variables, function(y) {
  map_df(names(pfas)[28:33], function(x) {
    formula <- as.formula(paste0(y, " ~ log(", x, ") + age1c + race + edu + bmiprepreg + bmi_midlf + hhinc1c + cig1c + nullip + dm_med + menopause"))
    model <- glm(formula, data = pfas)
    tidy(model, conf.int=T) %>% 
      mutate(outcome = y)})})%>% 
  filter(grepl('log', term)) %>% 
  left_join(rg) %>% 
  mutate(est=exp(iqr)*estimate,
         lci=exp(iqr)*conf.low,
         hci=exp(iqr)*conf.high) %>% 
  mutate(across(c('est', 'lci', 'hci'), round, 2)) %>% 
  unite("est.gmr", lci, hci, sep = ", ", remove = FALSE) %>%
  mutate("est.gmr" = paste0(est, " (", est.gmr, ")")) %>% 
  select(c(term, est.gmr, outcome))


##Model 2: Using covariates relevant to exposure and outcome timepoint----------
##Pregnancy PFAS and Midlife Glycemia
rg <- map_dfr(pfas[,c("pfos1c","pfoa1c","pfhxs1c","pfna1c","pfda1c", "mepfosa1c")],
              function(x) quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=T)) %>% 
  mutate(p10=log(`10%`), p25=log(`25%`), p50=log(`50%`), p75=log(`75%`), p90=log(`90%`), iqr=p75-p25,
         term=c("log(pfos1c)","log(pfoa1c)","log(pfhxs1c)","log(pfna1c)","log(pfda1c)","log(mepfosa1c)")) %>% 
  select(c(term, iqr))

names(pfas)
summary(pfas)

dependent_variables <- c("hba", "insulin") #n=399
dependent_variables <- c("homair", "glucose") #n=374
dependent_variables <- c("glucosetol") #n=176
a <- map_df(dependent_variables, function(y) {
  map_df(names(pfas)[22:27], function(x) {
    formula <- as.formula(paste0(y, " ~ log(", x, ") + age1c + race + edu + bmiprepreg + hhinc1c + cig1c + nullip"))
    model <- glm(formula, data = pfas)
    tidy(model, conf.int=T) %>% 
      mutate(outcome = y)})})%>% 
  filter(grepl('log', term)) %>% 
  left_join(rg) %>% 
  mutate(est=exp(iqr)*estimate,
         lci=exp(iqr)*conf.low,
         hci=exp(iqr)*conf.high) %>% 
  mutate(across(c('est', 'lci', 'hci'), round, 2)) %>% 
  unite("est.gmr", lci, hci, sep = ", ", remove = FALSE) %>%
  mutate("est.gmr" = paste0(est, " (", est.gmr, ")")) %>% 
  select(c(term, est.gmr, outcome))



##Midlife PFAS and Midlife Glycemia
rg <- map_dfr(pfas[,c("pfos2c","pfoa2c","pfhxs2c","pfna2c","pfda2c", "mepfosa2c")],
              function(x) quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=T)) %>% 
  mutate(p10=log(`10%`), p25=log(`25%`), p50=log(`50%`), p75=log(`75%`), p90=log(`90%`), iqr=p75-p25,
         term=c("log(pfos2c)","log(pfoa2c)","log(pfhxs2c)","log(pfna2c)","log(pfda2c)","log(mepfosa2c)")) %>% 
  select(c(term, iqr))

names(pfas)

summary(pfas)

dependent_variables <- c("hba", "insulin") #n=399
dependent_variables <- c("homair", "glucose") #n=374
dependent_variables <- c("glucosetol") #n=176
b <- map_df(dependent_variables, function(y) {
  map_df(names(pfas)[28:33], function(x) {
    formula <- as.formula(paste0(y, " ~ log(", x, ") + age2c + race + edu + bmi_midlf + hhinc2c + cig2c + nullip + dm_med + menopause"))
    model <- glm(formula, data = pfas)
    tidy(model, conf.int=T) %>% 
      mutate(outcome = y)})})%>% 
  filter(grepl('log', term)) %>% 
  left_join(rg) %>% 
  mutate(est=exp(iqr)*estimate,
         lci=exp(iqr)*conf.low,
         hci=exp(iqr)*conf.high) %>% 
  mutate(across(c('est', 'lci', 'hci'), round, 2)) %>% 
  unite("est.gmr", lci, hci, sep = ", ", remove = FALSE) %>%
  mutate("est.gmr" = paste0(est, " (", est.gmr, ")")) %>% 
  select(c(term, est.gmr, outcome))

##G-Function Bootstrap-Marginal Effects Model-----------------------------------
##baseline covariates: age, race, edu, hhinc, cig, nullip
##time-varying covariates: bmi_midlf, menopause, dm_meds
summary(lm_l4)
str(pfas)
gfun <- function(pfas){ 
  lm_l1 <-  lm(bmi_midlf ~  age1c + race + edu + hhinc1c + cig1c + bmiprepreg + #nullip + dm_med + menopause +
                 pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log, 
               dat = pfas)
  
  lm_l2 <-  glm(dm_med ~  age1c + race + edu + hhinc1c + cig1c + bmiprepreg + #nullip + bmi_midlf + menopause +
                  pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log, 
                dat = pfas, family  = binomial(link = "probit"))
  
  lm_l3 <-  glm(menopause ~  age1c + race + edu + hhinc1c + cig1c + bmiprepreg + #nullip + bmi_midlf + dm_med +
                  pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log,
                dat = pfas, family  = binomial(link = "probit"))
  
  lm_l4 <-  glm(parity2c ~  age1c + race + edu + hhinc1c + cig1c + bmiprepreg + #nullip + bmi_midlf + dm_med +
                  pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log,
                dat = pfas, family  = poisson(link = "log"))
  
  lm_y <-  lm(hba ~  age1c + race + edu + hhinc1c + cig1c + bmiprepreg +
                bmi_midlf + dm_med + menopause + parity2c + 
                pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log +
                pfos2c_log + pfoa2c_log + pfhxs2c_log + pfna2c_log + pfda2c_log + mepfosa2c_log, 
              dat = pfas)
  
  A <- dplyr::select(pfas, pfos1c_log, pfoa1c_log, pfhxs1c_log, pfna1c_log, pfda1c_log, mepfosa1c_log)##,
  #pfos2c_log, pfoa2c_log, pfhxs2c_log, pfna2c_log, pfda2c_log, mepfosa2c_log)
  
  B <- dplyr::select(pfas,pfos2c_log, pfoa2c_log, pfhxs2c_log, pfna2c_log, pfda2c_log, mepfosa2c_log)
  
  a <-   apply(A, 2, quantile, probs=0.50)        ##reference
  astar <-   apply(B, 2, quantile, probs=0.75)    ##intervention ##75th percentile at second timepoint
  
  
  boot <- foreach(i=1:1000, .combine='rbind') %dopar% {
    ind <- sample(1:391, replace=TRUE)
    
    # Calculate L_a and L_a*
    tmp_a <- tmp_astar <- model.frame(lm_l1)[ind,]
    tmp_a[,colnames(A)] <- a
    tmp_astar[,colnames(A)] <- astar
    L1a.true <- predict(lm_l1, newdata = tmp_a)  #la and la star for each L
    L1astar.true <- predict(lm_l1, newdata = tmp_astar)
    
    
    tmp_a <- tmp_astar <- model.frame(lm_l2)[ind,]
    tmp_a[,colnames(A)] <- a
    tmp_astar[,colnames(A)] <- astar
    L2a.true <- predict(lm_l2, newdata = tmp_a, type='response') #la and la star for each L
    L2astar.true <- predict(lm_l2, newdata = tmp_astar, type='response')
    
    tmp_a <- tmp_astar <- model.frame(lm_l3)[ind,]
    tmp_a[,colnames(A)] <- a
    tmp_astar[,colnames(A)] <- astar
    L3a.true <- predict(lm_l3, newdata = tmp_a, type='response')  #la and la star for each L
    L3astar.true <- predict(lm_l3, newdata = tmp_astar, type='response')
    
    tmp_a <- tmp_astar <- model.frame(lm_l4)[ind,]
    tmp_a[,colnames(A)] <- a
    tmp_astar[,colnames(A)] <- astar
    L4a.true <- predict(lm_l4, newdata = tmp_a, type='response')  #la and la star for each L
    L4astar.true <- predict(lm_l4, newdata = tmp_astar, type='response')
    
    # Predict the value of L2 and L3 using the predicted probabilities
    L2a.samp      <- rbinom(391, size=1, prob=L2a.true)
    L2astar.samp  <- rbinom(391, size=1, prob=L2astar.true)
    
    L3a.samp      <- rbinom(391, size=1, prob=L3a.true)
    L3astar.samp  <- rbinom(391, size=1, prob=L3astar.true)
    
    L4a.samp      <- rpois(391, lambda=L4a.true)
    L4astar.samp  <- rpois(391, lambda=L4astar.true)
    
    # Predict Y
    tmp_y_a <- tmp_y_astar <- model.frame(lm_y)[ind,]
    tmp_y_a[,colnames(A)] <- a
    tmp_y_astar[,colnames(A)] <- astar
    tmp_y_a[,'bmi_midlf'] <- L1a.true
    tmp_y_astar[,'bmi_midlf'] <- L1astar.true
    tmp_y_a[,'dm_med'] <- L2a.samp
    tmp_y_astar[,'dm_med'] <- L2astar.samp
    tmp_y_a[,'menopause'] <- L3a.samp
    tmp_y_astar[,'menopause'] <- L3astar.samp
    tmp_y_a[,'parity2c'] <- L4a.samp
    tmp_y_astar[,'parity2c'] <- L4astar.samp
    
    YaLa.pred   <- predict(lm_y, newdata = tmp_y_a, type = "response")
    YastarLastar.pred <- predict(lm_y, newdata = tmp_y_astar, type = "response")
    effect <- mean(YastarLastar.pred) - mean(YaLa.pred)
    effect
  }
  
  effect <- round(quantile(boot[,1], probs=0.5),3)
  CI_up <- round(quantile(boot[,1], probs=0.975),3)
  CI_low <- round(quantile(boot[,1], probs=0.025),3)
  #return(list(effect=effect, CI_up=CI_up, CI_low=CI_low))
  return(sprintf("%.2f (%.2f, %.2f)", effect, CI_low, CI_up))
}

gfun(pfas)

##G-Function Fixing Covariates-Conditional Model--------------------------------
##baseline covariates: age, race, edu, hhinc, cig, nullip
##time-varying covariates: bmi_midlf, menopause, dm_meds

gfun <- function(pfas){ 
  lm_l1 <-  lm(bmi_midlf ~  age1c + race + edu + hhinc1c + cig1c + nullip + bmiprepreg + #dm_med + menopause +
                 pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log, 
                 dat = pfas)
  
  lm_l2 <-  glm(dm_med ~  age1c + race + edu + hhinc1c + cig1c + nullip + bmiprepreg + #bmi_midlf + menopause +
                  pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log, 
               dat = pfas, family  = binomial(link = "probit"))
  
  lm_l3 <-  glm(menopause ~  age1c + race + edu + hhinc1c + cig1c + nullip + bmiprepreg + #bmi_midlf + dm_med +
               pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log,
               dat = pfas, family  = binomial(link = "probit"))
  
  #lm_l4 <-  glm(parity2c ~  age1c + race + edu + hhinc1c + cig1c + nullip + bmiprepreg + #bmi_midlf + dm_med +
                  #pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log,
                #dat = pfas, family  = binomial(link = "probit")) ##poisson
  
  lm_y <-  lm(hba ~  age1c + race + edu + hhinc1c + cig1c + nullip + bmiprepreg +
                 bmi_midlf + dm_med + menopause + 
                pfos1c_log + pfoa1c_log + pfhxs1c_log + pfna1c_log + pfda1c_log + mepfosa1c_log +
                pfos2c_log + pfoa2c_log + pfhxs2c_log + pfna2c_log + pfda2c_log + mepfosa2c_log, 
               dat = pfas)
  
  K = 1000
  
  A <- dplyr::select(pfas, pfos1c_log, pfoa1c_log, pfhxs1c_log, pfna1c_log, pfda1c_log, mepfosa1c_log)##,
                       #pfos2c_log, pfoa2c_log, pfhxs2c_log, pfna2c_log, pfda2c_log, mepfosa2c_log)
  
  B <- dplyr::select(pfas,pfos2c_log, pfoa2c_log, pfhxs2c_log, pfna2c_log, pfda2c_log, mepfosa2c_log)
  
  a <-   apply(A, 2, quantile, probs=0.50)        ##reference
  astar <-   apply(B, 2, quantile, probs=0.75)    ##intervention ##75th percentile at second timepoint
  
  # Calculate L_a*
  tmp = rbind(c(mean(pfas$age1c), 0, 1, 1, 0, 0, mean(pfas$bmiprepreg), a[1:6]), ##which group is reference
                                    c(mean(pfas$age1c), 0, 1, 1, 0, 0, mean(pfas$bmiprepreg), astar[1:6])) 
  
  
  newdat_l1 = as.data.frame(tmp) ##bootstrap from here select indices
  colnames(newdat_l1) <- c("age1c", "race", "edu", "hhinc1c", "cig1c", "nullip", "bmiprepreg",
                           "pfos1c_log", "pfoa1c_log", "pfhxs1c_log", "pfna1c_log", "pfda1c_log", "mepfosa1c_log")
  
  ##Predictions
  L1a.true <- predict(lm_l1, newdata = newdat_l1)[1]  #la and la star for each L
  L1astar.true <- predict(lm_l1, newdata = newdat_l1)[2]
  
  L2a.true <- predict(lm_l2, newdata = newdat_l1, type='response')[1] #la and la star for each L
  L2astar.true <- predict(lm_l2, newdata = newdat_l1, type='response')[2]
  
  L3a.true <- predict(lm_l3, newdata = newdat_l1, type='response')[1]  #la and la star for each L
  L3astar.true <- predict(lm_l3, newdata = newdat_l1, type='response')[2]
  
  #L4a.true <- predict(lm_l4, newdata = newdat_l1, type='response')[1]  #la and la star for each L
  #L4astar.true <- predict(lm_l4, newdata = newdat_l1, type='response')[2]
  
  sigma.samp1   <- summary(lm_l1)$sigma ##variance, binomial does not have variance
  
  L1a.samp      <- L1a.true     + rnorm(K, sd=sigma.samp1) 
  L1astar.samp  <- L1astar.true + rnorm(K, sd=sigma.samp1) 

  L2a.samp      <- rbinom(K, size=1, prob=L2a.true)
  L2astar.samp  <- rbinom(K, size=1, prob=L2astar.true)

  L3a.samp      <- rbinom(K, size=1, prob=L3a.true)
  L3astar.samp  <- rbinom(K, size=1, prob=L3astar.true)
  
  #L4a.samp      <- L4a.true     + rnorm(K, sd=sigma.samp4) 
  #L4astar.samp  <- L4astar.true + rnorm(K, sd=sigma.samp4) 
  
  ##L1asamp, L2asamp
  YaLa.samp         <- cbind(mean(pfas$age1c), 0, 1, 1, 0, 0, mean(pfas$bmiprepreg), 
                             L1a.samp, L2a.samp, L3a.samp,matrix(a, nrow=K, ncol=length(a), byrow=TRUE),
                             median(pfas$pfos2c_log),median(pfas$pfoa2c_log),median(pfas$pfhxs2c_log),
                             median(pfas$pfna2c_log),median(pfas$pfda2c_log),median(pfas$mepfosa2c_log))
  YastarLastar.samp <- cbind(mean(pfas$age1c), 0, 1, 1, 0, 0, mean(pfas$bmiprepreg), 
                             L1astar.samp, L2astar.samp, L3astar.samp, matrix(astar, nrow=K, ncol=length(a), byrow=TRUE),
                             median(pfas$pfos2c_log),median(pfas$pfoa2c_log),median(pfas$pfhxs2c_log),
                             median(pfas$pfna2c_log),median(pfas$pfda2c_log),median(pfas$mepfosa2c_log))
  
  newdat_YaLa <- as.data.frame(YaLa.samp)
  colnames(newdat_YaLa) <- c("age1c", "race", "edu", "hhinc1c", "cig1c", "nullip", "bmiprepreg", "bmi_midlf", "dm_med", "menopause",
                             "pfos1c_log", "pfoa1c_log", "pfhxs1c_log", "pfna1c_log", "pfda1c_log", "mepfosa1c_log",
                             "pfos2c_log", "pfoa2c_log", "pfhxs2c_log", "pfna2c_log", "pfda2c_log", "mepfosa2c_log")
  
  newdat_YastarLastar <- as.data.frame(YastarLastar.samp)
  colnames(newdat_YastarLastar) <- c("age1c", "race", "edu", "hhinc1c", "cig1c", "nullip", "bmiprepreg","bmi_midlf", "dm_med", "menopause",
                                     "pfos1c_log", "pfoa1c_log", "pfhxs1c_log", "pfna1c_log", "pfda1c_log", "mepfosa1c_log",
                                     "pfos2c_log", "pfoa2c_log", "pfhxs2c_log", "pfna2c_log", "pfda2c_log", "mepfosa2c_log")
  
  YaLa.pred   <- predict(lm_y, newdata = newdat_YaLa, type = "response")
  YastarLastar.pred <- predict(lm_y, newdata = newdat_YastarLastar, type = "response")
  
  # For parallel computing (useful for bootstrap)
  # Check how many cores your computer has
  # cl <- makeCluster(detectCores(logical = TRUE))
  registerDoParallel(10)
  
  # Bootstrap and save the results pasting them row by row (with rbind)
  boot <- foreach(i=1:1000, .combine='rbind') %dopar% {
    ind <- sample(1:389, replace=TRUE)  ##randomly sample N indexes
    #pred.data.astar.m1 <- pred.data.a.m1 <- model.frame(M[[1]])[ind,]  
    pred.data.astar.l1 <- pred.data.a.l1 <- newdat_l1[ind,]
    pred.data.astar.l2 <- pred.data.a.l2 <- model.frame(lm_l2[[1]])[ind,]
    pred.data.astar.l3 <- pred.data.a.l3 <- model.frame(lm_l3[[1]])[ind,]
    
    model.frame(lm_l1)[[1]][ind,]
    model.frame(lm_l1)
    
    pred.data.astar.m1[, treat] <- treat.value
    pred.data.a.m1[, treat] <- control.value
    
    m1mat.astar <- model.matrix(terms(lm_l1[[1]]), data = pred.data.astar.l1)
    m1mat.a <- model.matrix(terms(M[[1]]), data = pred.data.a.m1)
    
    PredictM_astar <- tcrossprod(MModel[[1]], m1mat.astar)
    PredictM_a <- tcrossprod(MModel[[1]], m1mat.a)
    
    effect <- PredictM_astar - PredictM_a
    effect
  }
  
  # Calculate the effects (quantile bootstrap)
  main_effect <- quantile(boot[,1], 0.5)
  lower_CI <- quantile(boot[,1], 0.025)
  upper_CI <- quantile(boot[,1], 0.975)
  
  
  
  # truth_gformula = quantile(YastarLastar.pred - YaLa.pred, probs=.5) 
  # CI_up = quantile(YastarLastar.pred - YaLa.pred, probs=0.975) 
  # CI_low = quantile(YastarLastar.pred - YaLa.pred, probs=0.025) 
  # 
  # 
  # return(truth_gformula)
  
}
gfun(pfas)

