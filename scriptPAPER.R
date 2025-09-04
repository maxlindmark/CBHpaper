# Herring (Clupea harengus) in subdivisions 25-29 and 32, excluding the Gulf of Riga (central Baltic Sea)

library(ggplot2)
library(plyr)
library(dplyr)
library(gridExtra)
library(lattice)
library(ggridges)
library(stargazer)
library(FSA)
library(readr)
library("nlme")
library(mixedup)
library(MuMIn)
library(PerformanceAnalytics)
library(ggridges)
library(ggrepel)
library(stringr)
library(plotrix)
library(psych)
library(tidyr)
rm(list=ls())  #clean the console

#######################################################################
#  Data exploration
setwd("C:/Users/frma6502/Desktop/Other/CBHpaper")
data_fin2 <- read_csv2("CBH_db_update2025.csv")

ggplot(data_fin2, aes(x = LENGTH)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue") +
  labs(title = "LENGTH Distribution of Fish",
       x = "LENGTH (mm)",
       y = "Count") +
  theme_minimal()+ facet_wrap(~Area)

ggplot(data_fin2, aes(x = agePlus)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue") +
  labs(title = "Age Distribution of Fish",
       x = "Age (years)",
       y = "Count") +
  theme_minimal()+ facet_wrap(~Area)

ggplot(data_fin2, aes(x = weight)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue") +
  labs(title = "Weight Distribution of Fish",
       x = "weight (gr)",
       y = "Count") +
  theme_minimal()+ facet_wrap(~Area)

ggplot(data_fin2, aes(x = Cohort)) +
  geom_histogram(binwidth = 1, color = "black", fill = "skyblue") +
  labs(title = "Distribution of Fish in Cohort",
       x = "Cohort (Year)",
       y = "Count") +
  theme_minimal()+ facet_wrap(~Area)

##############################################################################################
###########################################
###### length #############################
# checks in trend in length at age 
tmp_ageyr <- data_fin2 %>%
  group_by( Year,agePlus) %>%  
  dplyr::summarize(n=validn(LENGTH),
                   mean_len=round(mean(LENGTH),0) ) %>%
  # sd=round(sd(Length.class),1)) %>%
  as.data.frame()
tmp_ageyr 
jpeg(paste("mean len trend.jpeg"),width = 300, height = 200, units = "mm", res = 600)
ggplot(tmp_ageyr, aes(Year, mean_len,colour=as.factor(agePlus ))) + geom_point(size=1.3)  + geom_line(size=0.5)  +theme(legend.position="rigth")+ theme_bw()#+ geom_smooth(se = F)
dev.off()

# Next step: Checks in trend in length-at-age by SD...
# plot the data 
p1b <- ggplot(data_fin2) + geom_boxplot(aes(x = as.factor(agePlus), y = LENGTH, fill = as.factor(Area)))+theme(legend.position="bottom")
# checks in trend in length at age by area
tmp_agearea <- data_fin2 %>%
  group_by(Area,agePlus) %>%  
  dplyr::summarize(n=validn(LENGTH),
                   mean_len=round(mean(LENGTH),0) ) %>%
  # sd=round(sd(Length.class),1)) %>%
  as.data.frame()
tmp_agearea
p2b <- ggplot(tmp_agearea, aes(agePlus, mean_len,colour=as.factor(Area ))) + geom_point(size=1.3) + geom_line(size=0.9)# + geom_smooth(se =F) +theme(legend.position="")

jpeg(paste("mean len age by SD.jpeg"),width = 400, height = 200, units = "mm", res = 600)
gg <-grid.arrange( p1b, p2b,nrow = 1)
print(gg)
dev.off()

# ...and SD in time?
########################################
tmp_area_ageyr <- data_fin2 %>%
  group_by(Year,Area,agePlus) %>%  
  dplyr::summarize(n=validn(LENGTH),
                   mean_len=round(mean(LENGTH),0) ) %>%
  # sd=round(sd(Length.class),1)) %>%
  as.data.frame()
tmp_area_ageyr
jpeg(paste("mean len age by SD trends Yrs.jpeg"),width = 400, height = 200, units = "mm", res = 600)
ggplot(tmp_area_ageyr, aes(Year, mean_len,colour=as.factor(agePlus ))) + geom_point(size=1.3)+ geom_line(size=0.4)+ theme_bw() + facet_wrap(~Area) # + geom_smooth(se =F)
dev.off()



#######################################################################################################
#VBGP LENGHT by SD and Year (+ random effect)

# Apply nonlinear mixed effect model via nlme using SD as a random effect 
## model without random effect
mod1 <-nls(LENGTH ~ Linf*(1-exp(-K*(AgeFraction-t0))) , start=list(Linf=220,K=0.3, t0=-1.95),
          data=(data_fin2  ),lower=c(100,0.1, -4),
          upper=c(350,1,1),algorithm="port")
summary(mod1); AIC(mod1)
p3 <-plot(mod1, Area ~ resid(.), abline = c(0,0));p3
p4 <-plot(mod1, as.factor(Year) ~ resid(.), abline = c(0,0));p4

## by SD
mod2 <- nlme(LENGTH ~ Linf*(1-exp(-k0*(AgeFraction-(t0))))  ,fixed= Linf+k0+t0 ~1 ,random=  list(Area = Linf+k0~1) ,data=data_fin2, start=c(Linf=220,k0=0.3,t0=-1.95))
summary(mod2); AIC(mod2)
p3 <-plot(mod2, Area ~ resid(.), abline = c(0,0));p3
p4 <-plot(mod2, as.factor(Year) ~ resid(.), abline = c(0,0));p4

## by Year 
mod3 <- nlme(LENGTH ~ Linf*(1-exp(-k0*(AgeFraction-(t0))))  ,fixed= Linf+k0+t0 ~1 ,random=  list(Year = Linf+k0~1) ,data=data_fin2, start=c(Linf=220,k0=0.3,t0=-1.95)) 
summary(mod3); AIC(mod3)
p3 <-plot(mod3, Area ~ resid(.), abline = c(0,0));p3
p4 <-plot(mod3, as.factor(Year) ~ resid(.), abline = c(0,0));p4

## by Year in SD
mod4 <- nlme(LENGTH ~ Linf*(1-exp(-k0*(AgeFraction-(t0))))  ,fixed= Linf+k0+t0 ~1 ,random=  list(Area = Linf+k0~1,Year = Linf+k0~1) ,data=data_fin2, start=c(Linf=220,k0=0.3,t0=-1.95)) 
 summary(mod4); AIC(mod4)
 summarise_model(mod4, ci = FALSE)
p3 <-plot(mod4, Area ~ resid(.), abline = c(0,0));p3
p4 <-plot(mod4, as.factor(Year) ~ resid(.), abline = c(0,0));p4

## by Cohort 
mod5 <- nlme(LENGTH ~ Linf*(1-exp(-k0*(AgeFraction-(t0))))  ,fixed= Linf+k0+t0 ~1 ,random=  list(Cohort = Linf+k0~1) ,data=data_fin2, start=c(Linf=220,k0=0.3,t0=-1.95)) 
summary(mod5); AIC(mod5)
summarise_model(mod5, ci = FALSE)
p3 <-plot(mod5, Area ~ resid(.), abline = c(0,0));p3
p4 <-plot(mod5, as.factor(Cohort) ~ resid(.), abline = c(0,0));p4

## by Cohort in SD
mod6 <- nlme(LENGTH ~ Linf*(1-exp(-k0*(AgeFraction-(t0))))  ,fixed= Linf+k0+t0 ~1 ,random=  list(Area = Linf+k0~1,Cohort = Linf+k0~1) ,data=data_fin2, start=c(Linf=259.7278,k0=0.1492918,t0=-3.639129)) 
summary(mod6); AIC(mod6)
summarise_model(mod6, ci = FALSE)
p3 <-plot(mod6, Area ~ resid(.), abline = c(0,0));p3
p4 <-plot(mod6, as.factor(Cohort) ~ resid(.), abline = c(0,0));p4
p5 <-plot(mod6, as.factor(Year) ~ resid(.), abline = c(0,0));p5

## Year in SD & Cohort in SD
mod7 <- nlme(
  LENGTH ~ Linf * (1 - exp(-k0 * (AgeFraction - t0))),
  data  = data_fin2,
  fixed = Linf + k0 + t0 ~ 1,
  random = list(
    YearArea = Linf + k0 ~ 1 ,
    CohortArea = Linf + k0 ~ 1 
  ),
  start = c(Linf = 259.7278, k0 = 0.1492918, t0 = -3.639129)
)
summary(mod7); AIC(mod7)
summarise_model(mod7, ci = FALSE)
p3 <-plot(mod7, Area ~ resid(.), abline = c(0,0));p3
p4 <-plot(mod7, as.factor(Cohort) ~ resid(.), abline = c(0,0));p4
p5 <-plot(mod7, as.factor(Year) ~ resid(.), abline = c(0,0));p5
# residuals by genre
qqnorm(mod7, ~resid(.) | Cohort, pch = 20, col = "black" )
# residuals by genre
qqnorm(mod7, ~resid(.) | Area, pch = 20, col = "black" )

#compare AIC VB
AIC(mod1); AIC(mod2);AIC(mod3);AIC(mod4);AIC(mod5);AIC(mod6)

#######################
# BEST MODEL : mod4   # should be mod6 now!
find_typical(mod4, probs = c(.25, .50, .75))
extract_fixed_effects(mod4)
rf <- extract_random_effects(mod4)
ref <-extract_random_coefs(mod4  )
#library(stringr)
db_area <- ref %>% dplyr::filter(group_var == "Area")
db_areayr <- ref %>% dplyr::filter(group_var == "Year")
db_areayr$Year=(str_sub(db_areayr$group, start=-4))
db_areayr$Area =(str_sub(db_areayr$group, end = -6  )) 

# VB by year (area) 
ggplot(db_areayr  , aes(as.numeric(Year) ,value)) + geom_point(  aes(as.numeric(Year),value, colour= Area)) +geom_smooth()+facet_wrap(~effect, scales = "free")
ggplot(db_areayr  , aes(as.numeric(Year) ,value)) + geom_boxplot(aes(group =Year )) +geom_smooth()+facet_wrap(~effect, scales = "free")

# biplot by area
db_area_wide <- db_area %>%
  pivot_wider(
    names_from = effect,
    values_from = c(value, lower_2.5, upper_97.5),
    names_glue = "{effect}_{.value}"
  )
ggplot(db_area_wide, aes(x = Linf_value, y = k0_value, label = group, color = group)) +
  geom_point(size = 3) +
  # geom_text(size = 3,check_overlap = TRUE) +
  geom_text_repel(aes(label = group), size = 3, vjust=-0.5)+
  #geom_text(vjust = -0.5, hjust = 0.3) +
  labs(x = expression(L[infinity]), y = "k") +
  theme_bw() +
  theme(legend.position = "none")+ geom_vline(xintercept = 230.01, linetype = "dotted", color = "black")+ geom_hline(yintercept = 0.30, linetype = "dotted", color = "black")

# Linf area*year 
jpeg("timeseries_Linf_bySD.jpeg",width = 400, height = 200, units = "mm", res = 600)
ggplot(db_areayr %>% dplyr::filter(effect == "Linf") , aes(as.numeric(Year) ,value, colour= Area)) + geom_point(  aes(as.numeric(Year),value), size=2)+geom_smooth(method ="lm")+facet_wrap(~Area) + theme_bw() +theme(legend.position="NULL")+ ylab(expression(L[infinity]))+ xlab("Year")
dev.off()
# k area*year 
jpeg("timeseries_k_bySD.jpeg",width = 400, height = 200, units = "mm", res = 600)
ggplot(db_areayr %>% dplyr::filter(effect == "k0") , aes(as.numeric(Year) ,value, colour= Area)) + geom_point(  aes(as.numeric(Year),value), size=2)+geom_smooth(method ="lm")+facet_wrap(~Area)+theme_bw()+theme(legend.position="NULL")+ ylab("k")+ xlab("Year")
dev.off()


################ PREDICTION VB 
# create prediction grid
Year <- unique(data_fin2$Year)
Area <- unique(data_fin2$Area)
AgeFraction <- unique(data_fin2$AgeFraction)
ll <- list(Year=Year,Area=Area,AgeFraction=AgeFraction)
newdbpred <- expand.grid(ll)

predf <-as.data.frame(predict(mod4,newdbpred, type = "response")); predf
newdbpred$pred  <- predf$`predict(mod4, newdbpred, type = "response")`
findb <-left_join(data_fin2,newdbpred, by= c("Year","AgeFraction", "Area") ) 
tmp_area_ageyr_pred <- findb %>%
  group_by(Year,Area,agePlus ) %>%  
  dplyr::summarise(mean_pred=round(mean(pred),0), 
                   mean_len=round(mean(LENGTH),0)) %>%
  as.data.frame()
# plot RES by SD and Year
tmp_area_ageyr_pred$RES <-(tmp_area_ageyr_pred$mean_len- tmp_area_ageyr_pred$mean_pred)/tmp_area_ageyr_pred$mean_len * 100
jpeg("pred_discrep_VB.jpeg",width = 300, height = 200, units = "mm", res = 600)
ggplot(tmp_area_ageyr_pred, aes(Year, RES,colour=Area ))+ facet_wrap(~agePlus)+ theme_bw()  + geom_line()  +theme(legend.position="left") + ylab("Pred discrepancy (%)") + geom_hline(yintercept = 0, linetype="dotted")
dev.off()


# pred for GAM part
Year <- unique(data_fin2$Year)
Area <- unique(data_fin2$Area)
AgeFraction <- unique(data_fin2$agePlus+0.5) # using the middle year age
ll <- list(Year=Year,Area=Area,AgeFraction=AgeFraction)
newdbpred <- expand.grid(ll)
predf <-as.data.frame(predict(mod4,newdbpred, type = "response")); predf
newdbpred$pred  <- predf$`predict(mod4, newdbpred, type = "response")`
newdbpred$Age <- floor(newdbpred$Age)
newdbpred$Area <-  factor(newdbpred$Area, levels = c("27.3.d.25", "27.3.d.26", "27.3.d.27", "27.3.d.28.2","27.3.d.29","27.3.d.32"))
ggplot(newdbpred, aes(Year, pred,colour=as.factor(Age) ))+ geom_point()+ facet_wrap(~Area)+ theme_bw()  + geom_line(linetype= "dashed")  +theme(legend.position="left") 
newdbpred  <- newdbpred  %>%
  mutate(SD = recode(Area,
                     "27.3.d.25" = "SD25",
                     "27.3.d.26" = "SD26",
                     "27.3.d.27" = "SD27",
                     "27.3.d.28.2" = "SD28",
                     "27.3.d.29" = "SD29",
                     "27.3.d.32" = "SD32"))

predwider <-newdbpred %>% select(Year,SD,Age,pred) %>% pivot_wider(names_from = Age , values_from = pred) 
predwider <- predwider %>% arrange(Year, SD)
write.csv2(predwider, "pred_all_ages.csv")

age2to6 <- newdbpred %>% dplyr::filter(Age < 7 , Age > 1 )  %>%
  group_by(Year,SD  ) %>%  
  dplyr::summarise(mean_predTL=mean(pred)) %>%
  as.data.frame()

ggplot(age2to6, aes(Year, mean_predTL ))+ geom_point(size=2)+ facet_wrap(~SD)+ theme_bw()  + geom_line()  +theme(legend.position="bottom") 


#######################################################################################################
###### L-W #############################
########################################
# Log-transform variables
data_fin2 <- data_fin2 %>%
  mutate(logW = log(weight),
         logL = log(LENGTH))

# no random effect
LW_MODEL1 <-nls(logW ~ loga + b * logL , start=list(loga=log(0.0054),b=3),
          data=(data_fin2  ),lower=c(log(0.0000054), 0),
          upper=c(log(1),5),algorithm="port")
summary(LW_MODEL1); AIC(LW_MODEL1)
p3 <-plot(LW_MODEL1, Area ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p3
p4 <-plot(LW_MODEL1, as.factor(Year) ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p4

# only area
LW_MODEL2 <- nlme(
  logW ~ loga + b * logL,
  fixed = loga + b ~ 1,
  random = list(Area = loga + b ~ 1),
  data = data_fin2,
  start = c(loga = log(0.0054), b = 3)  # reasonable starting values
)
summary(LW_MODEL2); AIC(LW_MODEL2)
p3 <-plot(LW_MODEL2, Area ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p3
p4 <-plot(LW_MODEL2, as.factor(Year) ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p4

# only year
LW_MODEL3 <- nlme(
  logW ~ loga + b * logL,
  fixed = loga + b ~ 1,
  random = list(Year = loga + b ~ 1),
  data = data_fin2,
  start = c(loga = log(0.0054), b = 3)  # reasonable starting values
)
summary(LW_MODEL3); AIC(LW_MODEL3)
p3 <-plot(LW_MODEL3, Area ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p3
p4 <-plot(LW_MODEL3, as.factor(Year) ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p4

# ## by Year in SD
LW_MODEL4 <- nlme(
  logW ~ loga + b * logL,
  fixed = loga + b ~ 1,
  random = list(Area = loga + b ~ 1, Year = loga + b ~ 1),
  data = data_fin2,
  start = c(loga = log(0.0054), b = 3)  # reasonable starting values
)
summary(LW_MODEL4); AIC(LW_MODEL4)
p3 <-plot(LW_MODEL4, Area ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p3
p4 <-plot(LW_MODEL4, as.factor(Year) ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p4
p5 <-plot(LW_MODEL4, as.factor(Cohort) ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p5
summarise_model(LW_MODEL4, ci = FALSE)
exp(-12.45)

# ## by cohort
LW_MODEL5 <- nlme(
  logW ~ loga + b * logL,
  fixed = loga + b ~ 1,
  random = list( Cohort = loga + b ~ 1),
  data = data_fin2,
  start = c(loga = log(0.0054), b = 3)  # reasonable starting values
)
summary(LW_MODEL5); AIC(LW_MODEL5)
p3 <-plot(LW_MODEL5, Area ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p3
p4 <-plot(LW_MODEL5, as.factor(Year) ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p4

# ## by Cohort in SD
LW_MODEL6 <- nlme(
  logW ~ loga + b * logL,
  fixed = loga + b ~ 1,
  random = list(Area = loga + b ~ 1, Cohort = loga + b ~ 1),
  data = data_fin2,
  start = c(loga = log(0.0054), b = 3)  # reasonable starting values
)
summary(LW_MODEL6); AIC(LW_MODEL6)
p3 <-plot(LW_MODEL6, Area ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p3
p4 <-plot(LW_MODEL6, as.factor(Year) ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p4
p5 <-plot(LW_MODEL6, as.factor(Cohort) ~ resid(.), abline = c(0,0),xlim = c(-1.5, 1.5));p5
summarise_model(LW_MODEL6, ci = FALSE)

#compare AIC LW
AIC(LW_MODEL1); AIC(LW_MODEL2);AIC(LW_MODEL3);AIC(LW_MODEL4);AIC(LW_MODEL5);AIC(LW_MODEL6)

# BEST MODEL: mod4
rf <- extract_random_effects(LW_MODEL4)
ref <-extract_random_coefs(LW_MODEL4  )
db_area <- ref %>% dplyr::filter(group_var == "Area")
db_areayr <- ref %>% dplyr::filter(group_var == "Year")
db_areayr$Year=(str_sub(db_areayr$group, start=-4))
db_areayr$Area =(str_sub(db_areayr$group, end = -6  )) 

ggplot(db_areayr %>% dplyr::filter(effect == "loga") , aes(as.numeric(Year) ,exp(value), colour= Area)) + geom_point(  aes(as.numeric(Year),exp(value)))+geom_smooth(method ="lm")+facet_wrap(~Area)+ theme_bw()  +theme(legend.position="")+ ylab("a")+ xlab("Year")

ggplot(db_areayr %>% dplyr::filter(effect == "b") , aes(as.numeric(Year) ,value, colour= Area)) + geom_point(  aes(as.numeric(Year),value))+geom_smooth(method ="lm")+facet_wrap(~Area)+ theme_bw()+theme(legend.position="NULL")+ ylab("b")+ xlab("Year")

# biplot by area
db_area_wide <- db_area %>%
  pivot_wider(
    names_from = effect,
    values_from = c(value, lower_2.5, upper_97.5),
    names_glue = "{effect}_{.value}"
  )
ggplot(db_area_wide, aes(x = b_value, y = exp(loga_value), label = group, color = group)) +
  geom_point(size = 3) +
 # geom_text(size = 3,check_overlap = TRUE) +
   geom_text_repel(aes(label = group), size = 3, vjust=-0.5)+
  #geom_text(vjust = -0.5, hjust = 0.3) +
  labs(x = "b", y = "a") +
  theme_bw() +
  theme(legend.position = "none")+ geom_vline(xintercept = 3.082203, linetype = "dotted", color = "black")+ geom_hline(yintercept = exp(-12.411302), linetype = "dotted", color = "black")

################ PREDICTION LW 
# create prediction grid
Year <- unique(data_fin2$Year)
Area <- unique(data_fin2$Area)
logL <- unique(data_fin2$logL)
ll <- list(Year=Year,Area=Area,logL=logL)
newdbpred <- expand.grid(ll)
predf <-as.data.frame(predict(LW_MODEL4,newdbpred, type = "response")); predf
newdbpred$pred  <- predf$`predict(LW_MODEL4, newdbpred, type = "response")`
findb <-left_join(data_fin2,newdbpred, by= c("Year","logL", "Area") ) 
tmp_area_ageyr_pred <- findb %>%
  group_by(Year,Area, agePlus ) %>%  
  dplyr::summarise(mean_pred=(mean(exp(pred))), 
                   mean_w=(mean(weight))) %>%
  as.data.frame()
# plot RES by SD and Year
tmp_area_ageyr_pred$RES <- (tmp_area_ageyr_pred$mean_w- tmp_area_ageyr_pred$mean_pred)/tmp_area_ageyr_pred$mean_w * 100
jpeg("pred_discrep_meanW_YinA.jpeg",width = 300, height = 200, units = "mm", res = 600)
ggplot(tmp_area_ageyr_pred, aes(Year, RES,colour=Area ))+ facet_wrap(~agePlus)+ theme_bw()  + geom_line()  +theme(legend.position="left") + ylab("Pred discrepancy (%)") + geom_hline(yintercept = 0, linetype="dotted")
dev.off()

#########################################################
# pred per GAM (via model prediction)
Year <- unique(data_fin2$Year)
Area <- unique(data_fin2$Area)
logL <- log(165)
ll <- list(Year=Year,Area=Area,logL=logL)
newdbpred <- expand.grid(ll)
predf <-as.data.frame(predict(LW_MODEL4,newdbpred, type = "response")); predf
newdbpred$predCond  <- exp(predf$`predict(LW_MODEL4, newdbpred, type = "response")`)

condition_df  <- newdbpred  %>%
  mutate(SD = recode(Area,
                     "27.3.d.25" = "SD25",
                     "27.3.d.26" = "SD26",
                     "27.3.d.27" = "SD27",
                     "27.3.d.28.2" = "SD28",
                     "27.3.d.29" = "SD29",
                     "27.3.d.32" = "SD32")) %>% dplyr::select(-Area,-logL)
condition_df$SD <-  factor(condition_df$SD, levels = c("SD25", "SD26", "SD27", "SD28","SD29","SD32"))
ggplot(condition_df, aes(Year, predCond,colour=SD ))+ geom_point()+ facet_wrap(~SD)+ theme_bw()  + geom_line()  +theme(legend.position="left") 
# write.csv2(condition_df, "pred_conditionL165.csv")




######################################################
# final plot with proxy for GAMs 
proxydb <-left_join(age2to6,condition_df)
write.csv2(proxydb, "pred_GAM.csv")

jpeg("proxy_bySD.jpeg",width = 400, height = 200, units = "mm", res = 600)
ggplot(proxydb, aes(x = Year,colour=SD)) +
  geom_line(aes(y = mean_predTL)) +
  geom_line(aes(y = predCond * 8), linetype= "dashed") +  # scale for overlay
  geom_point(aes(y = mean_predTL)) +
  geom_point(aes(y = predCond * 8), shape = 1) + 
  scale_y_continuous(
    name = "Length-at-age (mm)",
    sec.axis = sec_axis(~ . /8, name = "Condition (gr)")
  )+ facet_wrap(~SD)+ theme_bw()   +theme(legend.position="") 
dev.off()

