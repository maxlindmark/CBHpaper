# Test difference in predicted length-at-age for juveniles and adults
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_light())
library(stringr)

d <- read.csv("pred_all_ages.csv", sep = ";", dec = ",")

# Scale lengths
d2 <- d |> 
  pivot_longer(X0:X9, names_to = "age", values_to = "length") |> 
  drop_na(length) |> 
  mutate(age = str_remove(age, "X"),
         age = as.numeric(age)) |> 
  mutate(l_sc = as.numeric(scale(length)),
         .by = c(SD, age)) |> 
  mutate(life_stage = ifelse(age <= 2, "adult", "juvenile")) |> 
  filter(age > 0 & age <= 6)
         
ggplot(d2, aes(Year, l_sc,
               color = as.factor(life_stage),
               linetype = as.factor(age))) + 
  geom_line() +
  labs(color = "") +
  facet_wrap(~SD) + 
  theme(legend.position = "bottom")

# How about averaging z-scored values
d2 |> 
  summarise(l_sc = mean(l_sc), .by = c(Year, life_stage, SD)) |> 
  ggplot(aes(Year, l_sc,
             color = as.factor(life_stage))) + 
  geom_line() +
  labs(color = "") +
  facet_wrap(~SD) + 
  theme(legend.position = "bottom")

# Is there a difference in z-scored trends between juveniles (0-2) and adults 
# Ignore SD
dm <- d2 |> 
  summarise(l_sc = mean(l_sc), .by = c(Year, life_stage, SD))

m <- lm(l_sc ~ Year*life_stage, data = dm)
summary(m)

# Super significant... 
# Predict?
nd <- expand.grid(Year = seq(min(dm$Year), max(dm$Year)),
                  life_stage = c("juvenile", "adult"))

nd$est <- as.numeric(predict(m, newdata = nd))

dm |> 
  ggplot(aes(Year, l_sc,
             color = as.factor(life_stage))) + 
  geom_point() +
  geom_line(data = nd, aes(Year, est)) +
  labs(color = "") +
  theme(legend.position = "bottom")

