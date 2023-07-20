## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
# install.packages("ipolygrowth")
library(ipolygrowth)
library(dplyr)
library(ggplot2)
library(kableExtra) # for RMarkdown table output

## -----------------------------------------------------------------------------
df.gr <- growthrates::bactgrowth
str(df.gr)

## -----------------------------------------------------------------------------
head(df.gr)

## ---- fig.height=10, fig.width=13---------------------------------------------
ggplot()+
  geom_point(data = df.gr, aes(x = time, y = value, color = factor(replicate)))+
  facet_wrap(~ strain + conc)+
  labs(color = "replicate")+
  theme_bw()

## -----------------------------------------------------------------------------
df.singlesample <- df.gr %>% filter(strain == "D", conc == 15.63) %>% select(-strain, -conc)
str(df.singlesample)
table(df.singlesample$replicate, df.singlesample$time)

## ---- results='asis'----------------------------------------------------------
out.singlesample <- ipg_singlesample(data = df.singlesample, time.name = "time", y.name = "value")

out.singlesample$estimates %>%
  kable %>% kable_styling("striped", full_width = F)  # table formatting for rmarkdown

## ---- fig.height=3, fig.width=4-----------------------------------------------
ggplot()+
  geom_point(data = df.singlesample, aes(x = time, y = value, color = factor(replicate)))+ 
  geom_line(data = out.singlesample$fitted, aes(x = time, y = fit))+ 
  labs(color = "replicate")+
  scale_x_continuous(n.breaks = 10)+
  scale_y_continuous(n.breaks = 7)+
  theme_bw()

## ---- results='asis'----------------------------------------------------------
out.singlesample2 <- ipg_singlesample(data = df.singlesample, time.name = "time", y.name = "value", epsilon = 1/100)
out.singlesample2$estimates %>%
  kable %>% kable_styling("striped", full_width = F)  # table formatting for rmarkdown

## -----------------------------------------------------------------------------
df.gr2 <- df.gr %>%  mutate(id = paste(strain, conc, sep = "-"))
str(df.gr2)
table(df.gr2$strain, df.gr2$conc)
unique(df.gr2$id)

## ---- results='asis'----------------------------------------------------------
out.multi.f <- ipg_multisample(data = df.gr2, id = "id", time.name = "time", y.name = "value", epsilon = 0.2/100)
out.multi.f$estimates %>%
  kable() %>% kable_styling("striped", full_width = F) %>% scroll_box(width = "800px", height = "300px")  # table formatting for rmarkdown

## ---- fig.height=10, fig.width=13---------------------------------------------
ggplot()+
  geom_point(data = df.gr2, aes(x = time, y = value, color = factor(replicate)))+ 
  geom_line(data = out.multi.f$fitted, aes(x = time, y = fit))+ 
  facet_wrap(~ id)+
  labs(color = "replicate")+
  theme_bw()

