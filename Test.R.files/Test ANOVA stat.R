invis <- read.csv("inv.simpson.index2023.01.01.csv")
invis

one.way <- aov(shannon_div ~ group, data = invis)
summary(one.way)
par(mfrow = c(1, 2)) # combine plots
library(car)
hist(one.way$residuals)
qqPlot(one.way$residuals,
       id = FALSE # id = FALSE to remove point identification
)
shapiro.test(one.way$residuals)

boxplot(shannon_div ~ group, data = invis)
aggregate(shannon_div ~ group,
          data = invis,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)


res_aov <- aov(shannon_div ~ group,
            data = invis
)

res_aov
summary(res_aov)
library("report") # Load the package every time you start R

report(res_aov)
TukeyHSD(aov(formula = shannon_div ~ group, data = invis))
library(multcomp)

# Tukey HSD test:
post_test <- glht(res_aov,
                  linfct = mcp(species = "Tukey")
)

summary(post_test)
