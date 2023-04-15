install.packages("palmerpenguins")
library(palmerpenguins)

library(tidyverse)
peg <- penguins
head(peg)

aov(flipper_length_mm ~ species,
    data = peg)

oneway.test(flipper_length_mm ~ species,
    data = peg)

res_aov <- aov(flipper_length_mm ~ species,
               data = peg)
res_aov

hist(res_aov$residuals)

library(car)
qqPlot(res_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)
shapiro.test(res_aov$residuals)

boxplot(flipper_length_mm ~ species,
        data = peg
)

aggregate(flipper_length_mm ~ species,
          data = peg,
          function(x) round(c(mean = mean(x), sd = sd(x)), 2)
)
