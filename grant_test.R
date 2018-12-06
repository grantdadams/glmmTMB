library(pscl)
library(glmmTMB)
data("bioChemists", package = "pscl")
zinb <- glmmTMB(art ~ fem + mar + kid5 + phd + ment, ziformula = ~ ., data =
                  bioChemists, family = nbinom2(link = "log"))

pred <- predict(zinb, type="like")


