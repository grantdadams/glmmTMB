library(pscl)
library(glmmTMB)

# Fit using glmmTMB update
data("bioChemists", package = "pscl")
zinb <- glmmTMB(art ~ fem + mar + kid5 + phd + ment, ziformula = ~ ., data =
                  bioChemists, family = nbinom2(link = "log"))
pred <- predict(zinb, type="like") # Get likelihood of each obs

# Fit using pscl
fm_zinb2 <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")
preds <- predict(fm_zinb2, type="prob")
pred2 <- c(preds[seq_along(bioChemists$art), bioChemists$art]) # Extract probs

# pred is about equal to pred2. Likely differences due to optimizer and maybe distributions

# What you want to do steps:
# 1. Create new data frame by copying each row for the number of unique values of the response
y_vals <- min(bioChemists$art, na.rm = T):max(bioChemists$art, na.rm = T)
new_data <- bioChemists[rep(1:nrow(bioChemists),each=length(y_vals)),]
new_data$art = rep(y_vals, nrow(bioChemists))

# 2. Get prob for each new obs
pred <- predict(zinb, newdata = bioChemists, type="like") # Get likelihood of each obs

