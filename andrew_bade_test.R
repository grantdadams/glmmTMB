# devtools::install_github("grantdadams/glmmTMB/glmmTMB")
# install.packages("pscl")
# Updated by Grant Adams: adamsgd@uw.edu

library(pscl)
library(glmmTMB)

# Fit using glmmTMB update
data("bioChemists", package = "pscl")
zinb <- glmmTMB(art ~ fem + mar + kid5 + phd + ment, ziformula = ~ ., data =
                  bioChemists, family = nbinom2(link = "log"))
pred <- predict(zinb, type="like") # Get likelihood of each obs

# Fit using pscl
fm_zinb2 <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin")
preds2 <- predict(fm_zinb2, type="prob")
pred2 <- c(preds2[seq_along(bioChemists$art), bioChemists$art]) # Extract probs

# pred is about equal to pred2. Likely differences due to optimizer and maybe distributions

# What you want to do steps:
# 1. Create new data frame by copying each row for the number of unique values of the response
y_vals <- min(bioChemists$art, na.rm = T):max(bioChemists$art, na.rm = T)
new_data <- bioChemists[rep(1:nrow(bioChemists),each=length(y_vals)),]
new_data$art = rep(y_vals, nrow(bioChemists))

# 2. Get prob for each new obs
pred3 <- predict(zinb, newdata = new_data, type="like") # Get likelihood of each obs
new_data$pred <- pred3

# 3. Make into matrix
preds3 <- matrix(pred3, ncol =  length(y_vals), byrow = TRUE)
colnames(preds3) = y_vals

# TEST
# Compare with https://stackoverflow.com/questions/22314921/no-zeros-predicted-from-zeroinfl-object-in-r
# Most pubs
most_pubs <- max(bioChemists$art)
most_pubs
extreme_biochemist <- bioChemists$art==most_pubs
which(extreme_biochemist)

par(mfrow = c(1,2))

# pscl
expected <- predict(fm_zinb2, type="response")[extreme_biochemist]
# barplot returns the midpoints for counts 0 up to 19
midpoints<-barplot(preds2[extreme_biochemist,],
                   xlab="Predicted #pubs", ylab="Relative chance among biochemists", main = "pscl")
# add 1 because the first count is 0
abline(v=midpoints[19+1],col="red",lwd=3)
abline(v=midpoints[round(expected)+1],col="yellow",lwd=3)

# glmmTMB
expected <- predict(zinb, type="response")[extreme_biochemist]
# barplot returns the midpoints for counts 0 up to 19
midpoints<-barplot(preds3[extreme_biochemist,],
                   xlab="Predicted #pubs", ylab="Relative chance among biochemists", main = "glmmTMB")
# add 1 because the first count is 0
abline(v=midpoints[19+1],col="red",lwd=3)
abline(v=midpoints[round(expected)+1],col="yellow",lwd=3)

