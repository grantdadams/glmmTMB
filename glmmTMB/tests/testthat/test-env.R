## check that everything works in weird environments

stopifnot(require("testthat"),
          require("glmmTMB"))

data(sleepstudy, cbpp,
     package = "lme4")

## need global env for test_that
sleepstudy <<- transform(sleepstudy, DaysFac = factor(Days))

context("basic examples")

test_that("basic example #1", {
    fitFun <- function(dat){
        glmmTMB(Reaction ~ Days + (1|Subject), data=dat)
    }
    f0 <- glmmTMB(Reaction ~ Days + (1|Subject), data=sleepstudy)
    f1 <- fitFun(sleepstudy)
    uncall <- function(x) {
        x$call <- NULL
        return(x)
    }
    expect_equal(uncall(f0),uncall(f1))
})