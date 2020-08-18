csfit1 <- function (cc, G, logRm = FALSE, logBase = 2)
{
     if (logRm == TRUE) {
          G = logBase^G
     }
     fit1 = lsfit(cc, G, intercept = FALSE)
     #se1 = ls.diag(fit1)$std.err
     if (logRm == TRUE) {
          ghat = log(fit1$coefficients, logBase)
          ghat[is.nan(ghat)] = 0
          #se = log(se1, logBase)
          return(list(ghat = ghat, residuals = fit1$residuals))
     }
     else {
          return(list(ghat = fit1$coefficients, residuals = fit1$residuals))
     }
}
