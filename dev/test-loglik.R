rstudioapi::documentSaveAll()
devtools::load_all()

# Simulate some data
set.seed(1)
data <- data.frame(
  x = sort(rweibull(100, shape = 2, scale = 5))
)
data$cy <- seq(1, nrow(data))/nrow(data)

# Prepare parameter info using the distribution's init function
param_info <- create_parameter_info(
  parameter_values = weibull_init(data$x, verbose=TRUE),
  estimate_flag = c(k = 1, lambda = 1, c = 0),
  formulas = list(k = ~1, lambda = ~1, c = ~1),
#  link=list(lambda='log'),
  model_data = data
)
#param_info$parameter_values$lambda <- log(param_info$parameter_values$lambda)
summary(param_info)

# Fit model
fit <- loglik(
  fun = weibull_pdf,
  fun_vars=c(x = "x"),
  parameter_info = param_info,
  data = data
)
parameter_info <- fit$parameter_info

summary(parameter_info)

# Inspect result
print(fit)
summary(fit)


# Plot fitted values, diagnostics, etc.
library(ggplot2)
ggplot() + geom_point(aes(x=x, y=cy), data=data) +
  geom_line(aes(x=seq(0,11,0.1), weibull_cdf(seq(0,11,0.1),
                                            parameter_info$parameter_values$k,
                                            parameter_info$parameter_values$lambda,
                                            parameter_info$parameter_values$c)))

vcov(fit)
logLik(fit)
AIC(fit)
BIC(fit)
predict(fit)
trace_history(fit)
control_settings(fit)
formula(fit)
coef(fit)
converged(fit)
terms(fit)
vcov_sandwich.loglik(fit)
sandwich::sandwich(fit)

estfun(fit)

methods(class='loglik')
# [1] AIC              anova            BIC              coef             control_settings converged
# [7] df.residual      fitted           formula          gradient_matrix  hessian          logLik
# [13] lrt              model.frame      nobs             parameter_info   predict          print
# [19] summary          terms            trace_history    vcov

rstudioapi::documentSaveAll()
.rs.restartR()
