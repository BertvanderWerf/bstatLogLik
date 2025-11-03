rstudioapi::documentSaveAll()
devtools::document()
devtools::load_all()

# Simulate some data
set.seed(1)
data <- data.frame(
  x = c(sort(rweibull(100, shape = 2, scale = 5)),
        sort(rweibull(100, shape = 3, scale = 3))),
  f = factor(rep(c("a","b"), c(100, 100))),
  cy = c(seq(1, 100)/100, seq(1, 100)/100)
)

pars <- as.list(as.data.frame(t(as.matrix(tapply(data$x, data$f, weibull_init)))))

# Prepare parameter info using the distribution's init function
param_info <- create_parameter_info(
  parameter_values = pars,
  estimate_flag = c(k = 1, lambda = 1, c = 0),
  formulas = list(k = ~0 + f, lambda = ~0+f, c = ~1),
#  link=list(lambda='log'),
  model_data = data
)
#param_info$parameter_values$lambda <- log(param_info$parameter_values$lambda)
summary(param_info)

# Fit model
fit4 <- loglik(
  fun = weibull_pdf,
  fun_vars=c(x = "x"),
  parameter_info = param_info,
  data = data
)

fit3_k <- loglik(
  fun = weibull_pdf,
  fun_vars=c(x = "x"),
  parameter_info = create_parameter_info(param_info, formula=list(lambda=~1)),
  data = data
)

fit3_lambda <- loglik(
  fun = weibull_pdf,
  fun_vars=c(x = "x"),
  parameter_info = create_parameter_info(param_info, formula=list(k=~1)),
  data = data
)

fit2 <- loglik(
  fun = weibull_pdf,
  fun_vars=c(x = "x"),
  parameter_info = create_parameter_info(param_info, formula=list(k=~1, lambda=~1)),
  data = data
)
lrt(fit4, fit2, fit3_lambda, fit3_k)


vcov(fit)
estfun(fit)
sandwich(fit)

parameter_info <- fit$parameter_info

coef(fit)

summary(parameter_info)

# Inspect result
print(fit4)
s=summary(fit4)


# Plot fitted values, diagnostics, etc.
library(ggplot2)
pars <- as.data.frame(parameter_info(fit4))
ggplot() + geom_point(aes(x=x, y=cy, col=f), data=data) +
  geom_line(aes(x=x, weibull_cdf(x, pars$k, pars$lambda, pars$c), col=f), data=data)

p1<-predict(fit4, se.fit = TRUE, nr_simulations=250, newfun=weibull_cdf)
ggplot() + geom_point(aes(x=x, y=cy, col=f), data=data) +
  geom_line(aes(x=x, weibull_cdf(x, pars$k, pars$lambda, pars$c), col=f), data=data)+
  geom_line(aes(x=x, y=fitted, col=f), data=cbind(data, p1)) +
  geom_ribbon(aes(x=x, ymin=lower.95, ymax=upper.95, fill=f), data=cbind(data,p1), alpha=0.2)

p<-predict(fit4, se.fit = TRUE, nr_simulations=250)
ggplot() + geom_line(aes(x=x, y=fitted, col=f), data=cbind(data, p)) +
  geom_ribbon(aes(x=x, ymin=lower.95, ymax=upper.95, fill=f), data=cbind(data,p), alpha=0.2)

ph <- predict(fit4, se.fit = TRUE, nr_simulations=250, newfun=weibull_hazard)
ggplot() + geom_line(aes(x=x, y=fitted, col=f), data=cbind(data, ph)) +
  geom_ribbon(aes(x=x, ymin=lower.95, ymax=upper.95, fill=f), data=cbind(data,ph), alpha=0.2)

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
sandwich(fit)

methods(class='loglik')
# [1] AIC              anova            BIC              coef             control_settings converged
# [7] df.residual      fitted           formula          gradient_matrix  hessian          logLik
# [13] lrt              model.frame      nobs             parameter_info   predict          print
# [19] summary          terms            trace_history    vcov

rstudioapi::documentSaveAll()
.rs.restartR()
