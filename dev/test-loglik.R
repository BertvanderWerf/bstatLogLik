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

# create a list with parameters for each level of f (as defined in formulas in create_paramter_info())
pars <- tapply(data$x, data$f, weibull_init)
pars <- as.list(as.data.frame(do.call(rbind, pars)))

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
compare_models(fit4, fit2, fit3_lambda, fit3_k)

vcov(fit4)
estfun(fit4)
vcov_sandwich(fit4)

coef(fit4)

s<-summary(fit4$parameter_info)
print(parameter_info(fit4))

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


# -------------------------------------------------------------------------


new_x <- seq(min(data$x), max(data$x), length.out = 50)
new_data <- expand.grid(x = new_x, f = factor(c("a", "b")))

p_survival <- predict(fit4,
                      se.fit = TRUE,
                      nr_simulations = 250,
                      newdata = new_data,
                      newfun = weibull_cdf)

# Plot
plot_surv <- cbind(new_data, p_survival)

new_data$y <- evaluate_function(weibull_cdf, fit4$vars, fit4$parameter_info, data=new_data)

ggplot() + geom_point(aes(x=x, y=y, fill=f), data=new_data)

ggplot() +
  geom_line(aes(x = x, y = fitted, col = f), data = plot_surv, linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = lower.95, ymax = upper.95, fill = f),
              data = plot_surv, alpha = 0.2) +
  labs(title = "Weibull Survival Function",
       x = "x", y = "Survival Probability",
       color = "Group", fill = "Group") +
  theme_minimal()


vcov(fit)

logLik(fit)
AIC(fit)
BIC(fit)
predict(fit)
trace_history(fit4)
control_settings(fit)
formula(fit)
coef(fit)
converged(fit)
terms(fit)
vcov_sandwich(fit)

as.data.frame(param_info)

methods(class='loglik')
# [1] AIC              anova            BIC              coef             control_settings converged
# [7] df.residual      fitted           formula          gradient_matrix  hessian          logLik
# [13] lrt              model.frame      nobs             parameter_info   predict          print
# [19] summary          terms            trace_history    vcov
methods(class='parameter_info')


object <- create_parameter_info(param_info, formulas=list(lambda=~0+f:x), model_data=data)
orthogonal_groups.parameter_info <- function(object) {
  L <- list()
  all_equal <- TRUE
  k <- 1
  for (i in seq_len(length(object$parameter_values))) {
    x <- object$design_matrices[[i]][,as.logical(param_info$estimate_flag[[i]]),
                                    drop=FALSE]
    if (length(x)>0) {
      L[[names(object$parameter_values)[i]]] <- (t(x) %*% x)>0
      n <- length(L)
      if (n>1) {
        all_equal <- equal(L[[n-1]], L[[n]])
      }
    }
  }

  if (all_equal==TRUE & length(L[[1]])>1) {
    # check if petrie matrix
    # sort matrix in order of number of 0 e.g. 1, 0, 0, 0 first
    if (equal(L[[1]], diag(L[[1]]))) { # groups found

    } else {
      n_zeros <- colSums(L[[1]]==0)
      n_zeros_unique <- unique(n_zeros)
      if (length(nZeros)==1 && n_zeros==0) { # no groups found
      } else {  # groups found by blocks
      }
    }


    groups <- colnames(L[[1]])
    return(groups) # to do
  }
  return(NULL)
}

rstudioapi::documentSaveAll()
.rs.restartR()
