library(lme4)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(contrasts=c("contr.sum","contr.poly"))
source("getStanDataInfo.R")
source("getSummary.R")
source("printSummary.R")

perception_data <- read.csv("perception_study_data.csv",header=TRUE)
perception_data$SYLL_MAXF0 <- scale(perception_data$SYLL_MAXF0)[,1]
perception_data$SYLL_MEAN_INT <- scale(perception_data$SYLL_MEAN_INT)[,1]
perception_data$SYLL_DUR <- scale(perception_data$SYLL_DUR)[,1]

max_form <- USER_RESP ~ PRIMARY_STRING + FUNCTION + EXPERIMENT +
	ONEBACK + SYLL_MAXF0 * SYLL_MEAN_INT * SYLL_DUR +
	(1 + PRIMARY_STRING + FUNCTION + ONEBACK +
	SYLL_MAXF0 * SYLL_MEAN_INT * SYLL_DUR | SUBJECT) +
	(1 + EXPERIMENT + ONEBACK | ITEM)

	
##### lme4
Full_Model_lme4 <- glmer(
	formula = max_form, data = perception_data, family = binomial)
# Model failed to converge with max|grad| = 0.0618226 (tol = 0.001, component 1)
# Time: 26 hours

Full_Model.bobyqa <- glmer(
	formula = max_form, data = perception_data, family = binomial,
	glmerControl(optimizer="bobyqa"))
# Model failed to converge with max|grad| = 0.00561514 (tol = 0.001, component 1)

Full_Model.Nelder_Mead <- glmer(
	formula = max_form, data = perception_data, family = binomial,
	glmerControl(optimizer="Nelder_Mead"))
# Model failed to converge: degenerate  Hessian with 2 negative eigenvalues

Intercepts_Only <- glmer(
	USER_RESP ~ PRIMARY_STRING + FUNCTION + EXPERIMENT + ONEBACK +
	SYLL_MAXF0 * SYLL_MEAN_INT * SYLL_DUR + (1|SUBJECT) + (1|ITEM),
	data = perception_data, family = binomial)
# Model failed to converge with max|grad| = 0.00746727 (tol = 0.001, component 1)

One_Intercept<- glmer(
	USER_RESP ~ PRIMARY_STRING + FUNCTION + EXPERIMENT + ONEBACK +
	SYLL_MAXF0 * SYLL_MEAN_INT * SYLL_DUR + (1|SUBJECT),
	data = perception_data, family = binomial)
# Converges


##### stan
di <- getStanDataInfo(formula = max_form, data = perception_data,
	subj = "SUBJECT", item = "ITEM")
Full_Model_stan <- stan(file = "glmer_logistic.stan", data = di$data,
	chains = 3, iter = 2000, warmup = 1000, refresh = 100,
	pars = di$info$keep, open_progress = TRUE,
	sample_file = "perception_study.stancsv")
Full_Model_stan_summ <- getSummary(Full_Model_stan,di)
printSummary(Full_Model_stan_summ)
# divergent transitions
# Time: 5.5 hours

Full_Model_stan <- stan(file = "glmer_logistic.stan", data = di$data,
	chains = 3, iter = 2000, warmup = 1000, refresh = 100,
	pars = di$info$keep, open_progress = TRUE, control = list(adapt_delta = 0.99),
	sample_file = "perception_study.stancsv")
Full_Model_stan_summ <- getSummary(Full_Model_stan,di)
printSummary(Full_Model_stan_summ)
# Converges
# Time: 10 hours
