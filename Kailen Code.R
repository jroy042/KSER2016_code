library(lme4)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(contrasts=c("contr.sum","contr.poly"))
source("getStanDataInfo.R")
source("getSummary.R")
source("printSummary.R")

behavioral_data <- read.csv("behavioral_task_data.csv",header=T)
behavioral_data$Syllables <- factor(behavioral_data$Syllables,ordered=T)
behavioral_data$Trial <- scale(behavioral_data$Trial)[,1]
behavioral_data$SUBTLEX_LogFrequency <- scale(behavioral_data$SUBTLEX_LogFrequency)[,1]

max_form <- Accuracy ~ (Trial + GoNoGo) * GoNoGo_Group * HandDecision +
	Gender + InitialSound + SUBTLEX_LogFrequency + Syllables +
	(1 + Trial * HandDecision + GoNoGo * HandDecision + Gender +
	InitialSound + SUBTLEX_LogFrequency + Syllables | Subject) +
	(1 + (Trial + GoNoGo) * GoNoGo_Group * HandDecision | Item)

	
##### lme4
Full_Model_lme4 <- glmer(formula = max_form, data = behavioral_data, family = binomial)
# Model failed to converge: degenerate  Hessian with 2 negative eigenvalues
# Time: 5.75 hours

Full_Model.bobyqa <- glmer(
	formula = max_form, data = behavioral_data, family = binomial,
	glmerControl(optimizer="bobyqa"))
# Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

Full_Model.Nelder_Mead <- glmer(
	formula = max_form, data = behavioral_data, family = binomial,
	glmerControl(optimizer="Nelder_Mead"))
# Model failed to converge: degenerate  Hessian with 15 negative eigenvalues

Intercepts_Only <- glmer(
	Accuracy ~ (Trial + GoNoGo) * GoNoGo_Group * HandDecision +
	Gender + InitialSound + SUBTLEX_LogFrequency + Syllables +
	(1|Subject) + (1|Item),
	family = binomial, data = behavioral_data)
# Model failed to converge with max|grad| = 0.00426312 (tol = 0.001, component 1)

One_Intercept<- glmer(
	Accuracy ~ (Trial + GoNoGo) * GoNoGo_Group * HandDecision +
	Gender + InitialSound + SUBTLEX_LogFrequency + Syllables + (1|Subject),
	family = binomial, data = behavioral_data)
# Converges


##### stan
di <- getStanDataInfo(formula = max_form, data = behavioral_data,
	subj = "Subject", item = "Item")
Full_Model_stan <- stan(file = "glmer_logistic.stan", data = di$data,
	chains = 3, iter = 2000, warmup = 1000, refresh = 100,
	pars = di$info$keep, open_progress = TRUE,
	sample_file = "behavioral_task.stancsv")
Full_Model_stan_summ <- getSummary(Full_Model_stan,di)
printSummary(Full_Model_stan_summ)
# divergent transitions
# Time: 1.5 hours

Full_Model_stan <- stan(file = "glmer_logistic.stan", data = di$data,
	chains = 3, iter = 2000, warmup = 1000, refresh = 100,
	pars = di$info$keep, open_progress = TRUE, control = list(adapt_delta = 0.99),
	sample_file = "behavioral_task.stancsv")
Full_Model_stan_summ <- getSummary(Full_Model_stan,di)
printSummary(Full_Model_stan_summ)
# Converges
# Time: 3 hours
