############################################
# Test script

############################################

library(rstanarm)
    cpus_per_task <- Sys.getenv('SLURM_CPUS_PER_TASK')
    if(cpus_per_task!=''){
        options(mc.cores=as.integer(cpus_per_task))
    } else {
        options(mc.cores=parallel::detectCores() - 1)
    }

############################################

set.seed(88)

df <- data.frame(x=rnorm(100, 8, 1), y=rnorm(100, 2*x, 3))

mod_test <- stan_glm(y ~ x, data=df, family=gaussian(link='identity'), algorithm='optimizing')

summary(mod_test)
