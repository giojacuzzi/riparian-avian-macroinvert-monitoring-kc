## Seasonal SEM 

library(piecewiseSEM)

data = read.csv('data/processed/site_ripoblig_validated.csv')

# Break down component regressions
model_bibi = lm(
  mean_BIBI ~ totalimp + Tree.Forest.High.Vegetation,
  data
)
model_early = glm(
  early_count ~ mean_BIBI + totalimp + stream_dist + Shrub.Low.Vegetation + Tree.Forest.High.Vegetation,
  data,
  family = poisson(link = "log")
)
model_late = glm(
  late_count ~ mean_BIBI + totalimp + stream_dist + Shrub.Low.Vegetation + Tree.Forest.High.Vegetation,
  data,
  family = poisson(link = "log")
)

# Use the `psem` function to create the SEMs
sem_early = psem(
  model_bibi,
  model_early
)
sem_late = psem(
  model_bibi,
  model_late
)

# Look at and plot objects, checking structure
sem_early
plot(sem_early)

sem_late
plot(sem_late)

# Conduct tests of directed separation (for each missing path)
# Establish the basis set & evaluate independence claims

# Use `dsep` function to perform the tests automagically
dSep(sem_early)
dSep(sem_late)

# Use `fisherC` function to evaluate claims
fisherC(sem_early) # P > 0.05 == model fits well
fisherC(sem_late)


summary(sem_early)
summary(sem_late)
