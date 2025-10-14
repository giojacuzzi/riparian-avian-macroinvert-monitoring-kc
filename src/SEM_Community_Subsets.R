## Simple Capstone SEM

library(piecewiseSEM)

data = read.csv('data/processed/site_ripoblig_validated.csv')

# Break down component regressions
model_bibi = lm(
  mean_BIBI ~ totalimp + Tree.Forest.High.Vegetation,
  data
)
model_alpha_total = glm(
  alpha ~ mean_BIBI + totalimp + stream_dist + Shrub.Low.Vegetation + Tree.Forest.High.Vegetation,
  data,
  family = poisson(link = "log")
)
model_alpha_riparian = glm(
  species_count_val ~ mean_BIBI + totalimp + stream_dist + Shrub.Low.Vegetation + Tree.Forest.High.Vegetation,
  data,
  family = poisson(link = "log")
)

# Use the `psem` function to create the SEMs
sem_alpha_total = psem(
  model_bibi,
  model_alpha_total
)
sem_alpha_riparian = psem(
  model_bibi,
  model_alpha_riparian
)

# Look at and plot objects, checking structure
sem_alpha_total
plot(sem_alpha_total)
sem_alpha_riparian
plot(sem_alpha_riparian)

# Conduct tests of directed separation (for each missing path)
# Establish the basis set & evaluate independence claims

# Use `dsep` function to perform the tests automagically
dSep(sem_alpha_total)
dSep(sem_alpha_riparian)

# Use `fisherC` function to evaluate claims
# A significant global Fisher’s C p-value (< 0.05) suggests that the modeled structure is statistically significantly different than the structure implied by the data, and that alternative pathways or causal links with missing variables warrant further exploration
fisherC(sem_alpha_total) # P > 0.05 => model fits well
fisherC(sem_alpha_riparian)

# Nagelkerke R2 describes proportion of variance explained by the model
summary(sem_alpha_total)
summary(sem_alpha_riparian)
