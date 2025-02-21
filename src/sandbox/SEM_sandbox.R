library(piecewiseSEM)

data(keeley)
head(keeley)

## Fit individual models
abiotic  = lm(abiotic ~ distance, keeley)
hetero   = lm(hetero ~ distance, keeley)
richness = lm(rich ~ abiotic + hetero, keeley)

# Wrap models into a SEM
model = psem(abiotic, hetero, richness, data = keeley)

# View structural equations and double-check it corroborates path diagram
model
plot(model)

## Evaluate fit

# Fit independence claims
dsep1 = lm(abiotic ~ hetero + distance, keeley)
dsep2 = lm(rich ~ distance + abiotic + hetero, keeley)

summary(dsep1)
summary(dsep2) # violation! there is a significant link between distance and richness

fisherC(model) # a low p-value is bad -- our model doesn't match the associations implied by our data

