## -----------------------------------------------------------------------------
data("birds", package = "FlyingR")
results <-
  FlyingR::migrate(data = birds,
                  method = "cmm",
                  speed_control = 1,
                  min_energy_protein = 0.05)

# extract range as a vector
results$range

## -----------------------------------------------------------------------------
results <-
  FlyingR::migrate(data = birds,
                  method = "cmm",
                  settings = list(ipf = 0.9),
                  speed_control = 1,
                  min_energy_protein = 0.05)

## -----------------------------------------------------------------------------
results <-
  FlyingR::migrate(data = birds,
                  method = "csw",
                  settings = list(ipf = 0.9),
                  speed_control = 1,
                  min_energy_protein = 0.05)


# obtain remaining body mass
results$bodyMass

# starting minimum power speed
results$startMinSpeed

# end of flight minimum power speed
results$endMinSpeed

