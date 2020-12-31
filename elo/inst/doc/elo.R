## -----------------------------------------------------------------------------
library(elo)

## -----------------------------------------------------------------------------
elo.A <- c(1500, 1500)
elo.B <- c(1500, 1600)
elo.prob(elo.A, elo.B)

## -----------------------------------------------------------------------------
wins.A <- c(1, 0)
elo.update(wins.A, elo.A, elo.B, k = 20)

## -----------------------------------------------------------------------------
elo.calc(wins.A, elo.A, elo.B, k = 20)

## -----------------------------------------------------------------------------
data(tournament)
str(tournament)

## -----------------------------------------------------------------------------
tournament$wins.A <- tournament$points.Home > tournament$points.Visitor
elo.run(wins.A ~ team.Home + team.Visitor, data = tournament, k = 20)
elo.run(score(points.Home, points.Visitor) ~ team.Home + team.Visitor, data = tournament, k = 20)

## -----------------------------------------------------------------------------
elo.run(score(points.Home, points.Visitor) ~ team.Home + team.Visitor +
        k(20*log(abs(points.Home - points.Visitor) + 1)), data = tournament)

## -----------------------------------------------------------------------------
k1 <- 20*log(abs(tournament$points.Home - tournament$points.Visitor) + 1)
elo.run(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + k(k1, k1/2), data = tournament)

## -----------------------------------------------------------------------------
elo.run(score(points.Home, points.Visitor) ~ adjust(team.Home, 10) + team.Visitor,
        data = tournament, k = 20)

## -----------------------------------------------------------------------------
tournament$elo.Visitor <- 1500
elo.run(score(points.Home, points.Visitor) ~ team.Home + elo.Visitor,
        data = tournament, k = 20)

## -----------------------------------------------------------------------------
tournament$elo.Visitor <- 1500
elo.run(score(points.Home, points.Visitor) ~ team.Home + elo.Visitor +
        regress(half, 1500, 0.2),
        data = tournament, k = 20)

## -----------------------------------------------------------------------------
e <- elo.run(score(points.Home, points.Visitor) ~ team.Home + team.Visitor,
             data = tournament, k = 20)
summary(e)
rank.teams(e)

## -----------------------------------------------------------------------------
head(as.matrix(e))

## -----------------------------------------------------------------------------
str(as.data.frame(e))

## -----------------------------------------------------------------------------
final.elos(e)

## -----------------------------------------------------------------------------
results <- elo.run(score(points.Home, points.Visitor) ~ adjust(team.Home, 10) + team.Visitor,
                   data = tournament, k = 20)
newdat <- data.frame(
  team.Home = "Athletic Armadillos",
  team.Visitor = "Blundering Baboons"
)
predict(results, newdata = newdat)

## -----------------------------------------------------------------------------
custom_update <- function(wins.A, elo.A, elo.B, k, adjust.A, adjust.B, ...)
{
  k*(wins.A - elo.prob(elo.A, elo.B, adjust.B = adjust.B,
                       adjust.A = ifelse(elo.A > 1500, adjust.A / 2, adjust.A)))
}
custom_prob <- function(elo.A, elo.B, adjust.A, adjust.B)
{
  1/(1 + 10^(((elo.B + adjust.B) - (elo.A + ifelse(elo.A > 1500, adjust.A / 2, adjust.A)))/400.0))
}
er2 <- elo.run2(score(points.Home, points.Visitor) ~ adjust(team.Home, 10) + team.Visitor,
                data = tournament, k = 20, prob.fun = custom_prob, update.fun = custom_update)
final.elos(er2)

## -----------------------------------------------------------------------------
er3 <- elo.run(score(points.Home, points.Visitor) ~ adjust(team.Home, 10) + team.Visitor,
               data = tournament, k = 20)
final.elos(er3)

## -----------------------------------------------------------------------------
dat <- data.frame(elo.A = c(1500, 1500), elo.B = c(1500, 1600),
                  wins.A = c(1, 0), k = 20)
form <- wins.A ~ elo.A + elo.B + k(k)
elo.prob(form, data = dat)
elo.update(form, data = dat)
elo.calc(form, data = dat)

## -----------------------------------------------------------------------------
elo.prob(~ elo.A + elo.B, data = dat)

## -----------------------------------------------------------------------------
elo.calc(wins.A ~ adjust(elo.A, 10) + elo.B + k(k), data = dat)

## -----------------------------------------------------------------------------
e.winpct <- elo.winpct(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + group(week), data = tournament,
                       subset = points.Home != points.Visitor) # to get rid of ties for now
summary(e.winpct)
rank.teams(e.winpct)
predict(e.winpct, newdata = data.frame(team.Home = "Athletic Armadillos", team.Visitor = "Blundering Baboons", stringsAsFactors = FALSE))

tournament$neutral <- replace(rep(0, nrow(tournament)), 30:35, 1)
summary(elo.winpct(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + neutral(neutral) + group(week),
                   data = tournament, subset = points.Home != points.Visitor))

## -----------------------------------------------------------------------------
e.winpct <- elo.winpct(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + group(week), data = tournament,
                       subset = points.Home != points.Visitor, running = TRUE, skip = 5)
summary(e.winpct)
predict(e.winpct, newdata = data.frame(team.Home = "Athletic Armadillos", team.Visitor = "Blundering Baboons", stringsAsFactors = FALSE)) # the same thing

## -----------------------------------------------------------------------------
results <- elo.glm(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + group(week), data = tournament,
                   subset = points.Home != points.Visitor) # to get rid of ties for now
summary(results)
rank.teams(results)
predict(results, newdata = data.frame(team.Home = "Athletic Armadillos", team.Visitor = "Blundering Baboons", stringsAsFactors = FALSE))

summary(elo.glm(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + neutral(neutral) + group(week),
                data = tournament, subset = points.Home != points.Visitor))

## -----------------------------------------------------------------------------
results <- elo.glm(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + group(week), data = tournament,
                   subset = points.Home != points.Visitor, running = TRUE, skip = 5)
summary(results)
predict(results, newdata = data.frame(team.Home = "Athletic Armadillos", team.Visitor = "Blundering Baboons", stringsAsFactors = FALSE)) # the same thing

## -----------------------------------------------------------------------------
mc <- elo.markovchain(score(points.Home, points.Visitor) ~ team.Home + team.Visitor, data = tournament,
                      subset = points.Home != points.Visitor, k = 0.7)
summary(mc)
rank.teams(mc)
predict(mc, newdata = data.frame(team.Home = "Athletic Armadillos", team.Visitor = "Blundering Baboons", stringsAsFactors = FALSE))
summary(elo.markovchain(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + neutral(neutral),
                        data = tournament, subset = points.Home != points.Visitor, k = 0.7))

## -----------------------------------------------------------------------------
mc <- elo.markovchain(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + group(week), data = tournament,
                      subset = points.Home != points.Visitor, k = 0.7, running = TRUE, skip = 5)
summary(mc)
predict(mc, newdata = data.frame(team.Home = "Athletic Armadillos", team.Visitor = "Blundering Baboons", stringsAsFactors = FALSE)) # the same thing

## ----eval=FALSE---------------------------------------------------------------
#  elo.markovchain(floor(wins.home) ~ team.home + team.visitor + k(ifelse(x > 0, rH(x), 1 - rH(x))))

## ----eval=FALSE---------------------------------------------------------------
#  elo.markovchain(ifelse(home.points - visitor.points > h, 1, 0) ~ team.home + team.visitor + k(pmax(rH(x), 1 - rH(x))))

## -----------------------------------------------------------------------------
co <- elo.colley(score(points.Home, points.Visitor) ~ team.Home + team.Visitor, data = tournament,
                 subset = points.Home != points.Visitor)
summary(co)
rank.teams(co)
predict(co, newdata = data.frame(team.Home = "Athletic Armadillos", team.Visitor = "Blundering Baboons", stringsAsFactors = FALSE))
summary(elo.colley(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + neutral(neutral),
                   data = tournament, subset = points.Home != points.Visitor))

## -----------------------------------------------------------------------------
co <- elo.colley(score(points.Home, points.Visitor) ~ team.Home + team.Visitor + group(week), data = tournament,
                      subset = points.Home != points.Visitor, running = TRUE, skip = 5)
summary(co)
predict(co, newdata = data.frame(team.Home = "Athletic Armadillos", team.Visitor = "Blundering Baboons", stringsAsFactors = FALSE)) # the same thing

## -----------------------------------------------------------------------------
summary(elo.glm(mov(points.Home, points.Visitor) ~ team.Home + team.Visitor, data = tournament,
                family = "gaussian"))

