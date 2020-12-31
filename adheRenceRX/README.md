
<!-- README.md is generated from README.Rmd. Please edit that file -->

# adheRenceRX

Check out our site [adheRenceRX](https://btbeal.github.io/adheRenceRX/)
<!-- badges: start --> <!-- badges: end -->

The goal of adheRenceRX is to provide a slightly opinionated set of
functions to allow researchers to assess medication adherence in the
most flexible way possible. The goal was (is) to write piping-friendly
verbs the “tidy” way to allow users to manipulate their data as they’d
like without storing data multiple times into their environment. In tidy
fashion, we aimed to create functions that did only one thing, ideally
that thing is obviated by the name of the function\! So, the package
makes assessing adherence as flexible as possible with some key things
left in the hands of the researcher. The final value is that functions
without vectorised solutions (`propagate_date()` and `rank_episodes()`)
are written with C++ allowing speed and performance when you’d rather do
research than run a function for an hour\!

This was a lot of fun to build but is still in production. If you find
errors, or know things you’d like to see done differently, reach out\!

## Installation

You can install the development version from
[GitHub](https://github.com/btbeal/adheRenceRX) with:

``` r
# install.packages("devtools")
devtools::install_github("btbeal/adheRenceRX")
```

## Overview

Much of the inspiration for this package came from conversations with
analysts who struggle to deal with the non-intuitive ways to deal with
medication adherence calculations from pharmaceutical claims data.

Our package is built around suggestions from Canfield and colleagues
(2019) who note that overlapping fill dates should be pushed forward and
never counted backwards, to assess adherence properly. For that reason,
our package revolves around the first step of creating adjusted dates
prior to any other calculation. Next, one can identify the gaps, rank
episodes of care, and calculate pdc. The purpose of the package was to
be as flexible as possible. So, there will be a lot left to be done by
the researcher (on purpose\!). For example, are there time periods
you’re particularly concerned with? Patient filters? Other groupings
(maybe episode of care?). Those are meant to be defined with dplyr verbs
outside of our functions.

Our verbs to date are:

  - `propagate_date()`
  - `identify_gaps()` or `summarise_gaps()`
  - `rank_episodes()`
  - `calculate_pdc()`

For the most part, our verbs assume that dates have been propagated
forward and gaps have been properly identified. This is on purpose but
is subject to change in the future.

## Examples

More examples of use can be found on within each functions
documentation; however, this should provide a decent overview of how the
package is to be used.

``` r
library(adheRenceRX)
library(dplyr)

# manipulate toy_claims, which has IDs based on the Canfield 2019 paper 
toy_claims %>% 
  # filter for some interesting IDs
  filter(ID %in% c("B", "D")) %>% 
  # Group by them (grouping not limited, of course)
  group_by(ID) %>% 
  # propagate the dates forward within those groups
  propagate_date(.date_var = date, .days_supply_var = days_supply)
#> # A tibble: 14 x 4
#> # Groups:   ID [2]
#>    ID    date       days_supply adjusted_date
#>    <chr> <date>           <dbl> <date>       
#>  1 B     2020-01-01          30 2020-01-01   
#>  2 B     2020-01-31          30 2020-01-31   
#>  3 B     2020-03-01          30 2020-03-01   
#>  4 B     2020-05-30          60 2020-05-30   
#>  5 B     2020-06-29          60 2020-07-29   
#>  6 B     2020-07-29          30 2020-09-27   
#>  7 B     2020-08-28          30 2020-10-27   
#>  8 B     2020-09-27          30 2020-11-26   
#>  9 D     2020-01-01          60 2020-01-01   
#> 10 D     2020-01-31          60 2020-03-01   
#> 11 D     2020-03-01          60 2020-04-30   
#> 12 D     2020-05-30          30 2020-06-29   
#> 13 D     2020-08-28          60 2020-08-28   
#> 14 D     2020-09-27          30 2020-10-27
```

Notice that several rows have been pushed forward to account for
overlaps in date. Also notice that the output changes the date and days
supply variable to `date` and `days_supply` while adding an
`adjusted_date` variable. The `adjusted_date` variable is used by some
of the other functions so it is important to complete this step first.

Once the dates have been adjusted, we can identify gaps in therapy with
`identify_gaps()` or summarise them with `summarise_gaps()`.

``` r
# The same code from above
toy_claims %>% 
  filter(ID %in% c("B", "D")) %>% 
  group_by(ID) %>% 
  propagate_date(.date_var = date, .days_supply_var = days_supply) %>% 
  # But now we can identify gaps
  identify_gaps()
#> # A tibble: 14 x 5
#> # Groups:   ID [2]
#>    ID    date       days_supply adjusted_date   gap
#>    <chr> <date>           <dbl> <date>        <dbl>
#>  1 B     2020-01-01          30 2020-01-01        0
#>  2 B     2020-01-31          30 2020-01-31        0
#>  3 B     2020-03-01          30 2020-03-01        0
#>  4 B     2020-05-30          60 2020-05-30       60
#>  5 B     2020-06-29          60 2020-07-29        0
#>  6 B     2020-07-29          30 2020-09-27        0
#>  7 B     2020-08-28          30 2020-10-27        0
#>  8 B     2020-09-27          30 2020-11-26        0
#>  9 D     2020-01-01          60 2020-01-01        0
#> 10 D     2020-01-31          60 2020-03-01        0
#> 11 D     2020-03-01          60 2020-04-30        0
#> 12 D     2020-05-30          30 2020-06-29        0
#> 13 D     2020-08-28          60 2020-08-28       30
#> 14 D     2020-09-27          30 2020-10-27        0


# Or, we could just summarise them all:
toy_claims %>% 
  filter(ID %in% c("B", "D")) %>% 
  group_by(ID) %>% 
  propagate_date(.date_var = date, .days_supply_var = days_supply) %>% 
  # Summarising gaps
  summarise_gaps()
#> # A tibble: 2 x 2
#>   ID    Summary_Of_Gaps
#>   <chr>           <dbl>
#> 1 B                  60
#> 2 D                  30
```

With the gaps identified, we can check for episodes of care using our
`rank_episodes()` functions. Note that this function assumes that you’ve
propagated your dates appropriately and identified all gaps. You can
then tell our function what can be considered a permissible gap, and
everything after a gap that large or more will be considered the next
episode\! Let me show you.

``` r
# The same code from above
toy_claims %>% 
  filter(ID %in% c("B", "D")) %>% 
  group_by(ID) %>% 
  propagate_date(.date_var = date, .days_supply_var = days_supply) %>% 
  # But now we can identify gaps
  identify_gaps() %>% 
  # say that anything over a 10 day gap should count as the next episode
  rank_episodes(.permissible_gap = 10)
#> # A tibble: 14 x 6
#> # Groups:   ID [2]
#>    ID    date       days_supply adjusted_date   gap episode
#>    <chr> <date>           <dbl> <date>        <dbl>   <dbl>
#>  1 B     2020-01-01          30 2020-01-01        0       1
#>  2 B     2020-01-31          30 2020-01-31        0       1
#>  3 B     2020-03-01          30 2020-03-01        0       1
#>  4 B     2020-05-30          60 2020-05-30       60       2
#>  5 B     2020-06-29          60 2020-07-29        0       2
#>  6 B     2020-07-29          30 2020-09-27        0       2
#>  7 B     2020-08-28          30 2020-10-27        0       2
#>  8 B     2020-09-27          30 2020-11-26        0       2
#>  9 D     2020-01-01          60 2020-01-01        0       1
#> 10 D     2020-01-31          60 2020-03-01        0       1
#> 11 D     2020-03-01          60 2020-04-30        0       1
#> 12 D     2020-05-30          30 2020-06-29        0       1
#> 13 D     2020-08-28          60 2020-08-28       30       2
#> 14 D     2020-09-27          30 2020-10-27        0       2
```

Finally, an actual adherence calculation. This is fairly straightforward
since the bulk of the work has been done adjusting your dates and then
appropriately identifying the gaps in therapy. Still, more functions =
more fun\!

``` r
toy_claims %>% 
  group_by(ID) %>% 
  propagate_date(.date_var = date, .days_supply_var = days_supply) %>% 
  identify_gaps() %>% 
  calculate_pdc()
#> # A tibble: 3 x 4
#>   ID    total_gaps total_days adherence
#>   <chr>      <dbl>      <dbl>     <dbl>
#> 1 A              0        270     1    
#> 2 B             60        330     0.818
#> 3 D             30        300     0.9
```

## Enjoy\!

That’s all we have for now. Again, this package is meant to provide some
helper functions with the meat of the project coming from our
`propagate_date()` and `rank_episodes()`. Notbaly, those tasks can’t be
accomplished with `dplyr` alone (as they do not have vectorised
solutions). For this reason, we’ve written some C++ functions to help
you speed up the task\!

## Citations

1.  Canfield SL, Zuckerman A, Anguiano RH, Jolly JA, DeClercq
    J.Navigating the wild west of medication adherence reporting in
    specialty pharmacy. J Manag Care Spec Pharm. 2019;25(10):1073-77.
