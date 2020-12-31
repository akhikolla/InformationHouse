
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FlyingR

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/BMasinde/flight.svg?branch=master)](https://travis-ci.org/BMasinde/flight)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/BMasinde/flight?branch=master&svg=true)](https://ci.appveyor.com/project/BMasinde/flight)
[![Coveralls test
coverage](https://coveralls.io/repos/github/BMasinde/flight/badge.svg)](https://coveralls.io/github/BMasinde/flight?branch=master)
<!-- badges: end -->

The package provides methods for predicting flight range of birds based
on their physiological characteristics. This is an R implementation of
Flight program provided by Pennycuick.

## Installation

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BMasinde/FlyingR")
```

## Examples

### Time-Marching computation

``` r
library(FlyingR)
#> Welcome to package FlyingR
## basic example code

## birds comes with the package
data("birds")

simulation <- migrate(data = birds,  method = "cmm", settings = list(airDensity = 0.905))

simulation$range
#>           Anser anser  Hydrobates pelagicus   Pachyptila desolata 
#>              3537.296              3060.701              4274.620 
#>       Regulus regulus      Calidris canutus     Aegypius monachus 
#>              1264.523              4464.400              3979.271 
#>      Limosa lapponica           Anas crecca       Hirundo rustica 
#>             16652.710              3981.401              3445.013 
#>         Cygnus cygnus          Sylvia borin     Luscinia luscinia 
#>              3489.445              2663.630              2162.845 
#>       Corvus monedula         Anas penelope   Fregata magnificens 
#>              2300.933              6436.829             11904.426 
#>      Larus ridibundus      Diomedea exulans   Phalacrocorax carbo 
#>              6211.954              6271.194              2948.927 
#>       Gyps rueppellii   Torgos tracheliotus         Ardeotis kori 
#>              7867.072              6901.371              4631.958 
#>      Sturnus vulgaris     Fringilla coelebs      Carduelis spinus 
#>              4390.668              2821.588              2625.428 
#>     Turdus philomelos Calidris tenuirostris     Buteo swainsoni M 
#>              3184.748              7060.782              5979.198 
#>     Buteo swainsoni F 
#>              7868.848
```

The function also returns the mechanical and chemical power during the
simulation

### Range estimation based on ODE

This function estimates the range based on Pennycuick (1975) Mechanics
of Flight where Breguet set of equations are used.

``` r
## when estimating range of a single bird
birds_range <- flysim(data = birds,  settings = list(airDensity = 0.905))

birds_range$range
#>           Anser anser  Hydrobates pelagicus   Pachyptila desolata 
#>                3193.8                3252.9                4192.0 
#>       Regulus regulus      Calidris canutus     Aegypius monachus 
#>                1521.6                4240.6                3951.6 
#>      Limosa lapponica           Anas crecca       Hirundo rustica 
#>               11209.3                3797.8                3898.0 
#>         Cygnus cygnus          Sylvia borin     Luscinia luscinia 
#>                3296.0                2801.7                2301.9 
#>       Corvus monedula         Anas penelope   Fregata magnificens 
#>                2452.7                5428.0               10527.4 
#>      Larus ridibundus      Diomedea exulans   Phalacrocorax carbo 
#>                6113.1                5436.8                2918.8 
#>       Gyps rueppellii   Torgos tracheliotus         Ardeotis kori 
#>                6808.3                6334.7                4211.5 
#>      Sturnus vulgaris     Fringilla coelebs      Carduelis spinus 
#>                4197.4                3036.6                2926.6 
#>     Turdus philomelos Calidris tenuirostris     Buteo swainsoni M 
#>                3243.8                6053.5                5671.6 
#>     Buteo swainsoni F 
#>                6994.8
```

## The data

*birds* definitions pulled from Flight program in-built datasets and fat
mass randomly generated where initially zero. In addition, by default
muscle mass was derived as 0.17 fraction of the all-up mass. User’s data
should have columns named appropriately. The package looks for columns
named *id, name or species.name*, *bodymass or allupmass*, *wingspan,
ws*, *wingarea*, *ordo, order* (which is a factored column with levels 1
or 2 passerines and non-passerines respectively) *fatmass, fat.mass,
fat\_mass* and lastly *muscle\_mass*.

``` r
birds
#>          Scientific.name Empty.mass Wing.span Fat.mass Order Wing.area
#> 1            Anser anser    3.77000     1.600  0.84641     2   0.33100
#> 2   Hydrobates pelagicus    0.02580     0.355  0.00591     2   0.01610
#> 3    Pachyptila desolata    0.15500     0.637  0.03886     2   0.04710
#> 4        Regulus regulus    0.00542     0.156  0.00112     1   0.00525
#> 5       Calidris canutus    0.12700     0.538  0.03500     2   0.03320
#> 6      Aegypius monachus    9.90000     3.040  2.02565     2   1.40000
#> 7       Limosa lapponica    0.36700     0.748  0.20112     2   0.05680
#> 8            Anas crecca    0.23500     0.582  0.06562     2   0.04580
#> 9        Hirundo rustica    0.01900     0.318  0.00570     1   0.01320
#> 10         Cygnus cygnus   12.50000     2.560  2.50000     2   0.75600
#> 11          Sylvia borin    0.02200     0.240  0.00660     1   0.01100
#> 12     Luscinia luscinia    0.02700     0.263  0.00675     1   0.01300
#> 13       Corvus monedula    0.18100     0.600  0.03620     1   0.06180
#> 14         Anas penelope    0.77000     0.822  0.28607     2   0.08290
#> 15   Fregata magnificens    1.67000     2.140  0.55799     2   0.37200
#> 16      Larus ridibundus    0.28500     0.967  0.07881     2   0.09920
#> 17      Diomedea exulans    9.57000     3.060  2.12836     2   0.64400
#> 18   Phalacrocorax carbo    2.56000     1.350  0.50768     2   0.22400
#> 19       Gyps rueppellii    7.30000     2.500  2.52588     2   0.89200
#> 20   Torgos tracheliotus    6.60000     2.640  2.01454     2   1.03000
#> 21         Ardeotis kori   11.90000     2.470  3.45889     2   1.06000
#> 22      Sturnus vulgaris    0.08190     0.384  0.02973     1   0.02530
#> 23     Fringilla coelebs    0.02300     0.264  0.00690     1   0.01310
#> 24      Carduelis spinus    0.01120     0.212  0.00336     1   0.00785
#> 25     Turdus philomelos    0.07160     0.361  0.02148     1   0.02250
#> 26 Calidris tenuirostris    0.23300     0.587  0.08970     2   0.03960
#> 27     Buteo swainsoni M    0.77500     1.250  0.22248     2   0.21000
#> 28     Buteo swainsoni F    1.06000     1.330  0.37117     2   0.24000
#>    Muscle.mass
#> 1    0.6409000
#> 2    0.0043860
#> 3    0.0263500
#> 4    0.0009214
#> 5    0.0215900
#> 6    1.6829999
#> 7    0.0623900
#> 8    0.0399500
#> 9    0.0032300
#> 10   2.1250000
#> 11   0.0037400
#> 12   0.0045900
#> 13   0.0307700
#> 14   0.1309000
#> 15   0.2839000
#> 16   0.0484500
#> 17   1.6268999
#> 18   0.4352000
#> 19   1.2410000
#> 20   1.1220000
#> 21   2.0229999
#> 22   0.0139230
#> 23   0.0039100
#> 24   0.0019040
#> 25   0.0121720
#> 26   0.0396100
#> 27   0.1317500
#> 28   0.1802000
```
