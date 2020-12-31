############################# CONTRIVED DATA ###################################
# LAST UPDATE: NA

#' Contrived Example Data
#' 
#' These data are a simulated point-referenced geospatial data that serve to provide 
#' a clean example of a kriging model. There are 500 observations with coordinates 
#' located on a unit square.
#' 
#' @name ContrivedData
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @format The \code{ContrivedData} dataset has 500 observations and 5 variables.
#' \describe{
#'   \item{\code{y}}{The outcome variable. Its true population functional form is 
#'     \eqn{y_s=0+1 x_{1s}+2 x_{2s}+\omega_{s}+\epsilon_{s}}. The true variance of 
#'     \eqn{\omega} is \eqn{\sigma^2=0.5} and of \eqn{\epsilon} is \eqn{\tau^2=0.5}. 
#'     The decay term that shapes spatial correlation levels is \eqn{\phi=2.5}.}
#'     \item{\code{x.1}}{A predictor with a standard uniform distribution.}
#'     \item{\code{x.2}}{A predictor with a standard normal distribution.}
#'     \item{\code{s.1}}{Coordinate in eastings for each observation, distributed 
#'       standard uniform.}
#'     \item{\code{s.2}}{Coordinate in northings for each observation, distributed 
#'       standard uniform.}
#' }
#' @examples
#' \dontrun{
#' # Summarize example data
#' summary(ContrivedData)
#' 
#' # Initial OLS model
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData)
#' # summary(contrived.ols)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' #For simple illustration, we set to few iterations.
#' #In this case, a 10,000-iteration run converges to the true parameters.
#' #If you have considerable time and hardware, delete the # on the next line.
#' #10,000 iterations took 39 min. with 8 GB RAM & a 1.5 GHz Quad-Core processor.
#' M <- 100
#' #M<-10000
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'    data = ContrivedData, n.iter = M, n.burnin=20, range.tol = 0.05)
#' # Alternatively, use burnin() after estimation  
#' #contrived.run <- burnin(contrived.run, n.burnin=20)
#' 
#' # Summarize the results and examine results against true coefficients	
#' summary(contrived.run)
#' (TRUTH<-c(0.5,2.5,0.5,0,1,2))
#' }
NULL

################################ NEW YORK STATE ################################
# LAST UPDATE: NA

#' New York State CCES Respondents in 2008
#' 
#' These data are a subset of the 2008 Cooperative Congressional Election Survey 
#' (CCES) Common Content. Only 1108 respondents from the state of New York are included, 
#' with predictors drawn from Gill's (2020) model of self-reported ideology. The 
#' CCES data are merged with predictors on geographic location based on ZIP codes 
#' (from ArcGIS & TomTom) and county ruralism (from the USDA).
#' 
#' @name NY_subset
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @format The \code{NY_subset} dataset has 1108 observations and 26 variables.
#' \describe{
#'   \item{\code{state}}{The state abbreviation of the respondent's residence.}
#'   \item{\code{zip}}{The respondent's ZIP code.}
#'   \item{\code{age}}{The age of the respondent in years.}
#'   \item{\code{female}}{An indicator of whether the respondent is female.}
#'   \item{\code{ideology}}{The respondent's self-reported ideology on a scale of 0 (liberal) to 100 (conservative).}
#'   \item{\code{educ}}{The respondent's level of education. 0=No Highschool, 
#'     1=High School Graduate, 2=Some College, 3=2-year Degree, 4=4-year degree, 5=Post-Graduate.}
#'   \item{\code{race}}{The respondent's race. 1=White, 2=African American, 3=Nonwhite & nonblack.}
#'   \item{\code{empstat}}{The respondent's employment status. 1=employed, 2=unemployed, 3=not in workforce.}
#'   \item{\code{ownership}}{Indicator for whether the respondent owns his or her own home.}
#'   \item{\code{inc14}}{The respondent's self reported income. 1=Less than $10,000, 
#'     2=$10,000-$14,999, 3=$15,000-$19,000, 4=$20,000-$24,999, 5=$25,000-$29,999, 
#'     6=$30,000-$39,999, 7=$40,000-$49,999, 8=$50,000-$59,999, 9=$60,000-$69,999, 
#'     10=$70,000-$79,999, 11=$80,000-$89,999, 12=$100,000-$119,999, 13=$120,000-$149,999, 
#'     14=$150,000 or more.}
#'   \item{\code{catholic}}{Indicator for whether the respondent is Catholic.}
#'   \item{\code{mormon}}{Indicator for whether the respondent is Mormon.}
#'   \item{\code{orthodox}}{Indicator for whether the respondent is Orthodox Christian.}
#'   \item{\code{jewish}}{Indicator for whether the respondent is Jewish.}
#'   \item{\code{islam}}{Indicator for whether the respondent is Muslim.}
#'   \item{\code{mainline}}{Indicator for whether the respondent is Mainline Christian.}
#'   \item{\code{evangelical}}{Indicator for whether the respondent is Evangelical Christian.}
#'   \item{\code{FIPS_Code}}{FIPS code of the repondent's state.}
#'   \item{\code{rural}}{Nine-point USDA scale of the ruralism of each county, with 
#'     0 meaning the most urban and 8 meaning the most rural.}
#'   \item{\code{zipPop}}{Indicates the population of the repondent's ZIP code.}
#'   \item{\code{zipLandKM}}{Indicates the land area in square kilometers of the repondent's ZIP code.}
#'   \item{\code{weight}}{Survey weights created by the CCES.}
#'   \item{\code{cd}}{The congressional district the respondent resides in.}
#'   \item{\code{fipsCD}}{Index that fuses the state FIPS code in the first two 
#'     digits and the congressional district number in the last two digits.}
#'   \item{\code{northings}}{Indicates the geographical location of the respondent in kilometer-based northings.}
#'   \item{\code{eastings}}{Indicates the geographical location of the respondent in kilometer-based eastings.}
#' }
#' 
#' @source {
#' Ansolabehere, Stephen. 2011. "CCES, Common Content, 2008." Ver. 4.
#' 
#' ArcGIS. 2012. "USA ZIP Code Areas."  \url{https://www.arcgis.com/home/item.html?id=8d2012a2016e484dafaac0451f9aea24}
#' 
#' United States Department of Agriculture. 2013. "2013 Rural-Urban Continuum Codes." \url{https://www.ers.usda.gov/data-products/rural-urban-continuum-codes.aspx}
#' }
#' 
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'     \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#' 
#' @examples 
#' \dontrun{
#' ny <- NY_subset
#' 
#' #data cleaning
#' ny$cathOrth<-ny$catholic+ny$orthodox
#' ny$consRelig<-ny$mormon+ny$evangelical
#' ny$jewMus<-ny$jewish+ny$islam
#' 
#' # Explanatory Variable Matrix
#' psrm.data <-cbind(ny$age, ny$educ, I(ny$age*ny$educ), as.numeric(ny$race==2), 
#'        as.numeric(ny$race==3), ny$female, I(as.numeric(ny$race==2)*ny$female), 
#'        I(as.numeric(ny$race==3)*ny$female), ny$cathOrth, ny$consRelig, 
#'        ny$jewMus, ny$mainline, ny$rural, ny$ownership, 
#'        as.numeric(ny$empstat==2), as.numeric(ny$empstat==3),ny$inc14)
#'
#' dimnames(psrm.data)[[2]] <- c("Age", "Education", "Age.education", 
#'                              "African.American", "Nonwhite.nonblack","Female", 
#'                              "African.American.female", "Nonwhite.nonblack.female", 
#'                              "Catholic.Orthodox", "Evang.Mormon", "Jewish.Muslim", 
#'                              "Mainline","Ruralism", "Homeowner", "Unemployed",
#'                              "Not.in.workforce","Income")
#' 
#' # Outcome Variable
#' ideo <- matrix(ny$ideology,ncol=1)
#' 
#' # Set Number of Iterations:
#' # WARNING: 20 iterations is intensive on many machines.
#' # This example was tuned on Amazon Web Services (EC2) over many hours
#' # with 20,000 iterations--unsuitable in 2020 for most desktop machines.
#' #M<-20000 
#' M<-100 
#' set.seed(1,kind="Mersenne-Twister")
#' 
#' # Estimate the Model
#' ny.fit <- metropolis.krige(formula = ideo ~ psrm.data, coords = cbind(ny$eastings, ny$northings),
#'           powered.exp=1, n.iter=M, spatial.share=0.31,range.share=0.23,beta.var=10,
#'           range.tol=0.01, b.tune=0.1, nugget.tune=20, psill.tune=5)		
#'       
#' # Discard first 20% of Iterations as Burn-In (User Discretion Advised).
#' ny.fit <- burnin(ny.fit, M/5)
#' 
#' # Summarize Results
#' summary(ny.fit)
#' 
#' #Convergence Diagnostics: Geweke and Heidelberger-Welch
#' geweke(ny.fit)
#' heidel.welch(ny.fit)
#' 
#' # Draw Semivariogram
#' semivariogram(ny.fit)
#' }     
NULL

################################ NEW YORK CITY ################################
# LAST UPDATE: NA

#' New York City CCES Respondents in 2008
#' 
#' These data are a subset of the 2008 Cooperative Congressional Election Survey 
#' (CCES) Common Content. Only 568 respondents from New York City are included, 
#' with predictors drawn from Gill's (2020) model of self-reported ideology. The 
#' CCES data are merged with predictors on geographic location based on ZIP codes 
#' (from ArcGIS & TomTom) and county ruralism (from the USDA).
#' 
#' @name NYcity_subset
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @format The \code{NYcity_subset} dataset has 568 observations and 26 variables.
#' \describe{
#'   \item{\code{state}}{The state abbreviation of the respondent's residence.}
#'   \item{\code{zip}}{The respondent's ZIP code.}
#'   \item{\code{age}}{The age of the respondent in years.}
#'   \item{\code{female}}{An indicator of whether the respondent is female.}
#'   \item{\code{ideology}}{The respondent's self-reported ideology on a scale of 0 (liberal) to 100 (conservative).}
#'   \item{\code{educ}}{The respondent's level of education. 0=No Highschool, 
#'     1=High School Graduate, 2=Some College, 3=2-year Degree, 4=4-year degree, 5=Post-Graduate.}
#'   \item{\code{race}}{The respondent's race. 1=White, 2=African American, 3=Nonwhite & nonblack.}
#'   \item{\code{empstat}}{The respondent's employment status. 1=employed, 2=unemployed, 3=not in workforce.}
#'   \item{\code{ownership}}{Indicator for whether the respondent owns his or her own home.}
#'   \item{\code{inc14}}{The respondent's self reported income. 1=Less than $10,000, 
#'     2=$10,000-$14,999, 3=$15,000-$19,000, 4=$20,000-$24,999, 5=$25,000-$29,999, 
#'     6=$30,000-$39,999, 7=$40,000-$49,999, 8=$50,000-$59,999, 9=$60,000-$69,999, 
#'     10=$70,000-$79,999, 11=$80,000-$89,999, 12=$100,000-$119,999, 13=$120,000-$149,999, 
#'     14=$150,000 or more.}
#'   \item{\code{catholic}}{Indicator for whether the respondent is Catholic.}
#'   \item{\code{mormon}}{Indicator for whether the respondent is Mormon.}
#'   \item{\code{orthodox}}{Indicator for whether the respondent is Orthodox Christian.}
#'   \item{\code{jewish}}{Indicator for whether the respondent is Jewish.}
#'   \item{\code{islam}}{Indicator for whether the respondent is Muslim.}
#'   \item{\code{mainline}}{Indicator for whether the respondent is Mainline Christian.}
#'   \item{\code{evangelical}}{Indicator for whether the respondent is Evangelical Christian.}
#'   \item{\code{FIPS_Code}}{FIPS code of the repondent's state.}
#'   \item{\code{rural}}{Nine-point USDA scale of the ruralism of each county, with 
#'     0 meaning the most urban and 8 meaning the most rural.}
#'   \item{\code{zipPop}}{Indicates the population of the repondent's ZIP code.}
#'   \item{\code{zipLandKM}}{Indicates the land area in square kilometers of the repondent's ZIP code.}
#'   \item{\code{weight}}{Survey weights created by the CCES.}
#'   \item{\code{cd}}{The congressional district the respondent resides in.}
#'   \item{\code{fipsCD}}{Index that fuses the state FIPS code in the first two 
#'     digits and the congressional district number in the last two digits.}
#'   \item{\code{northings}}{Indicates the geographical location of the respondent in kilometer-based northings.}
#'   \item{\code{eastings}}{Indicates the geographical location of the respondent in kilometer-based eastings.}
#' }
#' 
#' @source {
#' Ansolabehere, Stephen. 2011. "CCES, Common Content, 2008." Ver. 4.
#' 
#' ArcGIS. 2012. "USA ZIP Code Areas."  \url{https://www.arcgis.com/home/item.html?id=8d2012a2016e484dafaac0451f9aea24}
#' 
#' United States Department of Agriculture. 2013. "2013 Rural-Urban Continuum Codes." \url{https://www.ers.usda.gov/data-products/rural-urban-continuum-codes.aspx}
#' }
#' 
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'     \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#'  
#' @examples         
#' \dontrun{
#' nyc <- NYcity_subset
#' 
#' #data cleaning
#' nyc$cathOrth<-nyc$catholic+nyc$orthodox
#' nyc$consRelig<-nyc$mormon+nyc$evangelical
#' nyc$jewMus<-nyc$jewish+nyc$islam
#' 
#' # Explanatory Variable Matrix
#' psrm.data <-cbind(nyc$age, nyc$educ, I(nyc$age*nyc$educ), as.numeric(nyc$race==2), 
#'        as.numeric(nyc$race==3), nyc$female, I(as.numeric(nyc$race==2)*nyc$female), 
#'        I(as.numeric(nyc$race==3)*nyc$female), nyc$cathOrth, nyc$consRelig, 
#'        nyc$jewMus, nyc$mainline, nyc$rural, nyc$ownership, 
#'        as.numeric(nyc$empstat==2), as.numeric(nyc$empstat==3),nyc$inc14)
#'
#' dimnames(psrm.data)[[2]] <- c("Age", "Education", "Age.education", 
#'                              "African.American", "Nonwhite.nonblack","Female", 
#'                              "African.American.female", "Nonwhite.nonblack.female", 
#'                              "Catholic.Orthodox", "Evang.Mormon", "Jewish.Muslim", 
#'                              "Mainline","Ruralism", "Homeowner", "Unemployed",
#'                              "Not.in.workforce","Income")
#' 
#' # Outcome Variable
#' ideo <- matrix(nyc$ideology,ncol=1)
#' 
#' # WARNING: This example was tuned on Amazon Web Services (EC2) over many hours
#' # with 150,000 iterations--a strain in 2020 for most desktop machines.
#' # A test with few iterations allows illustration.
#' #M<-150000 
#' M<-150 
#' set.seed(1,kind="Mersenne-Twister")
#' 
#' # Estimate the Model
#' nyc.fit <- metropolis.krige(formula = ideo ~ psrm.data, coords = cbind(nyc$eastings, nyc$northings),
#'           powered.exp=1, n.iter=M, spatial.share=0.31,range.share=0.23,beta.var=10,
#'           range.tol=0.01, b.tune=0.1, nugget.tune=20, psill.tune=5)		
#'       
#' # Discard first 20% of Iterations as Burn-In (User Discretion Advised).
#' nyc.fit <- burnin(nyc.fit, M/5)
#' 
#' # Summarize Results
#' summary(nyc.fit)
#' 
#' #Convergence Diagnostics: Geweke and Heidelberger-Welch
#' geweke(nyc.fit)
#' heidel.welch(nyc.fit)
#' 
#' # Draw Semivariogram
#' semivariogram(nyc.fit)
#' }     
NULL


################################ STATE COMBINED ################################
# LAST UPDATE: NA

#' Congressional District Public Opinion Ideology in 2010
#' 
#' These data present measures of ideology in 2010 for 434 districts for the U.S. 
#'   House of Representatives, recorded as the variable \code{krige.cong}. Forecasts 
#'   are based on a kriging model fitted over the 2008 Cooperative Congressional 
#'   Election Survey (CCES), paired with predictive data from the 2010 Census. Each 
#'   district's public ideology is paired with the DW-NOMINATE common space score 
#'   of each of its representative in 2011 (update from McCarty, Poole and Rosenthal 
#'   1997). Eight districts have repeated observations in order to include the DW-NOMINATE 
#'   score when a member was replaced mid-term.
#' 
#' @name congCombined
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @format The \code{congCombined} dataset has 442 observations and 12 variables. 4
#'   34 out of 435 congressional districts are covered, with eight districts duplicated 
#'   when a member was replaced mid-term.
#'   
#' \describe{
#'   \item{\code{stateCD}}{Unique identifier for each congressional district by state. 
#'     The first two digits are \code{STATEA}, and the second two are \code{cd}.}
#'  \item{\code{krige.cong}}{The ideology of the average citizen in the congressional district.}
#'  \item{\code{krige.state.var}}{The variance of ideology among the district's citizens.}
#'  \item{\code{cong}}{The term of Congress studied--112 for this dataset.}
#'  \item{\code{idno}}{Identification number for the House member--ICPSR numbers 
#'    continued by Poole & Rosenthal.}
#'  \item{\code{state}}{The ICPSR code for the state.}
#'  \item{\code{cd}}{The congressional district number.}
#'  \item{\code{statenm}}{The first seven letters of the state's name.}
#'  \item{\code{party}}{Political party of the House member. 100=Democrat, 200=Republican.}
#'  \item{\code{name}}{Last name of the House member, followed by first name if ambiguous.}
#'  \item{\code{dwnom1}}{First dimension DW-NOMINATE common space score for the House member. Higher values are usually interpreted as more right-wing, with lower values as more left-wing.}
#'  \item{\code{STATEA}}{The FIPS code for the state.}
#' }
#' 
#' @source {
#' Ansolabehere, Stephen. 2011. "CCES, Common Content, 2008." Ver. 4.
#' 
#' McCarty, Nolan M., Keith T. Poole and Howard Rosenthal. 1997. \emph{Income 
#'   Redistribution and the Realignment of American Politics}. American Enterprise 
#'   Institude Studies on Understanding Economic Inequality. Washington: AEI Press.
#' 
#' Minnesota Population Center. 2011. \emph{National Historical Geographic Information 
#'   System: Version 2.0.} Minneapolis, MN: University of Minnesota. \samp{https://www.nhgis.org}
#' }
#' 
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'     \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#'     
#' @examples     
#' # Descriptive Statistics
#' summary(congCombined)
#' 
#' # Correlate House Members' DW-NOMINATE Scores with Public Opinion Ideology
#' cor(congCombined$dwnom1,congCombined$krige.cong)
#' 
#' # Plot House Members' DW-NOMINATE Scores against Public Opinion Ideology
#' plot(y=congCombined$dwnom1,x=congCombined$krige.cong,
#'    xlab="District Ideology (Kriging)", ylab="Legislator Ideology (1st Dim., Common Space)", 
#'    main="U.S. House of Representatives", type="n")
#' points(y=congCombined$dwnom1[congCombined$party==200],
#'    x=congCombined$krige.cong[congCombined$party==200],pch="R",col="red")
#'    points(y=congCombined$dwnom1[congCombined$party==100],
#'    x=congCombined$krige.cong[congCombined$party==100],pch="D",col="blue")
NULL

################################ STATE COMBINED ################################
# LAST UPDATE: NA

#' State Public Opinion Ideology in 2010
#' 
#' These data present measures of ideology in 2010 for the 50 American states, 
#'   recorded as the variable \code{krige.state}. Forecasts are based on a kriging 
#'   model fitted over the 2008 Cooperative Congressional Election Survey (CCES), 
#'   paired with predictive data from the 2010 Census. Each state is listed twice, 
#'   as each state's public ideology is paired with the DW-NOMINATE common space 
#'   score of each of its two senators in 2011 (update from McCarty, Poole and 
#'   Rosenthal 1997).
#' 
#' @name stateCombined
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @format The \code{stateCombined} dataset has 100 observations (2 each for 50 states) and 13 variables.
#' \describe{
#'   \item{\code{STATEA}}{The FIPS code for the state.}
#'   \item{\code{krige.state}}{The ideology of the average citizen in the state.}
#'   \item{\code{krige.state.var}}{The variance of ideology among the state's citizens.}
#'   \item{\code{cong}}{The term of Congress studied--112 for this dataset.}
#'   \item{\code{idno}}{Identification number for the senator--ICPSR numbers continued by Poole & Rosenthal.}
#'   \item{\code{state}}{The ICPSR code for the state.}
#'   \item{\code{cd}}{The congressional district number--0 for senators.}
#'   \item{\code{statenm}}{The first seven letters of the state's name.}
#'   \item{\code{party}}{Political party of the senator. 100=Democrat, 200=Republican, 328=Independent.}
#'   \item{\code{name}}{Last name of the senator, followed by first name if ambiguous.}
#'   \item{\code{dwnom1}}{First dimension DW-NOMINATE common space score for the senator. 
#'     Higher values are usually interpreted as more right-wing, with lower values as more left-wing.}
#'   \item{\code{stateCD}}{Combined index of \code{STATEA} followed by \code{cd}.}
#'   \item{\code{obama}}{Barack Obama's percentage of the two-party vote in the state in 2012.}
#' }
#' 
#' @source {
#' Ansolabehere, Stephen. 2011. "CCES, Common Content, 2008." Ver. 4.
#' 
#' McCarty, Nolan M., Keith T. Poole and Howard Rosenthal. 1997. \emph{Income 
#'   Redistribution and the Realignment of American Politics}. American Enterprise 
#'   Institude Studies on Understanding Economic Inequality. Washington: AEI Press.
#' 
#' Minnesota Population Center. 2011. \emph{National Historical Geographic Information 
#'   System: Version 2.0.} Minneapolis, MN: University of Minnesota. \samp{https://www.nhgis.org}
#' }
#' 
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'     \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#'
#' @examples 
#' # Descriptive Statistics
#' summary(stateCombined)
#' 
#' # Correlate Senators' DW-NOMINATE Scores with Public Opinion Ideology
#' cor(stateCombined$krige.state,stateCombined$dwnom1)
#' 
#' # Plot Senators' DW-NOMINATE Scores against Public Opinion Ideology
#' plot(y=stateCombined$dwnom1,x=stateCombined$krige.state,
#'      xlab="State Ideology (Kriging)", ylab="Legislator Ideology (1st Dim., Common Space)", 
#'      main="U.S. Senate", type="n")
#' points(y=stateCombined$dwnom1[stateCombined$party==200],
#'        x=stateCombined$krige.state[stateCombined$party==200],pch="R",col="red")
#' points(y=stateCombined$dwnom1[stateCombined$party==100],
#'        x=stateCombined$krige.state[stateCombined$party==100],pch="D",col="blue")     
NULL

################################# STATE LOWER ##################################
# LAST UPDATE: NA

#' State Legislative District (Lower Chambers) Public Opinion Ideology in 2010
#' 
#' These data present measures of ideology in 2010 for the districts for lower 
#' chambers of state legislatures, recorded as the variable \code{krige.lower}. 
#' 49 states' chambers are covered--the Nebraska Unicameral is omitted here to be 
#' included in the file \code{upperCombined}. Forecasts are based on a kriging model 
#' fitted over the 2008 Cooperative Congressional Election Survey (CCES), paired 
#' with predictive data from the 2010 Census. Each district's public ideology is 
#' paired with a measure of the ideology of the State House member (or members) 
#' from the district (update from Shor and McCarty 2011).
#' 
#' @name lowerCombined
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @format The \code{lowerCombined} dataset has 5446 observations and 10 variables.
#' \describe{
#'   \item{\code{st}}{Two-letter postal abbreviation for the state.}
#'   \item{\code{lower}}{The state legislative district number (lower chamber).}     		
#'   \item{\code{STATEA}}{The FIPS code for the state.}
#'   \item{\code{krige.lower}}{The ideology of the average citizen in the district.}
#'   \item{\code{lowerKluge}}{Combined index of \code{STATEA} followed by \code{lower}.}
#'   \item{\code{krige.lower.var}}{The variance of ideology among the district's citizens.}
#'   \item{\code{name}}{Last name of the state legislator, followed by first name and middle initial.}
#'   \item{\code{party}}{Political party of the legislator. D=Democrat, R=Republican, X=Other.}
#'   \item{\code{st_id}}{Temporary identifer variable. DO NOT USE.}
#'   \item{\code{np_score}}{Ideology score for the state legislator (lower chamber). 
#'     Higher values are usually interpreted as more right-wing, with lower values as more left-wing.}
#' }
#' 
#' @source {
#' Ansolabehere, Stephen. 2011. "CCES, Common Content, 2008." Ver. 4.
#' 
#' Minnesota Population Center. 2011. \emph{National Historical Geographic Information 
#'   System: Version 2.0.} Minneapolis, MN: University of Minnesota. \samp{https://www.nhgis.org}
#'   
#'  Shor, Boris and Nolan M. McCarty. 2011. "The Ideological Mapping of American Legislatures." 
#'  \emph{American Political Science Review} 105(3):530-551.
#' }
#' 
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'     \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#' 
#' @examples 
#' # Descriptive Statistics
#' summary(lowerCombined)
#' 
#' # Correlate Senators' DW-NOMINATE Scores with Public Opinion Ideology
#' cor(lowerCombined$np_score,lowerCombined$krige.lower,use="complete.obs")
#' 
#' # Plot Legislators' DW-NOMINATE Scores against Public Opinion Ideology
#' plot(y=lowerCombined$np_score,x=lowerCombined$krige.lower,
#'      xlab="District Ideology (Kriging)", ylab="Legislator Ideology (Shor & McCarty)", 
#'      main="State Legislatures: Lower Chambers", type="n")#
#' points(y=lowerCombined$np_score[lowerCombined$party=="R"],
#'        x=lowerCombined$krige.lower[lowerCombined$party=="R"],pch=".",col="red")
#' points(y=lowerCombined$np_score[lowerCombined$party=="D"],
#'        x=lowerCombined$krige.lower[lowerCombined$party=="D"],pch=".",col="blue")    
NULL

################################# STATE UPPER ##################################
# LAST UPDATE: NA

#' State Legislative District (Upper Chambers) Public Opinion Ideology in 2010
#' 
#' These data present measures of ideology in 2010 for the districts for upper 
#'   chambers of state legislatures, recorded as the variable \code{krige.upper}. 
#'   All 50 states' chambers are covered (including the Nebraska Unicameral). 
#'   Forecasts are based on a kriging model fitted over the 2008 Cooperative Congressional 
#'   Election Survey (CCES), paired with predictive data from the 2010 Census. Each 
#'   district's public ideology is paired with a measure of the ideology of the State 
#'   Senate member from the district (update from Shor and McCarty 2011).
#' 
#' @name upperCombined
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @format The \code{upperCombined} dataset has 1989 observations and 10 variables.
#' \describe{
#'   \item{\code{st}}{Two-letter postal abbreviation for the state.}
#'   \item{\code{upper}}{The state legislative district number (upper chamber).}     		
#'   \item{\code{STATEA}}{The FIPS code for the state.}
#'   \item{\code{krige.upper}}{The ideology of the average citizen in the district.}
#'   \item{\code{upperKluge}}{Combined index of \code{STATEA} followed by \code{upper}.}
#'   \item{\code{krige.upper.var}}{The variance of ideology among the district's citizens.}
#'   \item{\code{name}}{Last name of the state legislator, followed by first name and middle initial.}
#'   \item{\code{party}}{Political party of the legislator. D=Democrat, R=Republican, X=Other.}
#'   \item{\code{st_id}}{Temporary identifer variable. DO NOT USE.}
#'   \item{\code{np_score}}{Ideology score for the state legislator (upper chamber). 
#'     Higher values are usually interpreted as more right-wing, with lower values as more left-wing.}
#' }
#' 
#' @source {
#' Ansolabehere, Stephen. 2011. "CCES, Common Content, 2008." Ver. 4.
#' 
#' Minnesota Population Center. 2011. \emph{National Historical Geographic Information 
#'   System: Version 2.0.} Minneapolis, MN: University of Minnesota. \samp{https://www.nhgis.org}
#'   
#'  Shor, Boris and Nolan M. McCarty. 2011. "The Ideological Mapping of American Legislatures." 
#'  \emph{American Political Science Review} 105(3):530-551.
#' }
#' 
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'     \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#'     
#' @examples 
#' # Descriptive Statistics
#' summary(upperCombined)
#' 
#' # Correlate Senators' DW-NOMINATE Scores with Public Opinion Ideology
#' cor(upperCombined$np_score,upperCombined$krige.upper,use="complete.obs")
#' 
#' # Plot Legislators' DW-NOMINATE Scores against Public Opinion Ideology
#' plot(y=upperCombined$np_score,x=upperCombined$krige.upper,
#'      xlab="District Ideology (Kriging)", ylab="Legislator Ideology (Shor & McCarty)", 
#'           main="State Legislatures: Upper Chambers", type="n")
#' points(y=upperCombined$np_score[upperCombined$party=="R"],
#'        x=upperCombined$krige.upper[upperCombined$party=="R"],pch=".",col="red")
#' points(y=upperCombined$np_score[upperCombined$party=="D"],
#'        x=upperCombined$krige.upper[upperCombined$party=="D"],pch=".",col="blue")
NULL

################################### WV WELLS ###################################
# LAST UPDATE: NA

#' West Virginia Oil and Gas Production in 2012
#' 
#' These data are a subset of the West Virginia Geological and Economic Survey of 
#'   2014. They contain information on the coordinates of wells that yielded at 
#'   least some quantity of natural gas in 2012. In addition to coordinates, the 
#'   data contain information on well ownership and operation, rock pressure at 
#'   the well, elevation of the well, oil production, and gas production.
#' 
#' @name WVwells
#' 
#' @docType data
#' 
#' @keywords data
#' 
#' @format The \code{WVwells} dataset has 1949 observations and 18 variables.
#' \describe{
#'   \item{\code{APINum}}{A 10-digit number in the format assigned by the American 
#'     Petroleum Institute (API), consisting of a 2-digit state code, a 3-digit county 
#'     code with leading zeroes, and a 5-digit permit number with leading zeroes. 
#'     Data Source: West Virginia Department of Environmental Protection, Office 
#'     of Oil & Gas (WVDEP-OOG).}
#'   \item{\code{CntyCode}}{A 3-digit numeric code, assigned in numeric order by 
#'     county name. Data Source: The county code for a well is assigned by WVDEP-OOG, 
#'     based on the well location.}
#'   \item{\code{CntyName}}{The name of the county. Please see CntyCode (County Code) 
#'     for a list of all West Virginia county names. Data Source: The county code for 
#'     a well is assigned by WVDEP-OOG, based on the well location. The county name 
#'     is a translation of the county code.}
#'   \item{\code{Operator}}{The name of the operator who owns the well at the time 
#'     of reporting.  Data Source:  WVDEP-OOG plat; verified on the WR-35 completion record.}
#'   \item{\code{SurfaceOwn}}{The name of the owner of the surface property on which 
#'     the well is located. Data Source: WVDEP-OOG plat; verified on the WR-35 completion 
#'     record.}
#'   \item{\code{MineralOwn}}{Mineral Owner: The name of the owner of the mineral 
#'     rights where the well is located. Data Source: WVDEP-OOG plat.}
#'   \item{\code{CompanyNum}}{The operator's serial number for the well. Data Source: 
#'     WVDEP-OOG plat; verified on the WR-35 completion record.}
#'   \item{\code{WellNum}}{The operator's number for the well on the surface property 
#'     (farm). Data Source: WVDEP-OOG plat; verified on the WR-35 completion record.}
#'   \item{\code{UTMESrf}}{Surface Location--Universal Transverse Mercator, Easting: 
#'     The well location at the surface measured in meters to one decimal point, 
#'     east of the central meridian in UTM Zone 17; datum: NAD83. Data Source: 
#'     Taken directly from the plat if given as such. Otherwise, computed from 
#'     the location reported on the plat. Suspect locations may be adjusted using 
#'     various additional resources (e.g. topographic maps) if deemed necessary.}
#'   \item{\code{UTMNSrf}}{Surface Location--Universal Transverse Mercator, Northing: 
#'     The well location at the surface measured in meters to one decimal point, north 
#'     of the equator in UTM Zone 17; datum: NAD83. Data Source: Taken directly from 
#'     the plat if given as such. Otherwise, computed from the location reported on 
#'     the plat. Suspect locations may be adjusted using various additional resources 
#'     (e.g. topographic maps) if deemed necessary.}
#'   \item{\code{LonSrf}}{Surface Location--Longitude: The well location at the 
#'     surface measured to a precision of 6 decimal points, in degrees west of the 
#'     Prime Meridian. Data Source: Taken directly from the plat if given as such. 
#'     Otherwise, computed from the location reported on the plat. Suspect locations 
#'     may be adjusted using various additional resources (e.g. topographic maps) 
#'     if deemed necessary.}
#'   \item{\code{LatSrf}}{Surface Location--Latitude: The well location at the surface 
#'     measured to a precision of 6 decimal points, in degrees north of the equator. 
#'     Data Source: Taken directly from the plat if given as such. Otherwise, computed 
#'     from the location reported on the plat. Suspect locations may be adjusted using 
#'     various additional resources (e.g. topographic maps) if deemed necessary.}
#'   \item{\code{Elevation}}{Elevation: The height of the well in feet above mean 
#'     sea level. Data Source: WVDEP-OOG plat; verified on the WR-35 completion record.}
#'   \item{\code{RockPres}}{Formation Rock Pressure at Surface: The pressure measured 
#'     at the surface usually after stimulation, in pounds per square inch (psi). 
#'     Data Source: WVDEP-OOG WR-35 comp#'   letion record, submitted by the operator 
#'     to WVDEP-OOG.}
#'   \item{\code{GProd2012}}{2012 Gas Production Reported: The total gas production 
#'     for the well for 2012 in thousands of cubic feet (MCF); includes all pay zones. 
#'     Data Source: Production data reported by the operator to the State regulatory 
#'     authority for Oil and Gas (WVDEP-OOG); WVGES obtained the data from WVDEP-OOG.}
#'   \item{\code{OProd2012}}{2012 Oil Production Reported: The total oil production 
#'     for the well for 2012 in barrels (Bbl); includes all pay zones. Production 
#'     data reported by the operator to the State regulatory authority for Oil and 
#'     Gas (WVDEP-OOG); WVGES obtained the data from WVDEP-OOG.}
#'   \item{\code{logElevation}}{Logarithm of \code{Elevation}.}
#' }
#' 
#' @source {
#' 	West Virginia Geological and Economic Survey. 2014. "WVMarcellusWellsCompleted102014." Morgantown, WV. 
#' 	  \url{http://www.wvgs.wvnet.edu/www/datastat/devshales.htm} Accessed via: FracTracker. 2019. 
#' 	  "West Virginia Oil & Gas Activity." \url{https://www.fractracker.org/map/us/west-virginia/}
#' }
#' 
#' @references 
#'   Jason S. Byers & Jeff Gill. N.D. "Applied Geospatial Data Modeling in the Big 
#'     Data Era: Challenges and Solutions."
#'     
#' @examples
#' \dontrun{
#' # Descriptive Statistics
#' summary(WVwells)
#' 
#' # Record means of predictors: 
#' # These are used BOTH to eliminate the intercept and to recover predictions later.
#' mean.logGas<-mean(WVwells$logGProd2012);mean.logGas
#' mean.logElevation<-mean(WVwells$logElevation);mean.logElevation
#' mean.RockPres<-mean(WVwells$RockPres);mean.RockPres
#' 
#' # Outcome Variable, De-Meaned
#' WVwells$logGas <- WVwells$logGProd2012-mean.logGas
#' 
#' # Explanatory Variable: DE-MEANED PREDICTORS AND NO CONSTANT TERM
#' # Because we deducted the mean from all predictors and the outcome,
#' # it is valid to do regression through the origin. 
#' WVwells$LogElevation <- WVwells$logElevation-mean.logElevation
#' WVwells$RockPressure <- WVwells$RockPres-mean.RockPres
#' 
#' # OLS Model
#' fracking.ols<-lm(logGas~LogElevation+RockPressure-1, data = WVwells)
#' summary(fracking.ols)
#' 
#' intercept.mod<-lm(logGProd2012~ logElevation+RockPres,data=WVwells)
#' summary(intercept.mod)
#' 
#' # Set Number of Iterations:
#' # WARNING: 100 iterations is intensive on many machines.
#' # This example was tuned on Amazon Web Services (EC2) over many hours
#' # with 5,000 iterations--unsuitable in 2020 for most desktop machines.
#' #M<-5000
#' M<-20
#' 
#' set.seed(1000, kind="Mersenne-Twister")#SET SEED FOR CONSISTENCY
#' 
#' # Trial Run, Linear Model of Ideology with Geospatial Errors Using Metropolis-Hastings:
#' wv.fit <- metropolis.krige(logGas~LogElevation+RockPressure-1, coords = c("UTMESrf", "UTMNSrf"),
#'                            data = WVwells, n.iter=M, powered.exp=0.5, spatial.share=0.60, 
#'                            range.share=0.31, beta.var=1000, range.tol=0.1, b.tune=1, 
#'                            nugget.tune=1, psill.tune=30)
#' 
#' # Discard first 20% of Iterations as Burn-In (User Discretion Advised).
#' wv.fit <- burnin(wv.fit, M/5)
#' 
#' # Summarize Results
#' summary(wv.fit)
#' 
#' # Convergence Diagnostics
#' # geweke(wv.fit) Not applicable due to few iterations.
#' heidel.welch(wv.fit)
#' 
#' # Draw Semivariogram
#' semivariogram(wv.fit)
#' 
#' # Predictive Data for Two Wells Tapped in 2013
#' well.1<-c(log(1110)-mean.logElevation,1020-mean.RockPres)
#' well.2<-c(log(643)-mean.logElevation,630-mean.RockPres)
#' site.1<-c(557306.0, 4345265)
#' site.2<-c(434515.7, 4258449)
#' well.newdata <- as.data.frame(cbind(rbind(well.1,well.2),rbind(site.1,site.2)))
#' colnames(well.newdata)<-c("LogElevation", "RockPressure", "UTMESrf","UTMNSrf")
#' 
#' # Make predictions from median parameter values:
#' (median.pred <- predict(wv.fit, newdata = well.newdata))
#' 
#' # Prediction in thousands of cubic feet (MCF):
#' round(exp(median.pred+mean.logGas))
#' 
#' # Make predictions with 90\% credible intervals:
#' (cred.pred <- predict(wv.fit, newdata = well.newdata, credible = 0.9))
#' 
#' # Prediction in thousands of cubic feet (MCF) and the true yield in 2013:
#' Actual.Yield<-c(471171, 7211)
#' round(cbind(exp(cred.pred+mean.logGas),Actual.Yield))
#' }
NULL