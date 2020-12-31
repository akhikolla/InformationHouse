\name{wine}
\alias{wine}
\docType{data}
\title{
Wine data set
}
\description{
This data set is made of 178 observations (Italian wines)  described by 27 variables (physicochemical measures).
These wines come from from three different regions of Italy.
}
\usage{data("wine")}
\format{
  We have labels and data as follows :
 
  The last column of the data frame (wine[,28]): it indicates the class label 1,2 or 3.
  
  The data involving columns 1 to 27: 
   
   Alcohol                      
  
   Sugar-free_extract           
  
   Fixed_acidity               
  
   Tartaric_acid                
  
   Malic_acid                   
  
   Uronic_acids                
 
   pH                           
  
   Ash                          
  
   Alcalinity_of_ash           
  
   Potassium                    
  
   Calcium                      
  
   Magnesium                   
  
   Phosphate                    
  
   Chloride                     
  
   Total_phenols               
  
   Flavanoids                   
  
   Nonflavanoid_phenols         
  
   Proanthocyanins             

   Color_Intensity              
  
   Hue                          
  
   OD280/OD315_of_diluted_wines

   OD280/OD315_of_flavanoids    
  
   Glycerol                     
  
   2-3-butanediol              

   Total_nitrogen               
  
   Proline                      
  
   Methanol                    
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(wine)

head(wine)
}
\keyword{datasets}
