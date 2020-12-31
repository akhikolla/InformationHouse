\name{winratio}
\alias{winratio}
\title{Win Ratio for Prioritized Outcomes and 95\% Confidence Interval}
\usage{
winratio(id, trt, active = NULL, outcomes, fu, data, keep.matrix = FALSE)
}
\arguments{
\item{id}{
a string indicating the patient ID variable. The patient ID variable must not contain missing values or duplicates.
}

\item{trt}{
a string indicating the treatment variable. The variable can be of any type (binary, numeric, character, factor) and must imperatively have 2 unique values/levels.
}

\item{active}{
a numeric or string value used to define the active treatment group for the calculation of the win ratio. If \code{active = NULL} (default value), active group is automatically determined according to the type of treatment variable:
  \itemize{
    \item binary variable: active group = 1.
    \item numeric variable: active group = maximum value.
    \item character variable: active group = last string after sorting the character variable in ascending order.
    \item factor variable: active group = last level.
  }
}

\item{outcomes}{a list used to define all outcomes in order of priority. 
The first element to be defined must be the outcome considered as being of greater clinical importance, usually this outcome is a fatal event. 
Each element must be a character vector of length 3 or a character list of length 3 depending on the type of outcome (survival event, repeated survival event, continuous or ordinal event): 
  \itemize{
    \item{
    for survival events, the element is a character vector \code{c(x, "s", y)}, where \code{x} is a string indicating the name of the event variable, \code{"s"} indicating that this outcome is a "failure-time" event, and \code{y} is a string indicating the name of the time-to-event variable.
    } 
    \item{
    for repeated survival events, note \code{k} the maximum number of events observed during follow-up, the element is a character list \code{list(x, "r", y)}, where \code{x = c(x1, ..., xk)} is a character vector of length \code{k} indicating the name of event variables, \code{"r"} indicating that this outcome is a "failure-time" repeated event, and \code{y = c(y1, ..., yk)} is a character vector of length \code{k} indicating the name of time-to-event variables corresponding to the event variables.
    }
    \item{
    for continuous or ordinal events, the element is a character vector \code{c(x, "c", d)}, where \code{x} is a string indicating the name of the continuous or ordinal variable, \code{"c"} indicating that this outcome is a continuous or ordinal event, and \code{d} is a string with two options \code{"<"} or \code{">"} indicating the direction of deterioration or worsening.
    }
  }
}

\item{fu}{
a string indicating the name of the follow-up time variable.
}

\item{data}{
a data frame containing all the variables listed in the \code{id}, \code{trt}, \code{outcomes} and \code{fu} arguments.
}

\item{keep.matrix}{
a logical value indicating if the 'win-loss' matrix is kept. Default value is \code{FALSE}.
}

}

\value{
\item{call}{
a list with all arguments entered \code{list(id, trt, active, outcomes, fu, data, keep.matrix)}.
}
\item{group1}{
the value/level of the active group (group 1).
}
\item{group0}{
the value/level of the other group (group 0).
}
\item{n1}{
the number of subjects in group 1.}
\item{n0}{
the number of subjects in group 0.}
\item{n}{
the total number of subjects.}
\item{wins}{
a vector of numbers of 'winners' for each outcome.
}
\item{loss}{
a vector of numbers of 'losers' for each outcome.
}
\item{total.wins}{
the total number of 'winners'.
}
\item{total.loss}{
the total number of 'losers'.
}
\item{total.ties}{
the total number of ties.
}
\item{wr}{
the win ratio.
}
\item{v}{
the estimated variance of win ratio.
}
\item{z}{
the value of the test statistic.
}
\item{p.value}{
the p-value of the test.
}
\item{wr.lower}{
the lower end of the 95\% confidence interval.
}
\item{wr.upper}{
the upper end of the 95\% confidence interval.
}
\item{wr.matrix}{
the 'win-loss' matrix (only if \code{keep.matrix = T}). See details above.
}
}

\details{
\subsection{'win-loss' matrix}{
Let \code{n1} the number of patients in the active group and \code{n0} the number of patients in the other group. 
The 'win-loss' matrix is a matrix \code{M} of dimension \code{n1 x n0} with its element \code{M[i,j]} is defined as:

\code{M[i,j] = 1} if \code{i} wins against \code{j} on 1st outcome,

\code{M[i,j] = -1} if \code{i} loses against \code{j} on 1st outcome,

\code{M[i,j] = k} if \code{i} wins on outcome \code{k} after ties with \code{j} on outcomes \code{1} to \code{(k-1)},

\code{M[i,j] = -k} if \code{i} loses on outcome \code{k} after ties with \code{j} on outcomes \code{1} to \code{(k-1)},

\code{M[i,j] = 0} if \code{i} and \code{j} ties on all outcomes.
}
}

\description{
Calculate the win ratio for prioritized outcomes and the 95\% confidence interval based on Bebu and Lachin (2016). Three type of outcomes can be analyzed: survival "failure-time" events, repeated survival "failure-time" events and continuous or ordinal "non-failure time" events that are captured at specific time-points in the study.
}

\examples{
# Load survival package to use the dataset 'bladder1'
library(survival)

# Creation of dataset 'df' with 3 outcomes:
#   Outcome 1: death (survival event)
#   Outcome 2: cancer recurrence (repeated survival event)
#   Outcome 3: size of largest initial tumour (continuous event) 
data1 <- bladder1  \%>\%  
  mutate(trt = if_else(treatment == "placebo", "Placebo", "Treatment")) \%>\%
  group_by(id) \%>\% 
  mutate(death = if_else(max(status) \%in\% c(2, 3), 1, 0),
         t2death = max(stop)) \%>\% 
  ungroup() \%>\%
  select(id, trt, death, t2death, number, size) \%>\% 
 unique()

data2 <- bladder1 \%>\% 
  filter(status == 1) \%>\% 
  select(id, t2recurr = stop) \%>\%
  mutate(recurr = 1) \%>\% 
  arrange(id, t2recurr) \%>\% 
  group_by(id) \%>\% 
  mutate(nrecurr = row_number()) \%>\% 
  ungroup() \%>\% 
  full_join((data1 \%>\% select(id)), by = "id") \%>\% 
  complete(id, nrecurr) \%>\%
  filter(!is.na(nrecurr)) \%>\% 
  full_join((data1 \%>\% select(id, t2death)), by = "id") \%>\% 
  mutate(
    recurr = replace(recurr, is.na(recurr), 0),
    t2recurr = if_else(is.na(t2recurr), t2death, t2recurr)
  ) \%>\% 
  select(id, nrecurr, recurr, t2recurr)

data3 <- data2 \%>\% 
  pivot_wider(id_cols = "id", names_from =  nrecurr, 
              values_from = recurr, names_prefix  = "recurr") 

data4 <- data2 \%>\% 
  pivot_wider(id_cols = "id", names_from =  nrecurr, 
              values_from = t2recurr, names_prefix  = "t2recurr") 

data5 <- full_join(data3, data4, by = "id")

df <- full_join(data1, data5, by = "id")

# Calculate the win ratio
wr <- winratio(id = "id", trt = "trt", active = "Treatment", 
               outcomes = list(outc1 = c("death", "s", "t2death"),
                               outc2 = list(paste0("recurr", 1:9), "r", paste0("t2recurr", 1:9)),
                               outc3 = c("size", "c", ">")), 
               fu = "t2death", data = df)
summary(wr)
}
