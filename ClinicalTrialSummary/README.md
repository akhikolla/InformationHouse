---
output:
  pdf_document: default
  html_document: default
---
# ClinicalTrialSummary

ClinicalTrialSummary provides estimates of several summary measures of the treatment effect for design and analysis of clinical trials with survival outcomes, introduced in Yang (2018). These estimates are obtained under the short-term and long-term hazard ratio model (Yang and Prentice, 2005) which allows a range of time-varying hazard ratio shapes including crossing hazards situations.

## Summary Measures
Let $hr(x) = \lambda_{t}(x)/\lambda_{c}(x)$ be the hazard ratio function, where $\lambda_t(x)$ and $\lambda_c(x)$ are the hazard functions for the treatment group ($t$) and the control group ($c$), respectively.

* The average hazard ratio (AHR): $\int _{0}^{\tau} hr(x) dx$

* The weighted average hazard ratio (WAHR): $\int_{0}^{\tau} hr(x) dw(x)$ where $dw(x) = dF_c(x)/F_c(\tau)$

* The restricted superiority probability ratio (RSPR): $\frac{\int_{0}^{\tau} S_c(x) dF_t(x)}{\int_{0}^{\tau} S_t(x) dF_c(x)}$

* The restricted mean survival difference (RMSD): $\int_{0}^{\tau} S_t(x) dx - \int_{0}^{\tau} S_c(x) dx$

* The ratio of restricted mean times lost (RRMTL): $\frac{\tau - \int_{0}^{\tau} S_t(x) dx}{\tau - \int_{0}^{\tau} S_c(x) dx}$

Here, $\tau$ is the value less than or equal to the maximum follow-up duration of the trial. The asymtoptic results for the average hazard ratio and the restricted mean survival were established in Yang and Prentice (2011) and Yang (2013), respectively. The asymtoptic results for other measures were established in Yang (2018).

## Installation

```
install.packages("ClinicalTrialSummary")
```

## Example

```
library(ClinicalTrialSummary)
data(ggas)
result <- ypsummary(time=ggas$time, event=ggas$event, group=ggas$group, tau=8.2)
result
```
The data "ggas" is from Gastrointestinal Tumor Study Group (1982) and the value for `tau` must be user-specified. The object `result` can be formatted to a table using the function `summary`.

```
summary(result)
```

## Reference
Yang, S. (2018). Improving testing and description of treatment effect in clinical trials with survival outcomes. Statistics in medicine.

Yang S, and Ross L. Prentice (2005). Semiparametric analysis of short-term and long-term hazard ratios with two-sample survival data. Biometrika, 92.1:1-17.

Yang S, and Ross L. Prentice (2011). Estimation of the 2-sample hazard ratio function using a semiparametric model. Biostatistics, 12.2:354-368.

Yang S. (2013). Semiparametric inference on the absolute risk reduction and the restricted mean survival difference in clinical trials. Special issue on risk assessment. Lifetime Data analysis, 19:219-241.

Gastrointestinal Tumor Study Group (1982). A comparison of combination chemotherapy and combined modality therapy for locally advanced gastric carcinoma. Cancer.
