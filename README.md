# aneuQPCR
aneuQPCR is a collection of R functions developed for personal use to support the quantification of bulk aneuploidies using qPCR data. Although not structured as a formal R package, this repository can be helpful for anyone who needs to:

Estimate PCR amplification efficiency on a per-reaction basis and without the need for standard curves.

Infer the initial amount of a target DNA template in each reaction without relying on Ct values and taking amplification efficiency into account. 

These functions are particularly useful when studying aneuploidy, where differences in DNA copy number are often subtle (e.g., 1.5× or 2× higher, or 1× lower than a reference chromosome). In such cases, properly compensating for amplification efficiency is crucial to ensure accurate quantification. Although originally developed with aneuploidy detection in mind, the functions are applicable to any qPCR dataset. They can also be useful in other contexts where precise quantification matters, such as gene expression analysis, copy number variation studies, and general molecular diagnostics.

The calculations implemented here are based on the method described in:

Ramakers, C., Ruijter, J. M., Deprez, R. H. L., & Moorman, A. F. M. (2003).
Assumption-free analysis of quantitative real-time polymerase chain reaction (PCR) data.
Neuroscience Letters, 339(1), 62–66. https://doi.org/10.1016/S0304-3940(02)01423-4

Perfect! Here’s a clear and concise “How it works” section based on your explanation:

## How it works

To estimate PCR efficiency and infer initial target DNA concentrations, the following steps are done:

1. **Log2 Transformation**  
   Fluorescence values from each PCR reaction are log₂-transformed to linearize the exponential phase of the amplification curve.

2. **LL5 Model Fitting**  
   A 5-parameter logistic (LL5) model is fit to the log₂-transformed data for each PCR reaction. This model captures the full amplification curve and enables precise identification of curve features.

3. **Inflection Point and Linear Region Detection**  
   Using the LL5 model’s derivatives, the script identifies the linear amplification region and determines the inflection point—where the PCR efficiency begins to decline.

4. **Linear Model Fitting**  
   A linear model is fitted to the data points within the identified linear region, before the inflection point. This portion of the curve reflects the most consistent and efficient phase of DNA amplification.

5. **Efficiency and Initial Concentration Estimation**  
   - The **efficiency** of the PCR is calculated as 2^slope of the fitted linear model.  
   - The **initial amount of target DNA** (at cycle 0) is estimated as 2^intercept, representing the extrapolated fluorescence signal at the start of amplification.

All these steps are carried out by the `fit_model_LL5()` function, which takes as input a `qPCRobj` object. 

---
