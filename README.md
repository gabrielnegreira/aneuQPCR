# aneuQPCR
aneuQPCR is a collection of R functions developed for personal use to support the quantification of bulk aneuploidies using qPCR data. Although not structured as a formal R package, this repository can be helpful for anyone who needs to:

Estimate PCR amplification efficiency on a per-reaction basis and without the need for standard curves.

Infer the initial amount of a target DNA template in each reaction without relying on Ct values and taking amplification efficiency into account. 

These functions are particularly useful when studying aneuploidy, where differences in DNA copy number are often subtle (e.g., 1.5× or 2× higher, or 1× lower than a reference chromosome). In such cases, properly compensating for amplification efficiency is crucial to ensure accurate quantification. Although originally developed with aneuploidy detection in mind, the functions are applicable to any qPCR dataset. They can also be useful in other contexts where precise quantification matters, such as gene expression analysis, copy number variation studies, and general molecular diagnostics.

The calculations implemented here are based on the method described in:

Ramakers, C., Ruijter, J. M., Deprez, R. H. L., & Moorman, A. F. M. (2003).
Assumption-free analysis of quantitative real-time polymerase chain reaction (PCR) data.
Neuroscience Letters, 339(1), 62–66. https://doi.org/10.1016/S0304-3940(02)01423-4
