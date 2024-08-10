This repository includes R code and the poster presentation for the project titled "Bayesian Time-Dynamic Modeling of Post-Transplant Kidney Function in Chronic Kidney Disease Patients." The project was supported by the DS-PATH Summer Fellowship Program under the National Science Foundation Harnessing Data Revolution Data Science Corps Award #2123444, #2123271, #2123313.

**Team Members:** Liam Daly (UCR, Statistics), Anthony Gutierrez (Cal State San Bernardino, Computer Science), Mubarizuddin Mohammed (UCR, Computer Science), and Ali Syed (UCR, Data Science).  
**Mentor:** Dr. Esra Kurum (UCR, Statistics).

**Project Description:** The team analyzed data from chronic kidney disease (CKD) patients who underwent kidney transplantation, with the goal of modeling the progression of kidney function post-transplant and identifying risk factors that influence outcomes in patients with successful versus failed transplants. Time-dynamic models were employed for the analysis, and estimation was performed using Bayesian P-splines. The models were fit using the Bayesian software JAGS (version 4.3.0) via the R2jags package.

**Files in the Repository:**  
- `ckd.rdata`: Publicly available data from CKD patients who underwent a primary renal transplantation with a graft from a deceased or living donor at the University Hospital of the Catholic University of Leuven (Belgium) between January 21, 1983, and August 16, 2000.
- `CKD_data_description.pdf`: Detailed information on the dataset.
- `GFR_BayesianTDM_survived.R`: R code for fitting time-dynamic models to patients with successful transplants.  
- `GFR_BayesianTDM_failed.R`: R code for fitting time-dynamic models to patients with failed transplants.  
- `Poster.pdf`: Comprehensive project overview, including results and future directions.


