# Copula CHMM – Variational EM Algorithm Implementation in R

This repository contains an R implementation of a **Copula-based Coupled Hidden Markov Model (CHMM)** using a **Variational Expectation-Maximization (VEM)** algorithm.  
The model incorporates copula functions to capture dependence between multiple diseases (or variables) while estimating hidden state transitions and emission distributions.

---

## Features

- Initialization of CHMM parameters (transition, emission, and initial state probabilities).
- Variational EM algorithm for model parameter estimation.
- Log-forward and log-backward probability computation for numerical stability.
- **Binomial copula** integration for dependency modeling across dimensions. Or Uses **geometric copula functions** to model dependencies between multiple disease processes 
- Automatic calculation of AIC and BIC for model selection.
- Poisson emission distribution (adaptable for other types).

---

## Repository Structure

├── CHMM_VEM for Poisson.R # Main implementation of the Copula CHMM Binomial 

├── CHMM_VEM for Poisson Geometric.R # Main implementation of the Copula CHMM Geometric 

├── Init_Param.R #Parameter Initialization for the model 

├── Forward-Backward probabilities.R #to calculate Forward-backward

├── Data set up.R #To set up dataset for the model

├── Data.csv #Dataset example

├── README.md # Project documentation



---

## Requirements

- R (≥ 4.0)
- Required packages:
  - [`Rmpfr`](https://cran.r-project.org/package=Rmpfr)
  - [`dplyr`](https://cran.r-project.org/package=dplyr)

You can install them using:

```R
install.packages("Rmpfr")
install.packages("dplyr")

```
## Example Usage 

```R
# Load required libraries
library(Rmpfr)
library(dplyr)

# Example: Assuming you have data matrix X, number of states, and omega (copula parameters)
#In our case dataset up set these according to the dataset
nb.states <- 3
nbD <- ncol(X)            # Number of dimensions
omega <- matrix(1.5, nbD, nbD)  # Example copula matrix

# Run the model
mod <- CHMM_VEM(
X=X,
nb.states=3,
itmax = 100,
threshold = 0.00001,
nbD=nbD,
omega=omega
)

# View results
print(mod)

```
## Output includes:

-Posterior probabilities (lpostPr)

-Estimated emission parameters (esLambda)

-Transition matrix (transPr)

-Model fit criteria (AIC, BIC)

