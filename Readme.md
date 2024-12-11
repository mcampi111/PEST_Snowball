# Optimizing Snowball Referrals in Expert Surveys
This repository is linked to the methodology developed in the paper with title

**"Optimizing Snowball Referrals in Expert Surveys"**. The pdf for the paper is available at this url https://.. or provided in the folder "paper", where it is also possibile to find the Supplementary Materials.

## **Abstract**

Expert opinion is crucial for decision-making across many fields, from policy development to scientific assessments. However, identifying and recruiting appropriate experts presents significant challenges, particularly due to selection biases and limited visibility of expertise. This paper introduces a quantitative framework for optimizing snowball sampling—a method where experts help recruit other experts—by analyzing how different network structures affect sampling success. Through extensive simulations, we determine optimal numbers of initial experts (seeds) and recruitment waves needed to achieve representative samples. Our findings provide practical guidelines for researchers seeking expert opinions, showing how sampling efficiency varies with network structure and recruitment rates. This methodology helps minimize bias in expert selection and can be applied beyond expert surveys to other hard-to-reach populations, improving the reliability of expert-based research. 


## Contributions of the paper
The contribution of this work lies in three interrelated methodological advances aimed at defining a population and recruiting a representative sample, particularly when the population is partially hidden and small, and recruitment involves high costs in terms of time and resources. These are given as follows:
1. The first contribution is to treat expertise as a property of a community of experts rather than as an individual attribute. In stage one of our methodology, we construct a snowball sample from an initial population of identified experts, delegating the identification of expertise to those within the community itself. If feasible, multiple snowball waves are used to rank levels of expertise within the community. The success of this stage relies on identifying a sufficient number of gatekeepers (or seeds) who are well-connected within their field.
2. The second contribution is the use of snowball referrals to facilitate expert recruitment. Experts are recruited based on peer referrals, which are a reflection of their esteem within the community. The higher the number of referrals an individual receives, the more highly esteemed they are considered by their peers. This approach allows for a natural prioritization of experts who are recognized for their contributions, while also depending on the provision of appropriate incentives to continue recruitment. In our simulations, we assume a consistent response rate within the sampled population, though this may vary in practical settings.
3. The third contribution focuses on determining the optimal point at which to conclude data gathering. In the context of quantitative research, a representative sample is sought, while qualitative research aims for data saturation—the point at which no new information emerges. We adopt a mixed-methods approach, aiming for a level of sampling that ensures a sufficient portion of the expert population has been covered to achieve what we call pseudo-representativeness. This is crucial for ensuring the collected data is sufficiently rich while avoiding unnecessary costs associated with prolonged data collection.
4. 
## Motivations For the Development of the Paper

TOADD..

## Methodology of the Paper

TOADD

## Organization of the Repository
TOADD

```diff
+ 1) R_code 
```
This folder contains the R code developed for the analysis of the paper. This is structured within the following folders:

1.  **blabla**. This folder contains blabla.
2. **blabla**. This folder contains blabla.
3.  **blabla**. This folder contains blabla.
4.  **blabla**. This folder contains blabla.

```diff
+ 2) paper 
```
This folder contains the draft of the paper which is also available at https://...

```diff
+ 3) paper_figures 
```
This folder contains all the paper figures of the paper. For each figure produced, there is an Rcode associated within the R_code folder.


## Cite

If you use this code in your project, please cite:

@article{blablabla,
  title={Optimizing Snowball Referrals in Expert Surveys},
  author={Christopolou, Dimitris and Jose, Alex and Campi, Marta},
  journal={Available at ..},
  year={2024}
}

