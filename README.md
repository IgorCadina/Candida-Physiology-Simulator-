# CANDIDA ALBICANS PHYSIOLOGICAL SIMULATOR

## GENERAL DESCRIPTION

This project implements an interactive physiological simulator for *Candida albicans*, designed for educational and exploratory purposes. The system heuristically models the microorganism's response to various environmental and nutritional conditions, based on trends widely documented in scientific literature.

The simulator is not intended for diagnostic, clinical, or quantitative experimental purposes. Its goal is to assist in understanding biological and physiological relationships in a qualitative manner.

## PROGRAM OBJECTIVES

* Simulate microbial growth (OD600)
* Evaluate hyphae formation
* Estimate metabolic activity (glycolysis, TCA cycle, ribosomes)
* Evaluate cellular stress and response pathways
* Analyze the impact of pH, temperature, farnesol, and culture medium

## PROJECT ARCHITECTURE

Basic structure:

* **app.py**: Simulator core. Contains the physiological model, decision rules, and the Flask server.
* **index.html**: Web interface for interaction. Allows parameter adjustment and result visualization.
* **response.json**: Example of system output. Not involved in model execution.

## BIOLOGICAL MODEL

The implemented model is semi-realistic (heuristic), utilizing:

* Gaussian curves for pH and temperature response
* Logistic growth kinetics
* Farnesol inhibition thresholds
* Empirical rules for hyphae induction

The general behavior follows trends observed in classic studies, such as:

* Davis et al., 2000
* Hornby et al., 2001
* Lorenz et al., 2000

Numerical values were not adjusted via regression to actual experimental data.

## LIMITATIONS

* Does not quantitatively represent laboratory data
* Does not replace *in vitro* or *in vivo* experiments
* Results should be interpreted only qualitatively

## HOW TO RUN

**Minimum Requirements:**

* Python 3.9 or higher

**Install Dependencies:**

```bash
pip install -r requirements.txt
