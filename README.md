# Detecting Manipulation of Poverty Scores in Colombia's SISBEN

Stata code for my master's seminar paper in *Social Policy Evaluation in Latin America* (Freie Universität Berlin, 2026).

## Question

Colombia's welfare programmes use a poverty index score (SISBEN) to decide who is eligible: households scoring 47 or lower qualify. Camacho & Conover (2011) showed that after the cutoff became public in 1997–98, unusually many households appeared just below it, a sign that scores were manipulated.

This project replicates their discontinuity estimates and tests them with a modern density test.

## What the code does

1. **Sample construction:** rebuilds the score from its components, flags and removes inconsistent records (~0.03%), and collapses the data to the family level.
2. **Replication:** re-estimates the authors' yearly discontinuity at the cutoff (their Table 3) for individuals and families.
3. **Extension:** applies the Cattaneo–Jansson–Ma (CJM) density test (`rddensity`) per year, with automatic and fixed bandwidths, plus placebo cutoffs at 42 and 54.
4. **Municipality level:** compares the authors' bin-difference measure with CJM estimates for each municipality–year.
5. **Figures and tables:** score histograms, CJM estimates over time, and a combined results table.

## Results

<!-- Add 1–2 sentences with your main finding here, e.g. the size of the jump after 1998 -->

![CJM estimates by year](figures/fig2_cjm_by_year.png)

## How to run

1. Download the replication data for Camacho & Conover (2011), *American Economic Journal: Economic Policy* 3(2), from the [AEA website](https://www.aeaweb.org/articles?id=10.1257/pol.3.2.41) and put `sisben_aejep.dta` into `data/`.
2. Open Stata 17 or newer in the repository folder and run `replication.do`. Missing packages (`rddensity`, `lpdensity`, `estout`) install automatically.
3. Slow parts can be switched on or off at the top of the file (`RUN_A`, `RUN_YEARLY`, `RUN_CELLS`).

The data is not included in this repository.

## Reference

Camacho, A. & Conover, E. (2011). Manipulation of Social Program Eligibility. *American Economic Journal: Economic Policy*, 3(2), 41–65.
