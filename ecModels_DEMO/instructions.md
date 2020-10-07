# ecModels DEMO: FBA simulations with enzyme constrained models
The predictive capabilities of a purely metabolic model have been explored in **Exercise#3** in which the following limitations were identified:
1. A large number of unbounded reaction fluxes are obtained in a typical FBA simulation when just a few constraints are imposed.
2. GEMs predict a linear relation between maximum growth rates and glucose uptake rate, which might lead to erroneous predictions if high uptakes rates are provided.

Enzyme-constrained models of metabolism overcome these limitations by the addition of kinetic constraints to the reactions in a metabolic network together with a total protein pool that limits the amount of enzyme mass available for catalyzing such reactions. In this exercise the predictive advantages of ecModels are explored.

* The execution of the live script for this DEMO relies on some functions available in the **RAVEN toolbox**, the details for its installation can be found [here](https://github.com/SysBioChalmers/RAVEN/wiki)

1. Load both `yeastGEM` and `ecYeastGEM` models from the `models` subfolder in this repo. Explore the differences in their number of reactions, metabolites, genes, etc.

The `ecYeastGEM` is an enzyme-constrained model without any incorporated proteomics data, therefore all of its enzymes are conected and constrained by the total protein pool which is a function of the total protein content of the modelled cell.

2. Get the total number of promiscuous enzymes isoenzymes and rxns with Kcats present in `ecYeastGEM`.

3. Obtain a FBA solution for both `yeastGEM` and `ecYeastGEM` models with **growth rate** as an objective for maximization subject to a unit glucose uptake rate [1 mmol/gDw h]. Show both flux distributions as cumulative distributions.

- Explain the differences between both cumulative distributions (in terms of median and extreme values).

4. Get a GUR vs growth rate plot for both `yeastGEM` and `ecYeastGEM` models in the range 1-18 [mmol/gDw h] of glucose uptake rate.

a) Which of them looks more "realistic"?
b) Does any of the plots reach a saturation point? If so, explain why in biological terms.

5. Simulate batch growth on diverse carbon sources, with unconstrained carbon uptake, for:
a) Minimal media
b) Minimal media + aminoacids
c) Complex YPD media

6. Find the critical dilution rate at which *S. cerevisiae* switchs from pure respiration to fermentation (Crabtree effect) when growing on glucose minimal media. Assume growth on chemostats under the assumption that cells grow as efficiently as possible (in carbon usage terms) with the least amount of protein burden in order to adapt to any sudden perturbation.







