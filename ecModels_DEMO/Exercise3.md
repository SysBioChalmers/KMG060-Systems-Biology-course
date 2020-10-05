# Excercise 3: FBA simulations with enzyme constrained models
The predictive capabilities of a purely metabolic model have been explored with yeastGEM in [Exercise#2](https://github.com/SysBioChalmers/workshops/blob/master/Tartu_2020/instructions/Excercise2.md) in which the following limitations were identified:
1. A large number of unbounded reaction fluxes are obtained in a typical FBA simulation when just a few constraints are imposed.
2. GEMs predict a linear relation between maximum growth rates and glucose uptake rate, which might lead to erroneous predictions if high uptakes rates are provided.

Enzyme-constrained models of metabolism overcome these limitations by the addition of kinetic constraints to the reactions in a metabolic network together with a total protein pool that limits the amount of enzyme mass available for catalyzing such reactions. In this exercise the predictive advantages of ecModels are explored.

1. Load both `yeastGEM` and `ecYeastGEM` models from the `models` subfolder in this repo. 

The `ecYeastGEM` is an enzyme-constrained model without any incorporated proteomics data, therefore all of its enzymes are conected and constrained by the total protein pool which is a function of the total protein content of the modelled cell.

2. Get the total number of promiscuous enzymes isoenzymes and rxns with Kcats present in `ecYeastGEM`.

3. Obtain a FBA solution for both `yeastGEM` and `ecYeastGEM` models with *growth rate* as an objective for maximization subject to a unit glucose uptake rate [1 mmol/gDw h]. Show both flux distributions as cumulative distributions.

- Explain the differences between both cumulative distributions (in terms of median and extreme values).

4. Get a GUR vs growth rate plot for both `yeastGEM` and `ecYeastGEM` models in the range 1-20 [mmol/gDw h] of glucose uptake rate.

a) Which of the obtained plots agrees better with experimental evidence?
b) Does any of the plots reach a saturation point? If so, explain why in biological terms.

5. Simulate batch growth on diverse carbon sources, with unconstrained carbon uptake, for:
a) Minimal media
b) Minimal media + aminoacids
c) Complex YPD media

In the case of minimal media simulations calculate the biomass yield in terms of mass [g biomass / g carbon source]. To facilitate this task use the scripts `changeMedia_batch` and `Csources_simulations` provided in the `scripts/complementaryScripts` subfolder.

Compare your results with the experimental data provided in `data/growthRates_data_carbonSources.txt`.

6. Find the critical dilution rate at which **S. cerevisiae** switchs from pure respiration to fermentation (Crabtree effect) when growing on glucose minimal media.







