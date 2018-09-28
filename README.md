# Calculations for sub-GeV dark matter
--------------------------------------

This repository contains tools for computing amplitudes for models of sub-GeV dark matter that have been matched onto the chiral Lagrangian. Notes about the different files follow.

## `sub_GeV_DM_EFT_calcs/Spectra/charged_pion_decay.nb`

This notebook computes the contribution to the charged pion radiative decay spectrum coming from ISR off the pion, FSR off the muon or electron and IB from the pi-l-nu vertex.

## `sub_GeV_DM_EFT_calcs/`

Various FeynRules (FR) model files. These are definitions of models at the Lagrangian level, which can be processed into FeynArts (FA) model files. The FA models also contain code to copy the model file to the relevant Mathematica subdirectory: this may require changing `targetDir`.

### `sub_GeV_DM_EFT_calcs/FeynRules models/Chiral perturbation theory`

FR models for chPT including the rho and with only pions and the rho. Running `write_model.nb` generates the FA models found in `ChiPT/` and `ChiPT_expanded_pi_only/`.

### `sub_GeV_DM_EFT_calcs/FeynRules models/EFT of MeV DM`

Each of the subdirectories contains the FR model files and a `write_<mediator>_model.nb` file to generate the FA models. The notebook `correct_scalar_vev_calculation.nb` aims to compute the scalar's vev.

## `sub_GeV_DM_EFT_calcs/Simplified model amplitudes/`

Contains a directory for each mediator with a file to compute DM self-annihilation cross sections, mediator decay widths, and some FSR spectra, and to write them to files that can easily be reformatted into python functions.

TODO: Logan, could you put the FSR spectra in here?


