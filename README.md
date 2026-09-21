# Statistical patterns in baby name choice

## Paper
This repository contains code used for both versions of the age-structured neutral model described in the manuscript:

'Innovation, conservatism, and the need to be distinguishable shape the diversity of first names' by Anne Kandler,  Rafael D'Andrea, James O'Dwyer

> *Abstract*: While neutral models provide influential baselines in cultural evolution, identifying the mechanisms underlying observed deviations remains challenging. Here, we investigate these deviations using first-name data from multiple Western populations, focusing on two complementary measures: the variant abundance distribution (VAD, name diversity among living individuals) and the progeny distribution (PD, name diversity among newborns over time). Standard neutral theory predicts a power-law VAD with exponent of -1, yet 1930 US census data exhibit exponents of roughly -1.7 with elevated fractions of both rare and very common names. We show that age-constrained cultural transmission (preferential copying of recently-transmitted variants) combined with anti-novelty bias (preference for established over novel names) quantitatively reproduces these patterns. Despite these deviations in the VAD, we find that thresholded progeny distributions across eight datasets are well-described by an effective neutral model. Fitting this model reveals a striking scaling relationship: larger populations exhibit lower per-capita effective innovation rates. We explain this through a functional constraint: names must distinguish individuals within local social networks, not entire populations. Implementing this anti-dominance bias in our age-structured model reproduces the observed scaling. Our results demonstrate how age structure, cultural biases, and functional constraints interact to shape cultural diversity patterns.

## Code

All code is implemented in Matlab. The script `main_ageSim.m` can be used to
generate populations using the simulation models described in the paper. 


## License
[![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg

