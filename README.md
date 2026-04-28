# Statistical patterns in baby name choice

## Paper
This repository contains code used for both versions of the age-structured neutral model described in the manuscript:

'A balance of innovation, conservatism, and recent popularity explains statistical patterns in baby name choice' by Anne Kandler,  Rafael D'Andrea, James O'Dwyer

> *Abstract*: While neutral models have become influential baselines in cultural evolution, translating observed deviations from these baselines into specific underlying mechanisms remains challenging. Here, we investigate these deviations using first-name data from multiple Western populations, focusing on two complementary measures: the variant abundance distribution (VAD), which describes name diversity among living individuals, and the progeny distribution (PD), which characterizes names given to newborns over time. Standard neutral theory predicts a power-law VAD with exponent $\approx-1$, yet empirical data from the 1930 US census exhibit exponents $\approx−1.7$, alongside elevated fractions of both rare and common names. We show that age-constrained cultural transmission,where individuals preferentially copy recently-transmitted variants, combined with anti-novelty bias, a preference for established names over innovations, can quantitatively reproduce these patterns. Despite these deviations in the VAD, we find that thresholded progeny distributions across eight datasets are well-described by an effective neutral model. Fitting this model reveals a striking scaling relationship: larger populations exhibit lower per-capita effective innovation rates. We explain this inverse scaling through a functional constraint: names must distinguish individuals within local social networks, but not across entire populations. Implementing this as an anti-dominance bias in our age-structured model successfully reproduces the observed scaling. Our results demonstrate how age structure, cultural biases, and functional constraints interact to shape patterns of cultural diversity, and suggest that the function of cultural variants may serve as an organizing principle for understanding cultural change.


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

