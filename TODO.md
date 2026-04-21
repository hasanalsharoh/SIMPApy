# Changelog and Todo

## Changelog
#### 22.04.2026 - v1.1.2
1. Upgraded pillow to 12.2.0 to resolve security vulnerabilites.
2. Fix expected failures in unit tests
3. Bumped version to 1.1.2, which includes a breaking change ($MSD_{D_{x,s}}$ for RNAseq and DNAm is now not the default metric).
4. Added $MSD_{D_{x,s}^{robust}} as the default ranking metric for RNAseq and DNAm data types.

#### 02.04.2026
1. Upgraded pillow to resolve security vulnerabilities.

#### 01.09.2025
1. Updated plotting functions.

#### 18.06.2025
1. Added group trends analysis functions
2. Added new unit tests for new functions


## Todo
#### 22.04.2026
1. Include OLS for single sample single gene ranking.

#### 06.04.2026
1. Add arguments to sopa.ranking to allow for selection of group columns by adding parameter to replace .startswith('tw') natively for control group selection.
2. Upgrade and test non-static dependencies in higher python package versions.

#### 01.09.2025
1. Allow additional input arguments into different functions (font size, custom columns in .visualize)

#### 03.04.2025
~~1. Add additional error flairs for incorrect files~~
~~2. Fix expected failures in unit tests~~
~~3. Add additional functions for analysis of results obtained through SOPA and SIMPA~~
4. Add single -omic FDR to SIMPA
5. Add filtering based on defined FDR values in visualize module
6. Introduce example for SIMPA with clinical data.