# Changelog and Todo

## Changelog
#### 12.06.2026 - v1.5.1a
1. Unstable version with parallelization working on testing environment in Linux (Ubuntu v22.04). $\newline$
For a matrix of 9,142 genes x 1,968 samples, previous version on Windows took around 82-98 mins, current version takes approximately 16 mins.$\newline$
$\space$    *NOTE*: `kaleido==0.1.0.post1` was found to be broken in Linux. Linux users should pip install the git version of `kaleido==0.1.0` instead.  
2. `create_interactive_plot` in the `.visualize` module now no longer outputs image to desktop, and instead to local directory.
3. Updated dependencies to more flexible versioning.
4. `sopa` function in `.core` module is less verbose now by default.
5. Improved error handling in `.core.load_sopa`.
6. Migrated to a `src/` layout with `uv` project management and the `hatchling` build backend.

#### 12.06.2026 - v1.1.4
1. Added citation information to README.md.

#### 07.05.2026 - v1.1.3
1. Updated docstrings to accurately detail future iterations.
2. Updated README.md to remove requirements, which are exclusively now found in requirements.txt and pyproject.toml

#### 22.04.2026 - v1.1.2
1. Upgraded pillow to 12.2.0 to resolve security vulnerabilites.
2. Fix expected failures in unit tests
3. Bumped version to 1.1.2, which includes a breaking change ($MSD_{D_{x,s}}$ for RNAseq and DNAm is now not the default metric).
4. Added robust non-parametric $MSD_{D_{x,s}} as the default ranking metric for RNAseq and DNAm data types.

#### 02.04.2026
1. Upgraded pillow to resolve security vulnerabilities.

#### 01.09.2025
1. Updated plotting functions.

#### 18.06.2025
1. Added group trends analysis functions
2. Added new unit tests for new functions


## Todo
#### 26.05.2026
1. Test parallelization on Windows, provide comparison metrics in speed and function.
2. Test new functions to output static images from the `create_interactive_plot` function in the `.visualize` module.


#### 07.08.2026
1. Update MPES functions to return WOCS instead of WCOS to mirror the article.
2. Remove parametric $MSD_{D_{x,s}}$ in future iterations.
3. Add wrapped functions for GMM clustering
4. Update output file names from SOPA and SIMPA.
5. Create readthedocs for thorough documentation.
~~6. Implement parallelization features to improve speed.~~

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
~~6. Introduce example for SIMPA with clinical data.~~