## [0.2.0] - 2025-10-24

### Added

- Added a new `AUCuniform.2()` function to compute the uniform AUC (uAUC) and uniform sensitivity-star (uSe*) using the direct weighted trapezoidal estimation method. This function correctly deals with ties (instances of presence and absence with the same probability values), which were not correctly handled in `AUCuniform_trap()`.

### Deprecated

- `AUCuniform_trap()` is now deprecated. Please use `AUCuniform.2()` instead. This function will be removed in a future version.

## [0.1.0] - 2024-09-02

### Added

- Added a new `HSgraph()` function to visualize the distribution of suitability values for presences, absences, and all cases combined.

- Added a new `AUCuniform_trap()` function to compute the uniform AUC (uAUC) using the weighted trapezoidal method. 

### Documentation

- Updated the reference for Jiménez-Valverde (2022) in `AUCuniform.Rd` to include the volume number and page range.

