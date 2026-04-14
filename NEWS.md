## [0.2.2] - 2026-04-14

### Fixed

- Fix bug in function `AUCuniform.2` when ties occur at maximum value of the classification rule.

### Documentation

- Minor improvements and typo fixes.

## [0.2.1] - 2026-01-08

### Documentation

- Cleaned documentation: deprecated function `AUCuniform_trap()` is no longer indexed or shown in the PDF manual. The function is kept for backward compatibility.

- Updated the reference for Jiménez-Valverde (2025) in `AUCuniform.2.Rd` to include the journal name, volume number and article number.

- Clarified documentation of `AUCuniform()`: added note that it is retained for completeness and reproducibility, while `AUCuniform.2()` is recommended for most applications.

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

