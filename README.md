
# Resazurin Assay Development

This repository contains all code, data, and documentation for the development and application of the resazurin assay to measure oyster metabolic responses to environmental stress. The project integrates multiple experimental datasets, analysis scripts, and outputs, supporting research in oyster health, aquaculture, and genetics.

See our [landing page](https://robertslab.github.io/resazurin-assay-development/) for information on resazurin protocols, recent findings, experimental design tools, and data explor apps! 

## Repository Structure

```
resazurin-assay-development/
│
├── _quarto.yml                # Quarto website configuration
├── index.qmd                  # Technical documentation landing page
├── public-summary.qmd         # Public summary for general audiences
├── README.md                  # Project overview
├── resazurin-assay-development.Rproj # RStudio project file
│
├── data/                      # Raw and processed data for all projects
│   ├── 10k-seed/              # Data for 10K Seed Project
│   ├── oxygen-correlation/    # Data for oxygen-resazurin correlation
│   ├── spat-stress/           # Data for Spat Stress Project
│   ├── testing/               # Data for method validation/testing
│   ├── thermal-curve/         # Data for thermal curve project
│   └── usda-families/         # Data for USDA Families Project
│
├── scripts/                   # RMarkdown analysis scripts
│   ├── 10k-seed/              # Analyses for 10K Seed Project
│   ├── oxygen-correlation/    # Oxygen-resazurin analyses
│   ├── spat-stress/           # Analyses for Spat Stress Project
│   ├── testing/               # Method validation scripts
│   ├── thermal-curve/         # USDA thermal curve analyses
│   └── usda-families/         # USDA Families analyses
│
├── figures/                   # Figures and diagrams
│   ├── resazurin-applications-diagram.md # Visual schematic
│   └── ...                    # Project-specific figures
│
├── output/                    # Analysis outputs (tables, processed data)
│   ├── 10k-seed/
│   ├── oxygen-correlation/
│   ├── spat-stress/
│   ├── thermal-curve/
│   └── usda-families/
│
├── docs/                      # Quarto website output
│   ├── index.html
│   ├── public-summary.html
│   └── figures/
│
└── ...                        # Additional documentation and metadata
```

## Project Connections

- **10K Seed Project**: Examines oyster seed response to thermal stress, connecting acute temperature effects to resazurin measurements and mortality. [Related repo](https://github.com/RobertsLab/10K-seed-Cgigas)
- **Spat Stress Project**: Explores thermal stress response in spat. [Related repo](https://github.com/RobertsLab/polyIC-larvae)
- **USDA Families Project**: Investigates genetic variation in stress response among oyster families, comparing absorbance and fluorescence modes.
- **Oxygen Correlation**: Links resazurin assay results to direct oxygen consumption measurements for validation.
- **Thermal Curve**: Analyzes metabolic responses across temperature gradients in USDA families.
- **Testing**: Contains method validation and protocol optimization scripts.

## Purpose of Key Scripts

- `scripts/10k-seed/resazurin-analysis-10k-seed.Rmd`: Analyzes metabolic rate data for 10K Seed Project.
- `scripts/10k-seed/survival-analysis-10k-seed.Rmd`: Survival analysis based on resazurin data.
- `scripts/oxygen-correlation/oxygen-rate-extraction.Rmd`: Extracts respiration rates from oxygen data.
- `scripts/oxygen-correlation/oxygen-resazurin-correlations.Rmd`: Correlates oxygen and resazurin measurements.
- `scripts/oxygen-correlation/resazurin-analysis.Rmd`: Analyzes individual oyster resazurin data for oxygen comparison.
- `scripts/spat-stress/resazurin-analysis-spat-stress.Rmd`: Analyzes spat stress project data.
- `scripts/testing/shell-vs-open.Rmd`: Compares empty shell vs live animal metabolic rates.
- `scripts/thermal-curve/thermal-curve-analysis.Rmd`: Analyzes thermal curve data for USDA families.
- `scripts/usda-families/resazurin-analysis-usda-families.Rmd`: Analysis of metabolic rates between USDA families.

## Data Organization

- Each project folder in `data/` contains raw data files (CSV, XLSX) and metadata.
- `output/` contains processed data and results from analysis scripts.
- `figures/` contains visualizations and diagrams for each project.

## Documentation

- `index.qmd`: Technical documentation and protocol details.
- `public-summary.qmd`: Lay summary for non-specialists.
- `_quarto.yml`: Website configuration for navigation and theming.

## Visual Schematic

Below is a simplified flow diagram from `figures/resazurin-applications-diagram.md`:

```
OYSTER SAMPLE
	│
	▼
RESAZURIN ASSAY
 • Add dye
 • Apply stress
 • Measure signal
	│
	▼
METABOLIC ACTIVITY (Fluorescence Signal)
```

Applications: Breeding programs, health assessment, stress tolerance screening.

---

## Getting Started

1. Clone the repository and open in RStudio.
2. Explore `index.qmd` for technical documentation.
3. Run analysis scripts in `scripts/` for each project.
4. View outputs in `output/` and figures in `figures/`.

## Contact

For questions, see the [GitHub repo](https://github.com/RobertsLab/resazurin-assay-development) or contact the authors listed in `index.qmd`.

