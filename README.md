# Figure 2: Mach-5 collision-production validation

This archive contains only the numerical data and source code needed for the
Mach-5 validation plot (Figure 2). It contains no manuscript or Supplement.

## Contents

- `data/figure2_mach5_collision_production.csv`: the ten plotted estimates,
  their horizontal/vertical standard errors, shock-relative positions, and
  paired-difference diagnostics.
- `src/order8_m5_validation.cpp`: standalone hard-sphere DSMC and weak-form
  bracket evaluator. It compares the single-generator theorem prediction with
  an independent angular sampling of the conventional collision bracket for
  the two total-degree-eight Burnett probes `(j,l)=(2,4)` and `(3,2)`.
- `scripts/aggregate_blocks.py`: combines independent block outputs and
  calculates means, standard errors, paired differences, and paired z-scores.
- `scripts/plot_figure2.py`: redraws Figure 2 and prints its two global
  diagnostics.
- `scripts/run_mach5.sh`: compiles the C++ source, runs four independent
  Mach-5 realizations, aggregates them, and plots the regenerated result.
- `runs/`: the four raw 20-block CSV files and their run metadata.

The CSV retains the historical column name `tree`; in the paper this is named
the **single-generator theorem prediction**. The `direct` column is the
independently sampled hard-sphere weak-form collision bracket.

Both reported columns use the manuscript/PRL outgoing-direction normalization
`eta_HS = sqrt(pi)/2`. The DSMC collision dynamics retains the physical
diameter-one cross section `pi`; only the reported collision-production
estimators use `eta_HS`. Consequently this correction rescales both plotted
axes together and leaves the slope, normalized RMS mismatch, and paired
z-scores unchanged.

The labels `-6`, `-3`, `0`, `+3`, and `+6` are cell offsets from the
instantaneous density midpoint of the shock. Negative offsets are upstream
(cold/pre-shock); positive offsets are downstream (hot/post-shock).

## Reproduce the published plot from archived data

Install Python 3 and Matplotlib, then run:

```bash
python3 -m pip install -r requirements.txt
make figure
```

The output is `output/figure2_mach5.png`. The script reports

```text
slope_through_origin=1.001792740
normalized_rms_mismatch_percent=0.179595
```

which round in the figure to `1.001793` and `0.180%`.

## Regenerate collision-bracket estimates

A C++17 compiler is required. Run:

```bash
make validate
```

The simulation uses Mach 5, 50 cells, a target of 18,000 particles,
`dt=3e-4`, 1.5 time units of burn-in, 1.5 time units of sampling, 300 sampled
velocity pairs per station and stored sample, 12 angular samples per pair,
and 20 blocks per realization. Four seeds give 80 block estimates for every
plotted point. Monte Carlo regeneration will not be bitwise identical to the
archived CSV unless the same compiler and standard-library random-number
implementation are used.

Within `order8_m5_validation.cpp`, `active_24` and `active_32` are the explicit
single-generator theorem evaluations of the two order-eight collisional
increments. The angular loop evaluates the same brackets independently from
post-collision velocities. This is the code path used for the comparison; no
paper source files are needed.
