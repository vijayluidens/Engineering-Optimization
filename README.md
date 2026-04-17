# Engineering-Optimization
**F1 Front Wing Optimization Project (Composite Box Beam / Wing Spar)**

This repository contains MATLAB code for a structural optimization study of a thin‑walled **rectangular composite (CFRP) box beam**, motivated by a Formula 1 wing spar (the scripts also reference a rear wing spar in the console output). The design problem optimizes a 4‑variable cross‑section subject to multiple structural constraints.

---

## What the code does

### Optimization problem
- **Design variables:**  
  
  \[  x = [h,\; b,\; t,\; \theta]  \]
  where:
  - `h` = box height (mm)  
  - `b` = box width (mm)  
  - `t` = wall thickness (mm)  
  - `theta` = off‑axis ply angle (deg)

- **Objective:** minimize beam weight (computed by `beam_weight.m`)
- **Constraints:** computed by `beam_constraints.m` (stress / shear / deflection / twist / buckling style constraints; the scripts treat feasibility as `g(x) <= 0`)

---

## Repository layout (key files)

### Main “run everything” entry point
- **`run_all.m`**  
  Master driver script. Runs the workflow in report order:
  1. Initial investigation (`section_investigation.m`)
  2. Design space analysis (`section3_design_space.m`)
  3. Simplified optimization (`section4_simplified.m`)
  4. Full 4‑D optimization (`section5_optimization.m`)  
  Also creates an output folder called `figures/` and saves all open figures as **PNG + PDF**.

### Report / analysis sections
- **`section_investigation.m`**  
  Boundedness, monotonicity, convexity discussion + numerical noise checks + laminate property sweeps (via `compute_effective_properties.m`).
- **`section_monotonicity_test.m`**  
  Additional monotonicity-related testing (supporting script).
- **`section3_design_space.m`**  
  Generates 2‑D slices of the 4‑D design space and plots objective/constraint contours (multiple figures).
- **`section4_simplified.m`**  
  Reduced-dimension optimization:
  - 2‑D optimization in `(b, t)` with fixed `h` and `theta`
  - 1‑D optimization in `t` with penalties  
  Also compares custom 1‑D methods.
- **`section5_optimization.m`**  
  Full 4‑D constrained optimization and method comparisons:
  - `fmincon (SQP)` baseline (multi-start)
  - Augmented Lagrangian–style approach with **Steepest Descent** and **Conjugate Gradient**
  - **SLP with move limits** (multi-start)  
  Produces a comparison table in the console output.

### Beam model + composite property functions
- **`get_parameters.m`**  
  Central place for constants, bounds, material properties, loads, and allowables used throughout the project.
- **`beam_weight.m`**  
  Objective function: computes beam weight.
- **`beam_constraints.m`**  
  Computes constraint vector `g(x)` and (in some calls) returns extra diagnostic info.
- **`compute_effective_properties.m`**  
  Computes effective laminate properties as a function of ply angle `theta`.

### Custom optimization algorithms
- **`golden_section_search.m`**  
  1‑D golden section search (used for line searches and comparisons).
- **`custom_hybrid_search.m`**  
  Custom 1‑D hybrid optimizer:
  - Stage 1: golden section search (bracketing)
  - Stage 2: successive parabolic interpolation (local convergence)
- **`steepest_descent.m`**  
  Self‑implemented N‑D steepest descent with finite-difference gradients + line search.
- **`conjugate_gradient.m`**  
  Self‑implemented N‑D conjugate gradient optimizer.
- **`slp_move_limits.m`**  
  Sequential Linear Programming approach with move limits.

### Penalty / AL objective wrappers
- **`penalty_objective.m`**
- **`al_objective.m`**

### Outputs / figures
- **`Report Figures/`**  
  Pre-generated figures (repository folder).
- **`figures/`** *(generated at runtime)*  
  Created by `run_all.m` and populated with PNG/PDF outputs of figures.

---

## How to run

### Option A (recommended): run the full workflow
1. Open MATLAB in the repository folder.
2. Run:
   ```matlab
   run_all
   ```
3. Check generated plots in the `figures/` directory.

### Option B: run sections individually
You can run any of these (they expect a parameter struct `p` from `get_parameters.m`):
```matlab
p = get_parameters();
section_investigation(p);
section3_design_space(p);
section4_simplified(p);
section5_optimization(p);
```

---

## Requirements / notes
- MATLAB with Optimization Toolbox is required for `fmincon` usage.
- Many scripts use random multi-starts; seeds are set in places (e.g., `rng(42)`) for repeatability.

---

## License
See `LICENSE`.