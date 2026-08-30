# System Identification of Admittance Models for Large Real-World Objects

This repo identifies dynamic/admittance models for large real-world objects — a wheelbarrow, a door, and their handles — from VICON motion-capture, load-cell, and IMU data. Each object has its own `build_data.m` (raw data → processed `data` struct) and `run_regression.m` (loads that struct, fits the model via `lsqnonlin`, and plots validation results). Shared data-loading and geometry utilities live in `utils/` and are used by all three pipelines.

Developed and tested in **MATLAB R2021a**.

See `supplement.pdf` for detailed derivation notes covering the models fit by these pipelines.

## Running the wheelbarrow pipeline

From MATLAB, with the repo root as the working directory:

```matlab
run('wheelbarrow/build_data.m')      % processes raw data, saves wheelbarrow/results/wheelbarrow_data.mat
run('wheelbarrow/run_regression.m')  % loads that .mat, fits the model, plots results
```

`build_data.m` must be run first (or whenever the raw data or processing changes) — `run_regression.m` just loads its saved output.

## Running the door pipeline

From MATLAB, with the repo root as the working directory:

```matlab
run('door/build_data.m')      % processes raw data, saves door/results/door_data.mat
run('door/run_regression.m')  % loads that .mat, fits the model, plots results
```

As with the wheelbarrow pipeline, `build_data.m` must be run first — `run_regression.m` just loads its saved output.

The door model identifies hinge inertia, viscous/Coulomb friction, and tanh-blended backcheck/sweep/latch damper zones, with a separately fit quasistatic cubic spring (closer) model subtracted out of the regression reference.

## Running the handle pipelines

From MATLAB, with the repo root as the working directory:

```matlab
run('handles/build_data.m')      % processes both handles, saves handles/results/wb_handle_data.mat and handles/results/door_handle_data.mat
run('handles/run_regression.m')  % loads both, fits each via log-Cholesky lsqnonlin, plots results
```

Both scripts process the wheelbarrow handle and the door handle in a single run — there's no per-object toggle. `build_data.m` loads and tares the raw load-cell/IMU data per handle but does no train/test split; `run_regression.m` owns the train/test index selection and builds the stacked Newton-Euler regressor itself before fitting a 10-parameter log-Cholesky inertial model plus a 6-component bias wrench, regularized toward a CAD cylinder prior for each handle.
