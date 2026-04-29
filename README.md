# VoltControl — MPC-based Voltage Control for Renewable Energy Integration

Research project in partnership with RTE (French Transmission System Operator).

## Context

Renewable energy producers are sensitive to voltage variations and may automatically
disconnect from the grid. This project implements a secondary voltage control system
using Model Predictive Control (MPC) to regulate reactive power and keep load bus
voltages on reference, while preventing generators from hitting their reactive power limits.

## Repository Structure

- **`Main controller/`** — Primary controller: OLTC co-optimized with generator voltages in a single QP. Uses an augmented state formulation to handle the two-step tap delay.
- **`Alternative OLTC logic/`** — OLTC handled by a separate enumeration logic (`action_oltc.m`), independent from the MPC. Tap delay managed by explicit state correction.
- **`Naive controller/`** — Baseline controller: OLTC driven directly by the voltage reference, no cost function, no delay anticipation.
- **`Sensitivity test/`** — Standalone validation scripts comparing analytical sensitivity matrices against finite-difference approximations.

All three controllers share the same network definition (`case9xx_Bsh.m`), capacitor/reactor switching logic (`add_cap.m`), and sensitivity functions.

## Running a Simulation

Open MATLAB, navigate to the desired controller folder, and run `main.m`.

Simulation parameters (horizon, weights, noise, constraints, Jacobian update frequency, etc.) can be adjusted at the top of `main.m`.

## Tech Stack

- MATLAB 2025b
- MATPOWER

## Report

See [`rapport.pdf`](Rapport/rapport.pdf) for the full project report.
