# VoltControl — MPC-based Voltage Control for Renewable Energy Integration

Research project in partnership with RTE (French Transmission System Operator).

## Context

Renewable energy producers are sensitive to voltage variations and may automatically
disconnect from the grid. This project implements a secondary voltage control system
using Model Predictive Control (MPC) to regulate reactive power and keep load bus
voltages on reference, while preventing generators from hitting their reactive power limits.

## Controllers

Three MPC-based controllers are implemented, sharing the same network and
capacitor/reactor logic but differing in how OLTC transformers are handled:

- **Main controller** — OLTC co-optimized with generator voltages in a single QP
- **Alternative OLTC logic** — OLTC handled by a separate enumeration logic
- **Naive controller** — OLTC driven directly by the voltage reference (baseline)

## Tech Stack

- MATLAB 2025b
- MATPOWER

## Report

See [`rapport.pdf`](Rapport/rapport.pdf) for the full project report.
