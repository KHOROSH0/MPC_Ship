# MPC_Ship

## Overview
This project about application of the Model Predictive Control (MPC) using the example a linear ship. You can watch task and math model in task/Model_of_ship.pdf

## Prerequisites
Make sure you have Matlab installed

## Main Idea
Consider a discrete-time linear system described by

$$
x_{k+1} = Ax_k + Bu_k,
$$

where:

$x_k$ is the system state vector;

$u_k$ is the control input vector;

$A$ is the state transition matrix;

$B$ is the control input matrix.

The objective of MPC is to minimize a finite-horizon cost function:

$$
J =
\sum_{j=0}^{N-1}
\left[
(x_{k+j}-x_{k+j}^{\mathrm{ref}})^T
Q
(x_{k+j}-x_{k+j}^{\mathrm{ref}})
+
u_{k+j}^T R u_{k+j}
\right]
+
(x_{k+N}-x_{k+N}^{\mathrm{ref}})^T
P
(x_{k+N}-x_{k+N}^{\mathrm{ref}}),
$$

subject to the system dynamics and any imposed state and control constraints.

Here:

$N$ is the prediction horizon;

$x_k^{\mathrm{ref}}$ is the reference state;

$Q$ penalizes deviations from the desired state;

$R$ penalizes control effort;

$P$ is the terminal cost weighting matrix.

## 1. src/leinear_ship.m — MPC vs. LQR

This script compares Model Predictive Control with Linear Quadratic Regulator (LQR) control for a linear ship model.

The experiment investigates the relationship between the prediction horizon and the resulting control behavior. For a suitable unconstrained linear-quadratic formulation, as the MPC prediction horizon tends to infinity, the MPC feedback law approaches the infinite-horizon LQR controller, provided the standard stabilizability, detectability, and terminal-cost assumptions are satisfied.

## 2. src/leinear_ship2.m — Ship Heading Maneuver

This script simulates a ship turning through a heading angle of 10 degrees, subject to constraints based on the ship model described in task/Model_of_ship.pdf.
