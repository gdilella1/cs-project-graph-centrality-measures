# Limiting Behavior & Asymptotic Analysis of Parameter-Dependent Network Centrality Measures

This repository contains the MATLAB numerical framework, theoretical report, and empirical validation suite for analyzing the limiting behavior of parameter-dependent graph centrality metrics. The project tests and validates the theoretical results established by **C. Klymko and M. Benzi (SIAM J. Matrix Anal. Appl., 2015)** regarding the convergence of matrix-function-based centrality measures on large-scale real-world networks.

---

## 📌 Executive Summary

Parametric network centrality metrics — such as **Exponential Subgraph/Total Communicability** and **Resolvent/Katz Centralities** — depend on a tuning parameter $\beta \in (0, t^*)$ that controls the penalty applied to longer walks within a graph. While these measures depend continuously on $\beta$, the resulting node rankings exhibit non-continuous structural shifts.

This project rigorously investigates and empirically proves the asymptotic convergence of these measures at the spectrum boundaries:
1. **As $\beta \to 0^+$:** The rankings converge to simple local metrics (**Out-Degree / In-Degree Centrality**).
2. **As $\beta \to t^*$ (Spectral Radius Boundary):** The rankings converge to dominant global eigenvector measures (**Right / Left Eigenvector Centrality**) governed by Perron-Frobenius theory.

---

## 🧮 Mathematical Framework & Measures Evaluated

The codebase implements and evaluates the following broadcast and receive centrality formulations over a graph adjacency matrix $A$:

| Centrality Metric | Mathematical Formulation | Target Parameter Limit ($\beta \to 0^+$) | Target Parameter Limit ($\beta \to t^*$) |
| :--- | :--- | :--- | :--- |
| **Exponential Subgraph ($ESC$)** | $[e^{\beta A}]_{ii}$ | Out-Degree $d^{\text{out}}$ | Right Eigenvector $q_1$ |
| **Exponential Total Communicability ($ETC$)** | $e^{\beta A} \mathbf{1}$ | Out-Degree $d^{\text{out}}$ | Right Eigenvector $q_1$ |
| **Resolvent Subgraph ($RC$)** | $[(I - \beta A)^{-1}]_{ii}$ | Out-Degree $d^{\text{out}}$ | Right Eigenvector $q_1$ |
| **Katz Centrality ($K$)** | $(I - \beta A)^{-1} \mathbf{1}$ | Out-Degree $d^{\text{out}}$ | Right Eigenvector $q_1$ |

> Note: For directed graphs, receive-type metrics are evaluated using the transposed adjacency matrix $A^T$, converging to In-Degree and Left Dominant Eigenvectors respectively.
---

## 💡 Key Algorithmic Innovations & Visualizations

To overcome the inherent non-continuity of node ranking permutations, two custom MATLAB execution suites were engineered:

* **`testAngle.m` (Vector Alignment Analysis):**
  Calculates the inner product angle $\theta(\beta) = \arccos\left(\frac{\mathbf{v}(\beta) \cdot \mathbf{v}\_{\text{limit}}}{\|\mathbf{v}(\beta)\| \|\mathbf{v}\_{\text{limit}}\|}\right)$ between the parameter-scaled score vector and the asymptotic limit vector. Proves vector parallelism as $\theta \to 0$ near spectral boundaries.
* **`testColorful.m` (Chromatically Ordered Ranking Trajectories):**
  Maps node rankings dynamically across the parameter domain using the MATLAB `jet` colormap. Tracks rank swaps, ties, and order stabilization, visually proving the correctness of boundary limit theorems.

---

## 📂 Repository Structure

```text
├── assets/             # Plot outputs from numerical test runs
├── docs/               # Theoretical report (PDF) and reference paper
├── src/                # Core MATLAB numerical solvers and visualization engines
│   ├── testAngle.m     # Angle convergence test code
│   └── testColorful.m  # Dynamic chromatic ranking trajectory code
└── test-graph/         # Real-world network datasets (SNAP Stanford repository)
```

---

## 🔬 Experimental Datasets

All algorithms were benchmarked against real-world network topologies sourced from the **Stanford Large Network Dataset Collection (SNAP)**:
* **Undirected Networks:** Deezer User-User Social Network ($N = 287$).
* **Directed Networks:** EU Research Institution Temporal Email Network ($N = 181$).
* **Large Benchmark Graphs:** MathOverflow network, Autonomous Systems, and Wiki-Vote datasets.

---

## 📖 References

1. **C. Klymko, M. Benzi**, *"On the Limiting Behavior of Parameter-Dependent Network Centrality Measures"*, SIAM Journal on Matrix Analysis and Applications, 36(2), 686–706, 2015.
2. **J. Leskovec, A. Krevl**, *SNAP Datasets: Stanford Large Network Dataset Collection*, 2014.
