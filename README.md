# Stochastic Modeling of Bacteriophage Infection Cycles

This project provides a MATLAB-based implementation for simulating bacteriophage infection cycles, focusing on stochastic processes that distinguish between the **lytic** and **lysogenic** cycles. The model leverages **Gillespie’s Algorithm** to capture the random nature of chemical reactions occurring at low molecular counts, providing a more realistic simulation than deterministic approaches. This project was completed as part of a course I took at Stony Brook University, BIO/AMS332 with professor David Green.

## Table of Contents
1. [Project Overview](#project-overview)
2. [Key Concepts](#key-concepts)
   - [Lytic and Lysogenic Cycles](#lytic-and-lysogenic-cycles)
   - [Mutual Exclusivity and Stable Points](#mutual-exclusivity-and-stable-points)
   - [Stochastic vs. Deterministic Reactions](#stochastic-vs-deterministic-reactions)
3. [Modeling the Process](#modeling-the-process)
   - [Chemical Master Equation](#chemical-master-equation)
   - [Gillespie’s Algorithm](#gillespies-algorithm)
4. [Simulation Details](#simulation-details)
5. [Results and Observations](#results-and-observations)

## Project Overview
In this project, I developed a **stochastic model of bacteriophage infection cycles** to simulate how bacteriophages behave under different molecular conditions within a bacterial host. The model specifically looks at how bacteriophages randomly enter either the **lytic** or **lysogenic** cycle, based on molecular interactions in a low-copy environment.

## Key Concepts

### Lytic and Lysogenic Cycles
1. **Lytic Cycle**: The bacteriophage uses the host cell’s machinery to replicate its DNA. This aggressive replication can lead to cell death and, potentially, the eradication of the host bacteria population.
2. **Lysogenic Cycle**: The viral DNA integrates into the host's genome and replicates passively along with the host's DNA, without killing the cell. This allows the virus to propagate across generations until stress activates the lytic cycle.

### Mutual Exclusivity and Stable Points
The model assumes **mutual exclusivity** between the cycles. Two key transcription factors (`cI` for lysogeny and `cro` for lysis) inhibit each other’s expression, establishing two stable states:
- **Lysis State**: High `cro` and low `cI`.
- **Lysogeny State**: High `cI` and low `cro`.

### Stochastic vs. Deterministic Reactions
The model contrasts **stochastic processes** with **mass-action kinetics**. While deterministic approaches predict average behavior in large populations, this model simulates each reaction event randomly, accounting for low molecular copy numbers and their effects on reaction rates.

## Modeling the Process

### Chemical Master Equation
The **Chemical Master Equation (CME)** calculates the rate of change of probability distributions for reaction occurrences. It considers each molecular species’ state transitions, adapting to varying reaction propensities in the bacterial cell.

### Gillespie’s Algorithm
Gillespie’s Algorithm simulates individual reactions by:
1. Generating two random numbers.
2. Using these to determine **when** and **which** reaction occurs next.
3. Updating population vectors and reaction propensities accordingly.

This approach enables precise simulations at low copy numbers, where molecular collisions occur sporadically.

## Simulation Details
- The model defines the **population vector** as `x = {xp, xr}`, where `xp` is the count of protein molecules and `xr` is the count of RNA molecules.
- **Reactions** include protein synthesis, degradation, RNA synthesis, and RNA degradation, each with a unique propensity calculated based on population levels and reaction constants.
- **Reaction Propensities**:
  - Protein synthesis: `α1(x) = ω * xr`
  - Protein degradation: `α2(x) = χ_prot * xp`
  - RNA synthesis: `α3(x) = μ( xp^2 / (K1/2^2 + xp^2) )`
  - RNA degradation: `α4(x) = χ_RNA * xr`

## Results and Observations
Starting the simulation with no RNA or protein represents initial viral infection. This zero-state emphasizes **stochastic selection** of the cycle the virus will enter, as random chance dictates which gene (cI or cro) is expressed first. Typically:
- In stochastic models, lysogeny is favored due to faster synthesis and slower degradation of `cI`.
- Occasionally, stochastic fluctuations may lead to lysis, particularly if `cro` is synthesized faster than `cI` initially.

## Conclusion
This model illustrates the impact of **random molecular events** on viral life cycle choices, providing a robust alternative to deterministic models. It highlights how small molecular counts in cellular environments create significant variability in outcomes, a phenomenon essential to understanding viral behavior in real-world conditions.

---

> **Note**: This project was completed individually, with reference to course materials and select online resources.
