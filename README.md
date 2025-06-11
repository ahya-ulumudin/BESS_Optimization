# BESS Optimization Framework using Hybrid PSO–TS

This repository contains a MATLAB-based optimization framework developed as part of a Master's research project titled:

**"Hybrid PSO–TS Optimization for Optimal BESS Placement, Sizing, and Dispatch in PV-Rich Distribution Networks"**

The framework implements a two-stage hybrid metaheuristic approach combining **Particle Swarm Optimization (PSO)** and **Tabu Search (TS)** to improve the planning and operation of **Battery Energy Storage Systems (BESS)** in radial distribution networks under high photovoltaic (PV) penetration.

## 🌟 Key Features
- Supports both **IEEE 33-bus** and **69-bus** test systems.
- **Two-stage optimization**:
  - Stage 1: Optimal **placement** of BESS units.
  - Stage 2: Optimal **sizing** and **24-hour dispatch scheduling**.
- **Hybrid PSO–TS algorithm** designed for better global exploration and local refinement.
- Incorporates realistic **load and PV generation profiles**.
- Multi-objective fitness evaluation:
  - Voltage deviation minimization
  - Loss reduction
  - Reverse power flow mitigation
  - SOC balance and ramp rate smoothing
- Modular MATLAB scripts for extensibility and testing

## 📁 File Structure

| File/Folder | Description |
|-------------|-------------|
| `MAIN_PLACEMENT.m` | Main script for BESS placement optimization |
| `MAIN_SIZING.m` | Main script for BESS sizing and dispatch optimization |
| `Placement_Optimization_PSO.m` / `TS` / `PSO_TS` | Optimization algorithms for placement |
| `Sizing_Optimization_PSO.m` / `TS` / `PSO_TS` | Optimization algorithms for sizing and dispatch |
| `Placement_Objective.m` / `Sizing_Objective.m` | Multi-objective fitness functions |
| `HourlyLoadFlow.m` | Load flow engine using backward/forward sweep |
| `PV_out_profile.m` / `Residential_Load_Profile.m` | Customizable PV and load data profiles |
| `results/`, `figures_summary/`, `results_summary/` | Output folders for visualizations and analysis |
| `single_mode_fix.m` | Script to initialize and generate base case ranking |

## ⚙️ How to Use

1. **Base Case Initialization**  
   Run `single_mode_fix.m` to perform 24-hour load flow analysis and generate bus ranking.

2. **Placement Optimization**  
   Run `MAIN_PLACEMENT.m` to determine optimal BESS locations.  
   Manually input the resulting bus indices into the `Bus_Placement` field of `MAIN_SIZING.m`.

3. **Sizing and Dispatch Optimization**  
   Run `MAIN_SIZING.m` to perform 24-hour demand scheduling and optimal sizing based on the selected buses.

4. **Refinement (Optional)**  
   Set `use_previous = true` in `MAIN_SIZING.m` to refine results using a previous solution as the initial guess.

5. **Customization**  
   - Switch between IEEE 33-bus or 69-bus by modifying input files:
     - `linedataXXbus.m`, `loaddataXXbus.m`
   - Adjust PV and load profile:
     - `PV_out_profile.m`, `Residential_Load_Profile.m`
   - Output folder names automatically reflect the system size (e.g. `results_33bus`).

## 📊 Output and Visualization

- Results include:
  - Voltage profiles
  - State of Charge (SOC) trajectories
  - BESS dispatch curves
  - Loss and voltage deviation summaries
- Visualization scripts are included in:
  - `display_results.m`
  - `placement_display_results.m`

## 🛠 Requirements
- MATLAB R2023b or newer
- No special toolboxes required, but **Optimization Toolbox** may enhance performance

## 📖 Citation

If you use this framework in your work, please cite:

> Ahya Ulumudin (2025). Hybrid PSO–TS Optimization Framework for BESS Planning. GitHub Repository.  
> [https://github.com/ahya-ulumudin/BESS_Optimization](https://github.com/ahya-ulumudin/BESS_Optimization)

## 📝 License

This project is licensed under the [MIT License](LICENSE).  
You are free to use, modify, and distribute this code for academic or commercial purposes, provided that proper credit is given to the original author.
