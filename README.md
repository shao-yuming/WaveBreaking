# Wave Breaking Simulation 

- This repository provides a **wave breaking simulation framework** built on the open-source solver **[Basilisk](http://basilisk.fr/)**.  
- It is designed to simulate **air–water two-phase breaking waves**, generate visualization-ready outputs, and analyze vortical structures using Liutex.
- **Anyone can use our post-processing code or extend the wave breaking code to other two-phase flow problems**
---

## Features
- Implements **two-phase Navier–Stokes solver** with adaptive mesh refinement (AMR) in Basilisk.
- Supports **parameter modification** for different breaking conditions.
- Integrates **Liutex-Omega vortex identification** for 3D vortical structure analysis.
- Generates **surface VTP files** and **volume VTU files** for visualization in ParaView.
- Outputs droplet and bubble statistics, energy budget, and dissipation rates.

---

## File Descriptions

- **`wave.c`** – Main simulation file.  
  - Controls the setup (fluid properties, gravity, domain, initial conditions).  
  - The `ak` variable defines the **wave steepness** (modify it to change breaking intensity).  
  - Handles adaptive mesh, surface tracking, droplet/bubble counting, and energy/dissipation diagnostics.

- **`liutex_omega.h`** – Implements the **third-generation Liutex-Omega vortex identification method**.

- **`output_surfaces.h`** – Exports the free-surface field (`f`) as **VTP files**, which can be loaded into ParaView.

- **`output_vtu_foreach.h`** – Exports full-field **VTU files** (velocity, pressure, phase fraction, vorticity, Liutex vectors) for ParaView visualization.
- **`adapt_wavelet_leave_interface.h`** – Ensures that the highest mesh resolution is adopted near the free surface.

---

## Usage

1. Install **Basilisk** following instructions at [http://basilisk.fr/](http://basilisk.fr/).
2. Clone this repository:
   ```bash
   git clone https://github.com/shao-yuming/WaveBreaking.git
   cd WaveBreaking
3. Compile and Run ( MPI parallelism is used, and the case is 192 core, which can be changed according to requirements):
   ```bash
   CC99='mpicc -std=c99' qcc -O2 -Wall -D_MPI=1 wave.c -o wave -L$BASILISK/gl -lglutils -lfb_tiny -lm -disable-dimensions
   mkdir surface && mkdir vtu
   mpirun -np 192 ./wave
   
---

## Results
- Free surface files: surface/surface_*.vtp and .pvtp

- Full field VTU files: vtu/all_*.vtu and .pvtu

- Droplet and bubble statistics: droplets.dat and bubbles.dat

- Energy budget: budgetWater.dat and budgetAir.dat
     
---

## Visualization
Load the .pvtu or .pvtp files in ParaView for post-processing and visualization.


