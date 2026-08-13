# Aditya-L1 — Final Mission Digital Twin

A self-contained C++20/OpenGL visualization and educational physics model for the Sun–Earth Lagrange-point system, with a mission-oriented Aditya-L1 view.

## Final-revision goals

This last revision specifically addresses the visual problems in the previous build:

- **L1, L2, L3, L4 and L5 are all calculated from the Sun–Earth CR3BP and all five are shown.**
- L1/L2/L3 receive **readability markers and leader lines** because their true coordinates are extremely close to or inside the rendered Sun/Earth discs at overview scale. The leader line points back to the actual calculated location; the physics coordinates are not moved.
- The default camera is a **top-down system overview**, so all five Lagrange points are in view at startup.
- Right-mouse drag switches to a free **3-D camera orbit**; the mouse wheel zooms.
- The spacecraft is shown as **one spacecraft only**, moving slowly around an **illustrative 3-D halo** centered on the calculated L1 position. The period is tied to the ISRO reference value of 177.86 days, but the path is explicitly not the operational ephemeris.
- Sun and Earth use higher-resolution procedural 3-D spheres with animated shaders. The renderer has a more detailed spacecraft body and solar panels.
- Space-weather layers are visibly animated: **solar wind, solar radiation, solar flares/CME arcs, solar magnetic-field loops, and Earth's magnetic field/magnetotail**.
- The UI contains a complete **L1–L5 table** with normalized coordinates and Earth/Sun distances.
- Experimental spacecraft placement supports **L1, L2, L3, L4, L5, or a custom normalized position**.

## Scientific model

### Sun–Earth circular restricted three-body problem (CR3BP)

The normalized rotating-frame equations are:

```text
x¨ - 2 y˙ = ∂Ω/∂x
y¨ + 2 x˙ = ∂Ω/∂y
z¨       = ∂Ω/∂z
```

with

```text
Ω = 1/2 (x² + y²)
    + (1-μ)/r1
    + μ/r2
```

and

```text
μ = M_Earth / (M_Sun + M_Earth)
```

The implementation uses the Sun and Earth masses below and keeps the numerical physics in normalized CR3BP coordinates.

### Numerical integration

The experimental trajectory uses fourth-order Runge–Kutta (RK4). The Jacobi integral is evaluated as a consistency diagnostic:

```text
C = 2Ω - (vx² + vy² + vz²)
```

### Lagrange-point calculation

L1, L2 and L3 are obtained by solving the x-axis equilibrium equation numerically with safeguarded bracketing/bisection. L4 and L5 are the exact triangular CR3BP locations:

```text
L4 = (1/2 - μ, +√3/2, 0)
L5 = (1/2 - μ, -√3/2, 0)
```

The renderer preserves the calculated coordinates. Only the **display marker** for L1/L2/L3 may be offset for readability, with a leader line back to the true point.

## Constants used

```text
G            = 6.67430e-11 m^3 kg^-1 s^-2
M_Sun        = 1.98847e30 kg
M_Earth      = 5.9722e24 kg
1 AU         = 149,597,870.7 km
μ            ≈ 3.00340566519016e-06
```

For the simple Sun–Earth circular reference model:

```text
mean motion n        ≈ 1.99101677e-7 rad/s
velocity scale       ≈ 29.7852 km/s
characteristic time  ≈ 58.13 days
```

## Calculated Sun–Earth Lagrange points

With the constants above, the numerical model gives approximately:

```text
L1  x ≈ 0.99002667658     y = 0
L2  x ≈ 1.01003403275      y = 0
L3  x ≈ -1.00000125142     y = 0
L4  x ≈ 0.49999699659      y = +0.86602540378
L5  x ≈ 0.49999699659      y = -0.86602540378
```

The UI computes the Earth and Sun distances from these normalized coordinates at runtime rather than hard-coding them.

## Aditya-L1 mission reference

Official ISRO mission information used by the project/documentation includes:

- Launch: **2 September 2023** on PSLV-C57.
- Halo-orbit insertion: **6 January 2024**.
- Sun–Earth L1 region: roughly **1.5 million km from Earth** toward the Sun.
- Halo-orbit period: about **177.86 Earth days**.
- Nominal mission life: about **5 years**.
- Seven scientific payloads: VELC, SUIT, SoLEXS, HEL1OS, ASPEX, PAPA and MAG.

ISRO's first halo-orbit completion report states that the spacecraft completed its first halo orbit on 2 July 2024 and that station-keeping manoeuvres were required because perturbing forces cause the real spacecraft to depart from its nominal path.

ISRO's July 2026 Announcement of Opportunity says more than **30 TB** of Aditya-L1 scientific data are in the public domain. ISRO's July 2026 iron-fluorescence report describes observations associated with **47 X-class solar flares** during 2024.

Official sources:

- https://www.isro.gov.in/ISRO_EN/Aditya_L1.html
- https://www.isro.gov.in/ISRO_EN/halo-orbit-insertion-adtya-l1.html
- https://www.isro.gov.in/ISRO_EN/Aditya_L1_Mission_Completion_of_First_Halo_Orbit.html
- https://www.isro.gov.in/AdityaL1_Mission_Announcement_of_opportunity.html
- https://www.isro.gov.in/Iron_Fluorescence_on_Sun_Aditya-L1_Observation.html

## Space-weather visualizations

The renderer includes clearly separated **visualization models** for:

- Solar wind particles/streaks
- Solar electromagnetic radiation rays
- Solar-flare/CME arcs and particles
- Large-scale solar magnetic-field loops
- Earth's magnetosphere and magnetotail

These layers are **not operational space-weather forecasts** and should not be interpreted as real-time field or plasma solutions. They are designed to communicate the physical concepts Aditya-L1 studies.

## Camera and controls

```text
RMB drag        = rotate 3-D camera
Mouse wheel     = zoom
O               = toggle overview camera
R               = reset overview
P               = pause
I               = show/hide information panel
F11             = fullscreen
```

## UI controls

The left information panel provides:

- simulation speed
- pause/follow controls
- L1–L5 visibility
- world labels
- halo visibility
- experimental CR3BP path
- flare/CME visualization
- radiation visualization
- solar wind
- Sun magnetic field
- Earth magnetic field
- flare strength
- wind density
- spacecraft position mode
- custom normalized spacecraft position
- complete L1–L5 coordinate/distance table
- selected-point information
- live spacecraft distances and solar irradiance estimate
- physics constants and payload/missions sections

## Building on Windows

Open **Developer PowerShell for Visual Studio** in the project folder:

```powershell
cmake -S . -B build
cmake --build build --config Debug
.\build\Debug\AdityaL1.exe
```

If `build` already contains an incompatible old configuration, delete it first:

```powershell
Remove-Item -Recurse -Force .\build
cmake -S . -B build
cmake --build build --config Debug
```

## Building on Linux

Install a C++20 compiler, CMake, OpenGL development libraries and the GLFW/GLX/X11 development dependencies for your distribution. Then:

```bash
cmake -S . -B build
cmake --build build -j
./build/AdityaL1
```

## Building on macOS

Install Xcode Command Line Tools and CMake. Then:

```bash
cmake -S . -B build
cmake --build build -j
./build/AdityaL1
```

The CMake file links the platform OpenGL library and uses GLFW/GLM/ImGui fetched by CMake.

## GLAD / GLFW / GLM / Dear ImGui

The project contains the generated GLAD files locally in `external/glad/`. GLFW, GLM and Dear ImGui are fetched by CMake during configuration. A network connection is therefore required on the first configuration unless the corresponding CMake FetchContent sources are already cached.

## Architecture

The final project is intentionally compact:

```text
AdityaL1_Final/
├── CMakeLists.txt
├── README.md
├── external/
│   └── glad/
├── assets/
│   └── models/
└── src/
    └── main.cpp
```

The monolithic `main.cpp` keeps the final demo easy to copy, compile and present under time pressure while still separating the physics helpers, numerical integrator, rendering helpers and UI sections conceptually.


## Time-dependent Sun/Earth and L1-L5 coordinates

The visualization now uses the **circular Sun-Earth CR3BP consistently in both frames**:

- In the rotating CR3BP frame, L1-L5 are fixed at the numerically calculated coordinates.
- In the inertial display, the complete Sun-Earth configuration rotates with the CR3BP mean motion \(n\).
- The Sun and Earth therefore move continuously as `simulation day` advances.
- L1, L2, L3, L4 and L5 rotate with the same Sun-Earth angular motion, so their displayed positions remain mathematically consistent with the calculated rotating-frame coordinates.
- The mission spacecraft's illustrative halo is generated in the rotating frame around the exact L1 coordinate, then transformed into the inertial display frame. Its position therefore changes continuously with both the 177.86-day halo phase and the Sun-Earth orbital phase.

This is **exact for the circular restricted three-body model**. It is not a JPL/Horizons ephemeris and does not model the real eccentric Earth orbit or the operational Aditya-L1 station-keeping solution. A real mission ephemeris would require mission state history and a higher-fidelity force model.

The UI shows both the fixed rotating-frame L-point coordinates and the current inertial normalized coordinates so that the changing positions can be distinguished clearly.

## Scientific limitation

The moving halo envelope in the graphics is **illustrative**. It uses the official 177.86-day reference period and is centered on the correctly calculated Sun–Earth L1 point, but it is not a reconstructed flight ephemeris. Reproducing the operational Aditya-L1 trajectory requires mission state history plus higher-fidelity perturbation and station-keeping modelling.

The same distinction applies to the space-weather visualization layers.

## Presentation recommendation

For a competition/demo:

1. Start in **Overview Camera** so all L1–L5 points are visible.
2. Keep **Flares, Radiation, Solar Wind, Sun B and Earth B** enabled.
3. Enable **Follow spacecraft** for the dynamic mission view.
4. Select **L1** in the point table to show its exact calculated coordinates and distances.
5. Switch off Overview and use RMB drag to demonstrate the 3-D structure of the halo and fields.
