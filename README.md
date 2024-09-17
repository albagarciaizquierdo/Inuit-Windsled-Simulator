This project provides a MATLAB-based simulation tool developed as part of a thesis to model the dynamics and analyze the stability of the Inuit Windsled. The tool allows for the calculation of equilibrium positions, stability analysis, and parametric studies of the Windsled system.

**Features**
- Dynamic Model: Simulates the movement of the sled and kite using two rigid bodies connected by massless tethers.
- Equilibrium Analysis: Calculates the equilibrium positions of the system under specific wind conditions.
- Stability Analysis: Implements both linear stability theory and numerical integration to assess the stability of the equilibrium points.
- Parametric Study: Allows users to analyze how key parameters (wind speed, attachment points, kite and sled properties) influence the system's stability and equilibrium.

**Installation**
1. Clone this repository
2. Ensure you have MATLAB installed
3. Open the project in MATLAB and run the main script to get started

**Usage**

This simulation tool consists of three main scripts, each focusing on a different aspect of the Inuit Windsled's dynamics and stability:
- main.m: Calculates the equilibrium position of the system, selects the most appropriate ODE solver, and verifies the equilibrium through linear stability theory.
- main_ode.m: Runs numerical simulations to study the system’s stability by analyzing how the equilibrium responds to various perturbation scenarios.
- main_par.m: Conducts a parametric study, exploring how factors such as wind speed, attachment points, and kite/sled properties affect equilibrium and stability.
Before running any simulations, you can adjust the system parameters in the parameters.m file. This file contains all relevant settings, such as wind conditions and sled/kite characteristics, which can be modified to suit different test cases.
