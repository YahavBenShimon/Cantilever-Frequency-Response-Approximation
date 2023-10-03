# Cantilever Frequency Response Approximation

This repository contains a MATLAB code that approximates the frequency response of a HBN cantilever actuated via Kelvin polarization force.
The approximation utilizes the Galerkin method. The differential equations are transformed into a non-dimensional form to simplify computations.
The program then plots the frequency response and creates a video of the time response at resonance.

## Features

- Computes the natural frequency of the cantilever.
- Sweeps over a range of frequencies to simulate the response.
- Plots the frequency response of the cantilever.
- Generates a video showing the time response of the cantilever at a specific frequency close to its resonance.

## How to Use

1. Clone this repository.
2. Open the MATLAB environment.
3. Navigate to the directory containing the script.
4. Run the script.
5. The frequency response will be saved as 'Frq.png' and the video of the time response at resonance will be saved as 'my_video.mp4'.

## Constants & Parameters

- `epsilon_null`: Permittivity of free space.
- Geometrical properties such as width, length, height, and cross-sectional area.
- Actuation properties including AC and DC voltages.
- Material properties including relative permittivity, mass, density, bending stiffness, spring constant, quality factor, and damping ratio.

## Functions

- `phi`: Represents the mode shape function in a normalized form.
- `U_initial` and `U_dot_initial`: Defines the initial conditions for the system.
- `plotResponseAtResonance`: Plots the Time response of the cantilever at a specific frequency close to resonance and creates a video of the displacement over time.

