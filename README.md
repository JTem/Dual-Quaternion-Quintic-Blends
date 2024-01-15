# Dual Quaternion Quintic Blends (DQQB) Algorithm

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JTem/Dual-Quaternion-Quintic-Blends/main?labpath=interpolator_example.ipynb)

## Overview
This repository contains the implementation of the Dual Quaternion Quintic Blends (DQQB) algorithm, as detailed in our paper: "Dual Quaternion Quintic Blends: $\mathcal{C}^2$-Continuous, Time-Optimized Interpolation of Unit Dual Quaternions". The DQQB algorithm is a novel method for interpolating unit dual quaternions, offering $\mathcal{C}^2$-continuous trajectories, ideal for robotic applications requiring high precision and smooth motion.




## Features
- Implementation of the DQQB algorithm for trajectory planning.
- Examples demonstrating the use of the algorithm in various scenarios.
- Tools for evaluating and visualizing the trajectories.
- Comprehensive documentation for ease of use and customization.

## Examples
Two examples are given, 
the first one shows the execution of the trajectory depicted in the paper:


https://github.com/JTem/Dual-Quaternion-Quintic-Blends/assets/34837129/2ffc94f6-3057-4cf5-b16d-d00d01b11b46



The second video gives a small demonstration of a welding trajectory, which was executed by a LARA8 from Neura Robotics GmbH:


https://github.com/JTem/Dual-Quaternion-Quintic-Blends/assets/34837129/8c4ee954-03ad-400b-9add-ad02c1ef73a8


## Installation
Instructions to set up the environment and install necessary dependencies.

The quaternion and dual quaternion maths can be found in the package [neura_dual_quaternions](https://github.com/JTem/neura_dual_quaternions)

```bash
pip install neura-dual-quaternions
```

For interactive plots, the packages jupyter-matplotlib and jupyter-widgets are needed.

```bash
pip install ipympl
```

```bash
pip install ipywidgets
```

Other site-packages:

```bash
pip install numpy
```

```bash
pip install matplotlib
```

## Contact

For questions, please contact Jens Temminghoff at jens.temminghoff@neura-robotics.com.

## Authors

- Jens Temminghoff
- Marcel Huptych
- Jan Wiartalla
- Markus Schmitz
- Burkhard Corves
- Mathias Hüsing



