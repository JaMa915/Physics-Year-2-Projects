# Physics Year 2 Projects – University of Bath

Two numerical modelling projects from my MPhys Physics degree at the University of Bath.

## 1. Reflectance of Dielectric Multilayers (C)

Folder: [`dielectric-coating/`](dielectric-coating/)

C programs that use the **transfer matrix method** to calculate the reflectance spectrum (300–800 nm) of optical dielectric coatings on glass:

| File | What it models |
|------|----------------|
| `single_layer.c` | A single low-index dielectric layer |
| `bilayer.c` | One low-index / high-index bilayer |
| `multilayer.c` | A quarter-wave stack of *N* bilayers (N entered by the user) |

**Key results:** a single layer reflects only ~1.4–1.7% of visible light, while stacking quarter-wave bilayers (1 to 24 pairs) produces a high-reflectance stop band, the behaviour of a dielectric Bragg mirror.

The programs use `clib.h`, the complex-number and matrix library provided by the university, which is not included here.

Full write-up: [`Reflectance_of_Dielectric_Multilayers_report.pdf`](dielectric-coating/Reflectance_of_Dielectric_Multilayers_report.pdf)

## 2. Predator–Prey and Zombie Outbreak Population Dynamics (Python)

Folder: [`predator-prey/`](predator-prey/)

A Jupyter notebook that models population dynamics with the **Lotka–Volterra** equations (coupled nonlinear ODEs), using NumPy, SciPy and Matplotlib.

- **Rabbits and foxes:** compares Euler and RK4 numerical solutions. Euler trajectories spiral away from the true cycle, while RK4 preserves it, giving an oscillation period of 11.49 days.
- **Zombie outbreak:** varies the zombification rate from 0.002 to 0.02 and fits the outbreak time, showing it is inversely proportional to the rate.
- **Cure model:** adds a time-dependent cure term, under which the human population recovers and the zombies die out.

Notebook: [`predator_prey_zombie.ipynb`](predator-prey/predator_prey_zombie.ipynb) · Full write-up: [`Predator_Prey_Zombie_report.pdf`](predator-prey/Predator_Prey_Zombie_report.pdf)
