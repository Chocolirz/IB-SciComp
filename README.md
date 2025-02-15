# IB Scientific Computing

This repository contains the materials for the IB Computing course and some other computing-related problems during the 2024~2025 academic year at the University of Cambridge. 

Nothing hard, just some notes and exercises. Enjoy!

## Scientific Computing
+ [Session 1: Simple adder](https://github.com/Chocolirz/IB-SciComp/blob/main/SimpleAdder_plus_PolynomialSolver.py)
+ [Session 2: Oscillation](https://github.com/Chocolirz/IB-SciComp/blob/main/Oscillator_Visualiser.py)
+ [Session 3: Monte Carlo (MC)](https://github.com/Chocolirz/IB-SciComp/blob/main/MonteCarlo.py)
+ [Session 4: Method of Bisection](https://github.com/Chocolirz/IB-SciComp/blob/main/Method_of_Bisection.py)
+ [Session 5: Planets](https://github.com/Chocolirz/IB-SciComp/blob/main/Planet_Simulation.py)
+ [Session 6: Collisions](https://github.com/Chocolirz/IB-SciComp/blob/main/Collisions.py)

## Experimental Methods
+ [Q1: Fast Fourier Transform](https://github.com/Chocolirz/IB-SciComp/blob/main/Experimental_Methods_Problem_1.ipynb)
+ [Q17: Linear fit](https://github.com/Chocolirz/IB-SciComp/blob/main/Experimental_Methods_Problem_17.ipynb)
+ [Q18: The $\chi^2$ test](https://github.com/Chocolirz/IB-SciComp/blob/main/Experimental_Methods_Problem_18.ipynb)

## Oscillations, Waves, and Optics
+ [Q7: Simple Pendulum at large angles](https://github.com/Chocolirz/IB-SciComp/blob/main/OWO_Q7.ipynb)
+ [Q15: Reflection and transmission coefficients](https://github.com/Chocolirz/IB-SciComp/blob/main/OWO_Q15.ipynb)
+ [Q19: Propagation of Gausian Wave Packets](https://github.com/Chocolirz/IB-SciComp/blob/main/OWO_Q19.ipynb)
+ [Q31: Diffraction](https://github.com/Chocolirz/IB-SciComp/blob/main/OWO_Q31.ipynb)

## Physics Practicals
+ Experiment 6: Hysteresis used [this project](https://github.com/Chocolirz/IB-SciComp/tree/main/picoscope) to collect data from a picoscope (2206A).
+ Experiment 5: Faraday Rotation used [this project](https://people.phy.cam.ac.uk/db106/pub/) in the PSD folder during dythering. The project was almost flawless, so there's no need to apply any modifications.
+ Experiment 4: Funky pendulum. I didn't have time to look into this, so I just used [Zhen's original version](https://github.com/Zzzzhen1/Funky_Pendulum). 

## Additional
+ Yiran Ma (Peterhouse College) suggested that it's cool to do visualisation of collisions that can estimate $\pi$. I've done the basic calculations [here](https://github.com/Chocolirz/IB-SciComp/blob/main/pi_using_collisions.py). The topic is indeed interesting, but I don't think it's a good way to estimate $\pi$ (see the short notes at the end of the script). It's always better to use MC. 
+ Inspired by my supervisor Matthew Smith, I've done a [rough estimate](https://github.com/Chocolirz/IB-SciComp/blob/main/QM_Q1.ipynb) of the temperature of Millikan's lab in June 1917. This is about the first question in the quantum physics sheet. 
+ I've written two projects that might be insightful for future IB students choosing the topic of their posters. The first one is about [Dzhanibekov effect](https://github.com/Chocolirz/IB-SciComp/blob/main/Dzhanibekov/Dzhani_simulation.ipynb), which describes the free rotation about the principle axis with intermediate inertia. The second one is about [Fabry-Perot etalon](https://github.com/Chocolirz/IB-SciComp/blob/main/Fabry_Perot/FP_simulation.ipynb), which could give a finer spectrum than other interferometers. Note that they both have a variety of applications in different fields, for example, the code for Dzhanibekov effect can also be used to simulate free prosession of symmetric top, and Fabry-Perot etalon happens to support perfect transmission in quantum physics (from a purely mathematical point of view). 