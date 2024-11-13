# Session 5
# Ruize Li, rl737, Jesus College, Cambridge
# 2024/11/12
# Planet Simulation

'''
This is a programme simulating the motion of a planet in a central field of any order. 
It can be used to check the fact that only central fields of order 1 and -2 have closed orbits,
which is proved by the Lagrange-Kepler theorem.
The programme uses the Euler-Cromer method & the Leapfrog method & the Runge-Kutta method to solve the motion equations. 
The programme also plots the trajectory, velocity, acceleration, and field of the planet.
'''

# import libraries
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# define a class for the field
class CentralField:
    # the field is defined by its centre, charge, field constant, and field order
    def __init__(self, charge, field_constant, field_order) -> None:
        # the centre is set to the origin by default
        self.centre = np.array([0,0])
        # in gravitational fields, the charge is just the mass
        self.charge = charge
        self.field_constant = field_constant
        # the field is proportional to r ** field_order
        self.field_order = field_order

    def set_centre(self, centre) -> None:
        # if you want to change the centre to another position, use this function
        self.centre = centre

    def field_at(self, position) -> np.ndarray:
        # returns the field as a vector at a given position
        distance = np.linalg.norm(position - self.centre)
        return (-1) * self.field_constant * self.charge * distance**(self.field_order-1) * (position - self.centre)
    
    def __str__(self) -> str:
        return f'Central Field with Charge = {self.charge}, Field Constant = {self.field_constant}, Field Order = {self.field_order}'
    
# define a class for the motion of the planet
class PlanetSimulation:
    def __init__(self, central_field, planet_position, planet_velocity, time_step, time_duration) -> None:
        # the sign of the planet's charge determines whether the field is attractive or repulsive
        self.sign = 'minus' # by default, the planet is attractive
        self.central_field = central_field
        self.planet_position = planet_position # the initial position
        self.planet_velocity = planet_velocity # the initial velocity
        self.time_step = time_step # the time step
        self.time_duration = time_duration # the time duration we consider

    def set_sign(self, sign) -> None:
        # if you want to change the sign of the planet's charge, use this function
        self.sign = sign

    def acceleration_at(self, position) -> np.ndarray:
        # returns the acceleration as a vector at a given position
        if self.sign == 'plus' and self.central_field.charge > 0 or self.sign =='minus' and self.central_field.charge < 0:
            return -self.central_field.field_at(position)
        elif self.sign =='minus':
            return self.central_field.field_at(position)
        else:
            raise ValueError('Sign must be either "plus" or "minus".')
    
    def simulate(self) -> tuple:
        # the main function, simulates the motion of the planet and returns the time, position, and velocity
        time = np.arange(0, self.time_duration, self.time_step)
        position = np.zeros((len(time), 2))
        velocity = np.zeros((len(time), 2))
        position[0] = self.planet_position
        velocity[0] = self.planet_velocity

        '''
        # use the Euler-Cromer method to solve the motion equations
        for i in tqdm(range(1, len(time))):
            velocity[i] = velocity[i-1] + self.acceleration_at(position[i-1]) * self.time_step
            position[i] = position[i-1] + velocity[i] * self.time_step
        '''
            
        '''
        # use the Leapfrog method to solve the motion equations
        for i in tqdm(range(1, len(time))):
            velocity[i] = velocity[i-1] + self.acceleration_at(position[i-1]) * self.time_step/2
            position[i] = position[i-1] + velocity[i] * self.time_step/2
            velocity[i] = velocity[i] + self.acceleration_at(position[i]) * self.time_step/2
            position[i] = position[i] + velocity[i] * self.time_step/2
        '''

        #'''
        # use the Runge-Kutta method to solve the motion equations
        for i in tqdm(range(1, len(time))):
            k1 = self.acceleration_at(position[i-1])
            k2 = self.acceleration_at(position[i-1] + k1 * self.time_step/2)
            k3 = self.acceleration_at(position[i-1] + k2 * self.time_step/2)
            k4 = self.acceleration_at(position[i-1] + k3 * self.time_step)
            velocity[i] = velocity[i-1] + (k1 + 2*k2 + 2*k3 + k4) * self.time_step/6
            position[i] = position[i-1] + velocity[i] * self.time_step
        #'''

        return time, position, velocity
    
    def plot_trajectory(self) -> None:
        # the function plots the trajectory of the planet
        time, position, velocity = self.simulate()
        plt.figure(figsize=(8,8))
        plt.plot(position[:,0], position[:,1])
        plt.axis('equal')
        plt.grid(True)
        plt.xlabel(r'$x$ [m]')
        plt.ylabel(r'$y$ [m]')
        plt.title(f'Trajectory of Planet in Central Field of order {self.central_field.field_order}')
        plt.show()

    def plot_velocity(self) -> None:
        # the function plots the velocity of the planet
        time, position, velocity = self.simulate()
        plt.figure(figsize=(8,6))
        plt.plot(time, np.linalg.norm(velocity, axis=1))
        plt.grid(True)
        plt.xlabel(r'Time [s]')
        plt.ylabel(r'Velocity [${\rm ms^{-1}}$]')
        plt.title(f'Velocity of Planet in Central Field of order {self.central_field.field_order}')
        plt.show()

    def plot_acceleration(self) -> None:
        # the function plots the magnitude of the acceleration of the planet
        time, position, velocity = self.simulate()
        acceleration = np.zeros(len(time))
        for i in tqdm(range(1, len(time))):
            acceleration[i] = np.linalg.norm(self.acceleration_at(position[i]))
        plt.figure(figsize=(8,6))
        plt.plot(time, acceleration)
        plt.grid(True)
        plt.xlabel(r'Time [s]')
        plt.ylabel(r'Acceleration [${\rm ms^{-2}}$]')
        plt.title(f'Magnetude of Acceleration of Planet in Central Field of order {self.central_field.field_order}')
        plt.show()

    def __str__(self) -> str:
        return f'Planet Simulation with Central Field of order {self.central_field.field_order}.'

# an example is Earth's motion around the Sun
# Let's assume that the gravitational field is proportional to r**(-2+0.004)
# i.e. not exactly Newtonian, and see what happens. 
field = CentralField(1.989e30, 6.67430e-11, -2 + 0.04)
motion = PlanetSimulation(field, np.array([148.5e9, 0]), np.array([0, 30e3]), 3600, 15 * 365 * 24 * 3600)

# Make the plots
motion.plot_acceleration()
motion.plot_velocity()
motion.plot_trajectory()