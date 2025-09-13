### Main script for MFP_Project
#--------------------------------------------------------------
# This is a project to simulate the mean free path of a particle
# with a certain velocity, entering a medium (2D gas in a box) 
# in which molecules are performing Brownian motion. 
#--------------------------------------------------------------

# import basic libraries
import math 
import random

import numpy as np
import matplotlib.pyplot as plt

# introduce some new libraries
from dataclasses import dataclass # for structured data storage
from typing import List, Tuple, Optional # for type hinting
from matplotlib import animation # for animation (Optional)

'''
First, we need to think about what we need as a class, and what functions should be included in each class.

We need a class for the molecules in the medium, which includes 
    position, velocity, radius, and a function to update their position and velocity.
We need a class for the particle entering the medium, which includes 
    position, velocity, radius, and a function to update its position.
We need a class for the box, which includes width and height.
We need a class for the simulation, which includes all the parameters, 
    a list of molecules, the particle, and functions to run 
    the simulation, detect collisions, and record free paths.

Nothing else is needed.

Lets start coding.
'''

# important function: calculate distance between two points
def distance(a: np.ndarray, b: np.ndarray) -> float:
    '''
    Calculate Euclidean distance between two points with coordinates a and b.

    Parameters:
    a (np.ndarray): Coordinates of point a.
    b (np.ndarray): Coordinates of point b.

    Returns:
    float: Euclidean distance between points a and b.
    '''
    return np.linalg.norm(a - b)

# define a class for molecules in the medium
@dataclass
class Molecule:
    # This class includes everything we need to describe a molecule in our 2D box. 
    pos: np.ndarray # position, shape (2,)
    vel: np.ndarray # velocity, shape (2,)
    radius: float = 0.01 # radius of the molecule, set default value to 0.01

    def step(self, dt: float, thermal_kick_std: float, box_size: Tuple[float, float]) -> None:
        '''
        Update the position and velocity of the molecule for a time step dt.

        Parameters:
        dt (float): Time step for the update.
        thermal_kick_std (float): Standard deviation of the thermal kick.
        box_size (Tuple[float, float]): Size of the box (width, height).
        '''
        # Apply a random thermal "kick" to the velocity (Brownian motion)
        thermal_kick = np.random.normal(0, thermal_kick_std, size=2)
        self.vel += thermal_kick

        # Update position based on current velocity
        self.pos += self.vel * dt

        '''
        How to treat boundaries?
        
        As mentioned before, we can use periodic boundary conditions, or include reflective walls.
        Here, we will use the former. 
        '''
        # Apply periodic boundary conditions
        self.pos[0] = self.pos[0] % box_size[0]
        self.pos[1] = self.pos[1] % box_size[1]

# Define a class for the particle entering the medium
@dataclass
class Particle:
    pos: np.ndarray # position, shape (2,)
    vel: np.ndarray # velocity, shape (2,)
    radius: float = 0.02 # radius of the particle, set default value to 0.02

    # similarly, define a function for each step. 
    def step(self, dt: float, box_size: Tuple[float, float]) -> None:
        '''
        Update the position of the particle for a time step dt.

        Parameters:
        dt (float): Time step for the update.
        box_size (Tuple[float, float]): Size of the box (width, height).
        '''
        # Update position based on current velocity
        self.pos += self.vel * dt

        # We care very much about this particle, and we need to trace it.
        # Therefore, we cannot use periodic boundary conditions here.
        # Instead, we will use reflective boundary conditions.
        Lx, Ly = box_size
        # Reflective boundaries
        for i in (0, 1):
            if self.pos[i] - self.radius < 0:
                self.pos[i] = self.radius
                self.vel[i] *= -1
            elif (self.pos[i] + self.radius) > (Lx if i == 0 else Ly):
                if i == 0:
                    self.pos[i] = Lx - self.radius
                else:
                    self.pos[i] = Ly - self.radius
                self.vel[i] *= -1
        # Luckily, we don't need to apply brownian motion to this particle,
        # so the reflective conditions are easy to implement.

# The box is also an object, which can be defined as a class.
class Box:
    def __init__(self, width: float = 1.0, height: float = 1.0):
        '''
        Initialize the box with given width and height.

        Parameters:
        width (float): Width of the box.
        height (float): Height of the box.
        '''
        self.width = width
        self.height = height

    def size(self) -> Tuple[float, float]:
        '''
        Get the size of the box.

        Returns:
        Tuple[float, float]: Size of the box (width, height).
        '''
        return (self.width, self.height)
    
# We can also define a class for the entire simulation.
class Simulation:
    def __init__(
        self, 
        box: Box,
        n_molecules: int = 100,
        molecule_radius: float = 0.01,
        molecule_thermal_kick_std: float = 2.0,
        particle_radius: float = 0.02,
        particle_speed: float = 1.0,
        dt: float = 0.001,
        collision_distance: Optional[float] = None, # if None, use sum of radii
        max_steps: int = 10000,
    ):
        '''
        Initialize the simulation with given parameters.

        Parameters:
        box (Box): The box in which the simulation takes place.
        n_molecules (int): Number of molecules in the box.
        molecule_radius (float): Radius of each molecule.
        molecule_thermal_kick_std (float): Standard deviation of the thermal kick for molecules.
        particle_radius (float): Radius of the particle entering the medium.
        particle_speed (float): Speed of the particle entering the medium.
        dt (float): Time step for the simulation.
        collision_distance (Optional[float]): Distance threshold for collision detection. If None, use sum of radii.
        max_steps (int): Maximum number of steps to run the simulation.
        '''
        self.box = box
        self.n_molecules = n_molecules
        self.molecule_radius = molecule_radius
        self.molecule_thermal_kick_std = molecule_thermal_kick_std
        self.particle_radius = particle_radius
        self.particle_speed = particle_speed
        self.dt = dt

        # Set collision criteria
        self.collision_distance = collision_distance if collision_distance is not None else (molecule_radius + particle_radius)

        self.max_steps = max_steps

        # Initialise molecules at random positions with random velocities
        self.molecules: List[Molecule] = []
        for _ in range(n_molecules):
            pos = np.array([
                random.uniform(molecule_radius, box.width - molecule_radius),
                random.uniform(molecule_radius, box.height - molecule_radius)
            ])
            vel = np.random.normal(scale=0.0, size=2) # start with zero mean
            self.molecules.append(Molecule(pos=pos, vel=vel, radius=molecule_radius))

        # Initialise our particle, starting from the left, and moving to the right. 
        start_pos = np.array([particle_radius + 1e-5, box.height / 2.]) # start near the left wall, centered vertically
        particle_vel = np.array([particle_speed, 0.0]) # which will not change in magnitude. 
        self.particle = Particle(pos=start_pos, vel=particle_vel, radius=particle_radius)

        # tracking
        self.collision_count = 0 # count number of collisions
        self.free_paths: List[float] = [] # store lengths of free paths
        self._distance_since_last_collision = 0.0 # distance traveled since last collision
        self._last_collision_pos = self.particle.pos.copy() # position of last collision
        self.steps_taken = 0 # count steps taken

    def detect_collision(self) -> bool:
        '''
        Check if the particle has collided with any molecule.

        Returns:
        bool: True if a collision is detected, False otherwise.
        '''
        for molecule in self.molecules:
            if distance(self.particle.pos, molecule.pos) <= self.collision_distance:
                return True
        return False
    
    def collision_index(self) -> Optional[int]:
        '''
        Check if the particle has collided with any molecule and return the index of the molecule.

        Returns:
        Optional[int]: Index of the molecule if a collision is detected, None otherwise.
        '''
        for i, molecule in enumerate(self.molecules):
            if distance(self.particle.pos, molecule.pos) <= self.collision_distance:
                return i
        return None
    
    def handle_collision(self, molecule_index: int) -> None:
        '''
        Record free path and exchange momentum

        Parameters:
        molecule_index (int): Index of the molecule that collided with the particle.
        '''
        # record free path (distance between last collision and this one)
        dist = distance(self._last_collision_pos, self.particle.pos)
        self.free_paths.append(dist)
        self._last_collision_pos = self.particle.pos.copy()
        self.collision_count += 1

        # simple collision response: exchange momentum approximated by swapping velocities
        # or scatter particle isotropically; here we randomise the particle velocity direction
        # which is more physical due to quantum effects at small scales.
        speed = np.linalg.norm(self.particle.vel)
        theta = random.uniform(0, 2 * math.pi)
        self.particle.vel = speed * np.array([math.cos(theta), math.sin(theta)])

        # give molecule a small kick too
        m = self.molecules[molecule_index] # give kick to the collided molecule only. 
        m.vel += np.random.normal(scale=0.1, size=2)

    def step(self) -> None:
        # Perform a single time step for everything in the system. 
        
        # first, move each molecule
        for m in self.molecules:
            m.step(self.dt, self.molecule_thermal_kick_std * math.sqrt(self.dt), self.box.size())

        # then, move the particle
        prev_pos = self.particle.pos.copy() # record previous position
        self.particle.step(self.dt, self.box.size()) # move particle for one step
        self._distance_since_last_collision += distance(prev_pos, self.particle.pos) # update distance traveled for free path calculation

        # check for collision
        collision_molecule_index = self.collision_index()
        if collision_molecule_index is not None:
            self.handle_collision(collision_molecule_index)
        else:
            pass

        self.steps_taken += 1 # increment step count
    
    def run(self, stop_after_collisions: Optional[int] = None) -> None:
        '''
        Run until max_steps or until stop_after_collisions recorded.
        ''' 
        while self.steps_taken < self.max_steps:
            self.step()
            if stop_after_collisions is not None and self.collision_count >= stop_after_collisions:
                break
            # report progress every 10%
            if self.steps_taken % (self.max_steps // 10) == 0:
                print(f'Steps taken: {self.steps_taken}, Collisions: {self.collision_count}')


    '''
    Now we have finished everything about physics and simulation, 
    the rest is just about analysing data and visulising them. 
    '''

    def mean_free_path(self) -> float:
        '''
        Calculate the mean free path from recorded free paths.
        '''
        if len(self.free_paths) == 0:
            return float('nan')
        return float(np.mean(self.free_paths))


    def std_free_path(self) -> float:
        '''
        Calculate the standard deviation of the free paths.
        '''
        if len(self.free_paths) == 0:
            return float('nan')
        return float(np.std(self.free_paths))
    
    def plot_histogram(self, bins: int = 100, show: bool = True) -> None:
        '''
        Plot histogram of free path lengths.

        Parameters:
        bins (int): Number of bins for the histogram.
        show (bool): Whether to display the plot immediately.
        '''
        if len(self.free_paths) == 0:
            print("No free paths recorded.")
            return
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.hist(self.free_paths, bins = bins, density = False, color = 'royalblue', histtype='step')
        ax.set_xlabel('Free Path Length [arbitrary units]')
        ax.set_ylabel('Entries')
        ax.set_yscale('log')
        if show:
            plt.show()
        return fig, ax
    
    def plot_trajectories(self, show: bool = True, draw_molecules: bool = True) -> None:
        '''
        Plot the trajectories of the particle and molecules.

        Parameters:
        show (bool): Whether to display the plot immediately.
        draw_molecules (bool): Whether to draw the molecules in their final positions.
        '''
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlim(0, self.box.width)
        ax.set_ylim(0, self.box.height)
        ax.set_aspect('equal', 'box')
        ax.set_xlabel('$x$ [arbitrary units]')
        ax.set_ylabel('$y$ [arbitrary units]')

        # draw molecules as circles at their current positions
        if draw_molecules:
            for m in self.molecules:
                circ = plt.Circle(tuple(m.pos), m.radius, alpha=0.5)
                ax.add_patch(circ)
        
        # draw particle
        circ_p = plt.Circle(tuple(self.particle.pos), self.particle.radius, color='red')
        ax.add_patch(circ_p)
        ax.set_title('Snapshot of particle + molecules')

        if show:
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()
        return fig, ax

    # (Optional) Animation of the process
    def animate(self, frames: int = 500, interval: int = 20, save_path: Optional[str] = None):
        # Make a lightweight animation showing positions; runs step() each frame
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_xlim(0, self.box.width)
        ax.set_ylim(0, self.box.height)
        ax.set_aspect('equal', adjustable='box')

        mol_patches = [plt.Circle(tuple(m.pos), m.radius, alpha=0.6) for m in self.molecules]
        for p in mol_patches:
            ax.add_patch(p)
        particle_patch = plt.Circle(tuple(self.particle.pos), self.particle.radius, color='red')
        ax.add_patch(particle_patch)

        def init():
            return []

        def update(frame):
            self.step()
            for p, m in zip(mol_patches, self.molecules):
                p.center = tuple(m.pos)
            particle_patch.center = tuple(self.particle.pos)
            ax.set_title(f'Step {self.steps_taken} Collisions: {self.collision_count}')
            return mol_patches + [particle_patch]

        anim = animation.FuncAnimation(fig, update, frames=frames, init_func=init, blit=False, interval=interval)
        if save_path is not None:
            anim.save(save_path, writer='ffmpeg')
        plt.show()
        return anim
    

#--------------------------------------------------------------
# Now let's do an example run
#--------------------------------------------------------------
if __name__ == '__main__':
    # Parameters (tweak as you like)
    box = Box(width=1.0, height=1.0)
    sim = Simulation(
        box = box,
        n_molecules = 200,
        molecule_radius = 0.01,
        molecule_thermal_kick_std = 2.0, # controls Brownian activity
        particle_speed = 2.0,
        particle_radius = 0.02,
        dt = 0.0005,
        max_steps = 200000,
    )

    print('Starting simulation...')
    sim.run(stop_after_collisions = 30000)

    print('Calculating mean free path...')
    print(f'Collisions recorded: {len(sim.free_paths)}')
    print(f'Mean free path: {sim.mean_free_path():.5f}')
    print(f'Std free path: {sim.std_free_path():.5f}')

    # Plot results
    print('Generating plots...')
    sim.plot_histogram(bins=100)
    sim.plot_trajectories()

    # Optional: animate (can be very slow)
    # sim.animate(frames=500, interval=30)

    print('Simulation complete.')

#--------------------------------------------------------------
# Think: What would the histogram be like, theoretically?
# Observe: Does the histogram look like what you expected?
# Improve: How can we improve the simulation? for example, can you 
#          include temperature, pressure, or other factors? Can you
#          make the code more concise by putting classes and methods
#          into separate files?
#--------------------------------------------------------------

# End of file. 