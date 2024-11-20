### Collisions.py 
### Ruize Li, Jesus College Cambridge, 20th Nov 2024. 

'''
This is a program that simulates the collisions of gas molecules (assuming they are spheres of constant density) in a square container. 
It also includes the brownian motion of the gas molecules and the collision of the gas molecules with the walls of the container.
Temperature was calculated (in 2D) and plotted in the animation.
'''

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import itertools

# Define the Canvas class
class Canvas():
    def __init__(self) -> None:
        self.size = np.array([20., 20.]) # in kilometers
        self.blocks = []
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot()

    def add_block(self, block) -> None:
        # This function adds a block to the canvas
        self.blocks.append(block)

    def update_blocks(self) -> None:
        # This function updates the position of all blocks on the canvas
        self.ax.clear()
        for i, block in enumerate(self.blocks):
            block.move()
            block.draw()

    def fix_axes(self) -> None:
        # This function fixes the axes and its labels of the canvas
        self.ax.set_xlim((-self.size[0]/2, self.size[0]/2))
        self.ax.set_ylim((-self.size[1]/2, self.size[1]/2))
        self.ax.set_aspect('equal')
        self.ax.set_ylabel('$y$ [km]')
        self.ax.set_xlabel('$x$ [km]')
        self.ax.set_title(f'Collision simulation, $T = ${self.temperature():.9f} K')
    
    def check_collisions(self) -> None:
        # This function checks for collisions between all pairs of blocks on the canvas
        # For more information, see the collide() function of the Block class
        combinations = list(itertools.combinations(range(len(self.blocks)), 2))
        for pair in combinations:
            self.blocks[pair[0]].collide(self.blocks[pair[1]])
        for block in self.blocks:
            block.bounce(canvas)
    
    def temperature(self) -> float:
        # This function calculates the temperature of the system
        k_B = 1.38064852e-23 # Boltzmann constant
        temperature = 0
        velocity_rms = np.sqrt(np.sum(np.linalg.norm(block.velocity)**2 for block in self.blocks) / len(self.blocks))
        mass_mean = np.mean([block.mass for block in self.blocks])
        temperature = mass_mean * 1E-22 * velocity_rms**2 / (2 * k_B)
        return temperature

# Define the Block class
class Block():
    def __init__(self, canvas, mass, position = np.array([0., 0.]), velocity = np.array([0., 0.])) -> None:
        self.canvas = canvas
        self.mass = mass # in 1E-25 kg
        self.position = position # in meters
        self.velocity = velocity # in kilometers per second
        self.canvas.add_block(self)
        self.color = 'black'
    
    def move(self) -> None:
        # This function updates the position of the block
        self.velocity = self.velocity + 0.01 * np.random.randn(2) - np.array([0.005, 0.005])
        self.position = self.position + self.velocity

    def draw(self) -> None:
        # This function draws the block on the canvas
        canvas.ax.plot(self.position[0], self.position[1], 'o')

    def collide(self, other) -> None:
        # This function checks for a collision between two blocks
        if np.linalg.norm(self.position - other.position) < 0.1:
            # change into self's frame
            velocity2 = other.velocity - self.velocity
            velocity1 = 0
            # collide in self's frame
            velocity1 = 2 * velocity2 * other.mass / (self.mass + other.mass)
            velocity2 = (self.mass - other.mass) * velocity2 / (self.mass + other.mass)
            # change back to lab frame
            other.velocity = velocity2 + self.velocity
            self.velocity = velocity1 + self.velocity

    def bounce(self, canvas) -> None:
        # This function checks for a bounce off the walls of the canvas
        # This is done by checking if the block will go outside the walls in the next time step
        for i in range(2):
            if (self.position[i] - canvas.size[i]/2) * (self.position[i] + self.velocity[i] - canvas.size[i]/2) < 0 or (self.position[i] + canvas.size[i]/2) * (self.position[i] + self.velocity[i] + canvas.size[i]/2) < 0:
                self.velocity[i] = -self.velocity[i]

# Create the canvas and add some blocks
canvas = Canvas()
blocks = []
for i in range(500):
    blocks.append(Block(canvas, mass = 1, position = 20*np.random.rand(2) - np.array([10., 10.]), velocity = np.random.rand(2) - np.array([0.5, 0.5])))

# Define the animation function
def animate(i) -> None:
    print("The frame is:", i)
    canvas.update_blocks()
    canvas.check_collisions()
    canvas.fix_axes()

# Create the animation object and save it as a gif
anim = animation.FuncAnimation(canvas.fig, animate, frames = 300, interval = 10)

# Save the animation as a gif
writervideo = animation.FFMpegWriter(fps=60)
anim.save('collisions.gif', writer='Pillow', fps=60)