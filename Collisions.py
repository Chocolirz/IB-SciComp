### Collisions.py 
### Ruize Li, Jesus College Cambridge, 20th Nov 2024. 

'''
This is a program that simulates the collisions of gas molecules (assuming they are spheres of constant density) in a cubic container. 
'''
# Importing necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import itertools

class Canvas():
    def __init__(self) -> None:
        self.size = 20
        self.blocks = []
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot()

    def add_block(self, block) -> None:
        self.blocks.append(block)

    def update_blocks(self) -> None:
        self.ax.clear()
        for i, block in enumerate(self.blocks):
            block.move()
            block.draw()

    def fix_axes(self) -> None:
        self.ax.set_xlim((-self.size/2, self.size/2))
        self.ax.set_ylim((-1, 1))
    
    def check_collisions(self) -> None:
        combinations = list(itertools.combinations(range(len(self.blocks)), 2))
        for pair in combinations:
            self.blocks[pair[0]].collide(self.blocks[pair[1]])

class Block():
    def __init__(self, canvas, mass, position = 0, velocity = 0) -> None:
        self.canvas = canvas
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.canvas.add_block(self)
        self.color = 'black'
    
    def move(self) -> None:
        self.position += self.velocity

    def draw(self) -> None:
        canvas.ax.plot(self.position, 0, 'o')

    def collide(self, other) -> None:
        if abs(self.position - other.position) < 0.1:
            self.velocity *= -1
            other.velocity *= -1

canvas = Canvas()
block1 = Block(canvas, mass = 1, position = -2, velocity = 0.07)
block2 = Block(canvas, mass = 1, position = 2, velocity = -0.07)
block3 = Block(canvas, mass = 1, position = 4, velocity = -0.05)
block4 = Block(canvas, mass = 1, position = -5, velocity = 0.05)

def animate(i) -> None:
    print("The frame is:", i)
    canvas.update_blocks()
    canvas.check_collisions()
    canvas.fix_axes()

anim = animation.FuncAnimation(canvas.fig, animate, frames = 500, interval = 10)

writervideo = animation.FFMpegWriter(fps=30)
anim.save('collisions.gif', writer='imagemagick', fps=60)