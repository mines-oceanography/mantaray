'''ode solution plotting

A script to plot the output from the ode_solvers, stepper.y_out(), after
stepper.integrate() is called.

This script currently only works and was tested for space separated files, where
each row is the state of the wave at the given timestep. The columns  of the
file are t, x, y, kx, ky.

This script requires matplotlib and numpy.

The path to the file is set by changing line 22, and the option to use
animations is on line 44. But both of these will change to be command line
arguments in the future.
'''

import matplotlib.pyplot as plt
import matplotlib.animation as mpla
import numpy as np

# Read the values from the file
with open('../y_out.txt', 'r') as file:
    data = file.readlines()

t = []
x = []
y = []
kx = []
ky = []

for i, line in enumerate(data):
    if i == 0: continue
    values = line.strip().split()
    t.append(float(values[0]))
    x.append(float(values[1]))
    y.append(float(values[2]))
    kx.append(float(values[3]))
    ky.append(float(values[4]))

# Create a figure and axis object
fig, ax = plt.subplots()

# If wanting an animation, change below to True:
animating = True

# Graph the vectors
ax.quiver(x, y, kx, ky, color='red')

if not animating:
    # Plot x(t) and y(t) on the axis
    ax.plot(x, y, label='x(t), y(t)', linewidth=2)
    ax.scatter(x, y)
else:
    # Create animation of points x(t),y(t) and line
    line, = ax.plot([], [], label="x(t),y(t)")
    scatter = ax.scatter([], [])

    def animate(i):
        line.set_data(x[:i],y[:i])
        scatter.set_offsets(np.column_stack((x[:i], y[:i])))
        return line, scatter

    anim = mpla.FuncAnimation(
        plt.gcf(), animate, interval=100, frames=len(t)+1, repeat=True, repeat_delay=500,
    )

# Set the scale on the horizontal and vertical axes
ax.set_xlim(min(x) - 1, max(x) + 1)
ax.set_ylim(min(y) - 1, max(y) + 1)

# Customize the graph
ax.set_xlabel('x(t)')
ax.set_ylabel('y(t)')
ax.set_title('Graph of x(t) and y(t)')
ax.legend()

# Turn on minor ticks
plt.minorticks_on()

# Display the graph
plt.show()
