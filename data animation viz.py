import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LinearSegmentedColormap

# Read the data from the file
with open('dades.txt', 'r') as file:
    data = file.readlines()

# Parse the data into time steps and grids
time_steps = []
grids = []
current_grid = []

for line in data:
    line = line.strip()
    if line.startswith('# t'):
        if current_grid:
            grids.append(np.array(current_grid))
            current_grid = []
        time = float(line.split()[2])  # Extract time value
        time_steps.append(time)
    else:
        if line:
            parts = list(map(float, line.split()))
            current_grid.append(parts)

if current_grid:
    grids.append(np.array(current_grid))

# Determine grid dimensions
nx = len(np.unique(grids[0][:, 0]))  # Unique x values
ny = len(np.unique(grids[0][:, 1]))  # Unique y values

# Reshape each grid into a 2D array
temperature_grids = []
for grid in grids:
    T = grid[:, 2]
    T_grid = T.reshape(ny, nx)  # Reshape to (y, x) for imshow
    temperature_grids.append(T_grid)

# Set fixed temperature scale (0 to 60°C)
vmin, vmax = 0, 60

# Create a high-resolution colormap (e.g., 'jet' or 'viridis')
cmap = plt.get_cmap('jet', 256)  # 256 color levels for smooth gradients

# Create the animation
fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(temperature_grids[0], cmap=cmap, origin='lower', 
               extent=[0, 1, 0, 1], vmin=vmin, vmax=vmax)
cbar = plt.colorbar(im, label='Temperature (°C)', extend='both')
ax.set_xlabel('X Position')
ax.set_ylabel('Y Position')
ax.set_title(f'Time: {time_steps[0]:.2f}')

def update(frame):
    im.set_array(temperature_grids[frame])
    ax.set_title(f'Time: {time_steps[frame]:.2f}')
    return im,

# Speed up the animation (interval=50ms per frame)
ani = FuncAnimation(fig, update, frames=len(temperature_grids), 
                    interval=50, blit=False)

plt.tight_layout()
plt.show()

# Save the animation (uncomment if needed)
# ani.save('temperature_animation.mp4', writer='ffmpeg', fps=20, dpi=300)