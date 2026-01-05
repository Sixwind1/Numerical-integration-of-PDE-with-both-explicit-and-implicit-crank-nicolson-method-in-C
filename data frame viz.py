import numpy as np
import matplotlib.pyplot as plt
import os

def read_data(filename):
    """Read the data file and organize it by time steps"""
    data = {}
    current_time = None
    current_data = []
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('#'):
                # This is a time marker line
                if current_time is not None:
                    data[current_time] = np.array(current_data)
                parts = line.split()
                current_time = float(parts[2])  # Get the time value
                current_data = []
            else:
                # This is a data line (x, y, temperature)
                parts = list(map(float, line.split()))
                current_data.append(parts)
                
        # Add the last time step
        if current_time is not None and current_data:
            data[current_time] = np.array(current_data)
    
    return data

def visualize_temperature(data, time, output_filename):
    """Visualize temperature at a specific time and save to file"""
    if time not in data:
        print(f"Error: Time {time} not found in data. Available times: {sorted(data.keys())}")
        return
        
    # Extract data for the requested time
    time_data = data[time]
    
    # Create grid coordinates
    x_coords = np.unique(time_data[:, 0])
    y_coords = np.unique(time_data[:, 1])
    
    # Create temperature grid
    temp_grid = np.zeros((len(y_coords), len(x_coords)))
    
    # Fill the temperature grid
    for point in time_data:
        x, y, temp = point
        x_idx = np.where(x_coords == x)[0][0]
        y_idx = np.where(y_coords == y)[0][0]
        temp_grid[y_idx, x_idx] = temp
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    plt.imshow(temp_grid, cmap='hot', origin='lower', 
               extent=[x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()],
               vmin=0, vmax=70)  # Uniform scale from 0 to 70
    
    plt.colorbar(label='Temperature')
    plt.title(f'Temperature Distribution at t = {time}')
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    
    # Save the figure
    plt.savefig(output_filename)
    plt.close()
    print(f"Visualization saved as {output_filename}")

def main():
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_file = os.path.join(script_dir, 'dades.txt')
    
    # Read the data
    try:
        data = read_data(data_file)
    except Exception as e:
        print(f"Error reading data file: {e}")
        return
    
    # Get available times
    available_times = sorted(data.keys())
    print(f"Available times in data: {available_times}")
    
    # Get user input for time
    while True:
        try:
            time = float(input("Enter the time you want to visualize: "))
            if time in data:
                break
            else:
                print(f"Time {time} not found in data. Please choose from available times.")
        except ValueError:
            print("Please enter a valid number for time.")
    
    # Create output filename
    output_file = os.path.join(script_dir, f'temperature_at_time_{time}.png')
    
    # Visualize and save
    visualize_temperature(data, time, output_file)

if __name__ == "__main__":
    main()