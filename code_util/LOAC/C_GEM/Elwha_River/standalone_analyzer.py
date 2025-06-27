"""
Standalone C-GEM Data Analyzer
Can analyze any .dat file from any directory
Usage: 
1. Command line: python standalone_analyzer.py <path/to/filename.dat>
2. Direct call: analyze_dat_file("path/to/filename.dat")
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def analyze_dat_file(filepath):
    """
    Analyze any .dat file from C-GEM model output
    """
    if not os.path.exists(filepath):
        print(f"Error: {filepath} not found!")
        return None
    
    # Read the binary data
    print(f"Loading {filepath}...")
    data = np.fromfile(filepath, dtype=np.float64)
    
    # Try to determine grid size from file size
    # Common grid sizes for C-GEM models
    possible_grid_sizes = [40, 80, 160, 320, 640, 1280]  # Common M+1 values
    
    n_grid_points = None
    for grid_size in possible_grid_sizes:
        if len(data) % grid_size == 0:
            n_grid_points = grid_size
            break
    
    if n_grid_points is None:
        # If no exact match, try to estimate
        n_grid_points = int(np.sqrt(len(data)))
        print(f"Warning: Could not determine exact grid size. Using estimate: {n_grid_points}")
    
    n_time_steps = len(data) // n_grid_points
    
    print(f"Data shape: {len(data)} total values")
    print(f"Estimated grid points: {n_grid_points}")
    print(f"Time steps: {n_time_steps}")
    
    # Reshape data
    data_reshaped = data[:n_time_steps * n_grid_points].reshape(n_time_steps, n_grid_points)
    
    # Create distance array (assuming typical DELXI = 2000m)
    DELXI = 2000  # Default grid spacing in meters
    distance_km = np.arange(n_grid_points) * DELXI / 1000
    
    # Create time array (assuming typical DELTI = 150s)
    DELTI = 150  # Default time step in seconds
    time_days = np.arange(n_time_steps) * DELTI / (24 * 3600)
    
    # Get variable name from filename
    variable_name = os.path.splitext(os.path.basename(filepath))[0]
    
    # Create visualizations
    create_plots(data_reshaped, distance_km, time_days, variable_name, filepath)
    
    # Print statistics
    print_statistics(data_reshaped, variable_name)
    
    return data_reshaped, distance_km, time_days, variable_name

def create_plots(data, distance_km, time_days, variable_name, original_filepath):
    """
    Create comprehensive plots for the data
    """
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    
    # Plot 1: Spatial distribution at different times
    ax1 = plt.subplot(2, 3, 1)
    times_to_plot = [0, len(time_days)//4, len(time_days)//2, -1]
    for t in times_to_plot:
        if t < 0:
            t = len(time_days) + t
        if t < len(time_days):
            ax1.plot(distance_km, data[t, :], 
                    label=f'Day {time_days[t]:.1f}', linewidth=2)
    
    ax1.set_xlabel('Distance from Mouth (km)')
    ax1.set_ylabel(f'{variable_name}')
    ax1.set_title(f'{variable_name} Spatial Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Time series at different locations
    ax2 = plt.subplot(2, 3, 2)
    n_grid = data.shape[1]
    locations = [0, n_grid//4, n_grid//2, 3*n_grid//4, n_grid-1]
    for loc in locations:
        if loc < n_grid:
            ax2.plot(time_days, data[:, loc], 
                    label=f'{distance_km[loc]:.1f} km', linewidth=2)
    
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel(f'{variable_name}')
    ax2.set_title(f'{variable_name} Time Series')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Mean and std along estuary
    ax3 = plt.subplot(2, 3, 3)
    mean_data = np.mean(data, axis=0)
    std_data = np.std(data, axis=0)
    
    ax3.plot(distance_km, mean_data, 'b-', linewidth=2, label='Mean')
    ax3.fill_between(distance_km, mean_data - std_data, mean_data + std_data, 
                    alpha=0.3, color='blue', label='±1 Std Dev')
    ax3.set_xlabel('Distance from Mouth (km)')
    ax3.set_ylabel(f'{variable_name}')
    ax3.set_title(f'{variable_name} Statistics Along Estuary')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Heatmap (sampled)
    ax4 = plt.subplot(2, 3, 4)
    sample_rate = max(1, len(time_days) // 50)
    sampled_data = data[::sample_rate, :]
    sampled_time = time_days[::sample_rate]
    
    im = ax4.imshow(sampled_data.T, aspect='auto', 
                   extent=[sampled_time[0], sampled_time[-1], distance_km[0], distance_km[-1]],
                   origin='lower', cmap='viridis')
    plt.colorbar(im, ax=ax4, label=f'{variable_name}')
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Distance from Mouth (km)')
    ax4.set_title(f'{variable_name} Heatmap')
    
    # Plot 5: Min/Max range
    ax5 = plt.subplot(2, 3, 5)
    min_data = np.min(data, axis=0)
    max_data = np.max(data, axis=0)
    
    ax5.fill_between(distance_km, min_data, max_data, alpha=0.3, color='green', label='Range')
    ax5.plot(distance_km, mean_data, 'k-', linewidth=2, label='Mean')
    ax5.set_xlabel('Distance from Mouth (km)')
    ax5.set_ylabel(f'{variable_name}')
    ax5.set_title(f'{variable_name} Range Along Estuary')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Time series of estuary-wide mean
    ax6 = plt.subplot(2, 3, 6)
    estuary_mean = np.mean(data, axis=1)
    ax6.plot(time_days, estuary_mean, 'purple', linewidth=2)
    ax6.set_xlabel('Time (days)')
    ax6.set_ylabel(f'Mean {variable_name}')
    ax6.set_title(f'Estuary-wide Mean {variable_name} Over Time')
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot in the same directory as the original file
    output_dir = os.path.dirname(original_filepath)
    if output_dir == '':
        output_dir = '.'
    output_filename = os.path.join(output_dir, f'{variable_name}_analysis.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Saved plot as: {output_filename}")
    
    plt.show()

def print_statistics(data, variable_name):
    """
    Print comprehensive statistics
    """
    print(f"\n=== {variable_name} Statistics ===")
    print(f"Data shape: {data.shape}")
    print(f"Total values: {data.size}")
    print(f"Global mean: {np.mean(data):.6f}")
    print(f"Global std: {np.std(data):.6f}")
    print(f"Global min: {np.min(data):.6f}")
    print(f"Global max: {np.max(data):.6f}")
    
    # Check for anomalies
    if np.any(np.isnan(data)):
        print("Warning: NaN values found in data!")
    if np.any(np.isinf(data)):
        print("Warning: Infinite values found in data!")
    
    # Spatial statistics
    spatial_mean = np.mean(data, axis=0)
    spatial_std = np.std(data, axis=0)
    print(f"\nSpatial Statistics:")
    print(f"  Mean along estuary: {np.mean(spatial_mean):.6f}")
    print(f"  Std along estuary: {np.mean(spatial_std):.6f}")
    
    # Temporal statistics
    temporal_mean = np.mean(data, axis=1)
    temporal_std = np.std(data, axis=1)
    print(f"\nTemporal Statistics:")
    print(f"  Mean over time: {np.mean(temporal_mean):.6f}")
    print(f"  Std over time: {np.mean(temporal_std):.6f}")

def main():
    """
    Main function - can be called with command line argument or used interactively
    """
    if len(sys.argv) > 1:
        # Command line usage
        filepath = sys.argv[1]
        analyze_dat_file(filepath)
    else:
        # Interactive usage - you can modify these paths directly
        print("No command line argument provided.")
        print("You can:")
        print("1. Call analyze_dat_file() directly with a file path")
        print("2. Modify the filepath below and run the script")
        
        # Example usage - modify this path as needed
        filepath = r"C:\Users\xaric\OneDrive\Documents\visual studio code projects\python\nasa_jpl_cgem\ecco_darwin\code_util\LOAC\C_GEM\code_python_FCO2\FCO2.dat"
        
        # Uncomment the line below to run with the example path
        # analyze_dat_file(filepath)

if __name__ == "__main__":
    main() 