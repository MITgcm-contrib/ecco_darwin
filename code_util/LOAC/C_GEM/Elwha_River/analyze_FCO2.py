"""
FCO2 Data Analysis and Visualization Script
Analyzes the FCO2.dat file from the C-GEM model output
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from config import M, MAXT, DELTI, EL, DELXI
import os

def load_FCO2_data(filename="FCO2.dat"):
    """
    Load FCO2 data from the .dat file
    Returns: data array with shape (time_steps, grid_points)
    """
    if not os.path.exists(filename):
        print(f"Error: {filename} not found!")
        return None
    
    # Read the binary data
    data = np.fromfile(filename, dtype=np.float64)
    
    # Calculate dimensions based on model configuration
    n_grid_points = M + 1  # Number of grid points
    n_time_steps = len(data) // n_grid_points
    
    print(f"Data shape: {len(data)} total values")
    print(f"Grid points: {n_grid_points}")
    print(f"Time steps: {n_time_steps}")
    print(f"Simulation duration: {n_time_steps * DELTI / (24*3600):.1f} days")
    
    # Reshape data to (time_steps, grid_points)
    data_reshaped = data[:n_time_steps * n_grid_points].reshape(n_time_steps, n_grid_points)
    
    return data_reshaped

def create_distance_array():
    """Create distance array along the estuary"""
    return np.arange(M + 1) * DELXI

def plot_time_series(data, grid_points_to_plot=[0, M//4, M//2, 3*M//4, M], 
                    filename="FCO2_time_series.png"):
    """
    Plot FCO2 time series at selected grid points
    """
    plt.figure(figsize=(12, 8))
    
    time_days = np.arange(data.shape[0]) * DELTI / (24 * 3600)  # Convert to days
    
    for i, grid_point in enumerate(grid_points_to_plot):
        if grid_point < data.shape[1]:
            distance_km = grid_point * DELXI / 1000  # Convert to km
            plt.plot(time_days, data[:, grid_point], 
                    label=f'Distance: {distance_km:.1f} km', linewidth=2)
    
    plt.xlabel('Time (days)')
    plt.ylabel('FCO2 (μmol/m³)')
    plt.title('FCO2 Time Series at Different Estuary Locations')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def plot_spatial_distribution(data, time_steps_to_plot=[0, -1], 
                            filename="FCO2_spatial.png"):
    """
    Plot FCO2 spatial distribution at selected time steps
    """
    plt.figure(figsize=(12, 8))
    
    distance_km = create_distance_array() / 1000  # Convert to km
    
    for i, time_step in enumerate(time_steps_to_plot):
        if time_step < 0:
            time_step = data.shape[0] + time_step
        if time_step < data.shape[0]:
            time_days = time_step * DELTI / (24 * 3600)
            plt.plot(distance_km, data[time_step, :], 
                    label=f'Day {time_days:.1f}', linewidth=2)
    
    plt.xlabel('Distance from Mouth (km)')
    plt.ylabel('FCO2 (μmol/m³)')
    plt.title('FCO2 Spatial Distribution Along Estuary')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def plot_heatmap(data, filename="FCO2_heatmap.png"):
    """
    Create a heatmap showing FCO2 variation over time and space
    """
    plt.figure(figsize=(14, 8))
    
    # Sample data for heatmap (every 10th time step to avoid memory issues)
    sample_rate = max(1, data.shape[0] // 100)
    sampled_data = data[::sample_rate, :]
    
    time_days = np.arange(0, data.shape[0], sample_rate) * DELTI / (24 * 3600)
    distance_km = create_distance_array() / 1000
    
    # Create the heatmap
    im = plt.imshow(sampled_data.T, aspect='auto', 
                   extent=[time_days[0], time_days[-1], distance_km[0], distance_km[-1]],
                   origin='lower', cmap='RdBu_r')
    
    plt.colorbar(im, label='FCO2 (μmol/m³)')
    plt.xlabel('Time (days)')
    plt.ylabel('Distance from Mouth (km)')
    plt.title('FCO2 Heatmap: Time vs Distance')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def plot_statistics(data, filename="FCO2_statistics.png"):
    """
    Plot statistical summary of FCO2 data
    """
    distance_km = create_distance_array() / 1000
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # Mean FCO2 along estuary
    mean_fco2 = np.mean(data, axis=0)
    ax1.plot(distance_km, mean_fco2, 'b-', linewidth=2)
    ax1.set_xlabel('Distance from Mouth (km)')
    ax1.set_ylabel('Mean FCO2 (μmol/m³)')
    ax1.set_title('Mean FCO2 Along Estuary')
    ax1.grid(True, alpha=0.3)
    
    # Standard deviation
    std_fco2 = np.std(data, axis=0)
    ax2.plot(distance_km, std_fco2, 'r-', linewidth=2)
    ax2.set_xlabel('Distance from Mouth (km)')
    ax2.set_ylabel('Std Dev FCO2 (μmol/m³)')
    ax2.set_title('FCO2 Variability Along Estuary')
    ax2.grid(True, alpha=0.3)
    
    # Min and Max
    min_fco2 = np.min(data, axis=0)
    max_fco2 = np.max(data, axis=0)
    ax3.fill_between(distance_km, min_fco2, max_fco2, alpha=0.3, color='green')
    ax3.plot(distance_km, mean_fco2, 'k-', linewidth=2, label='Mean')
    ax3.set_xlabel('Distance from Mouth (km)')
    ax3.set_ylabel('FCO2 (μmol/m³)')
    ax3.set_title('FCO2 Range Along Estuary')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Time series of estuary-wide mean
    time_days = np.arange(data.shape[0]) * DELTI / (24 * 3600)
    estuary_mean = np.mean(data, axis=1)
    ax4.plot(time_days, estuary_mean, 'purple', linewidth=2)
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Mean FCO2 (μmol/m³)')
    ax4.set_title('Estuary-wide Mean FCO2 Over Time')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def analyze_data(data):
    """
    Print basic statistics about the FCO2 data
    """
    print("\n=== FCO2 Data Analysis ===")
    print(f"Data shape: {data.shape}")
    print(f"Total values: {data.size}")
    print(f"Global mean: {np.mean(data):.4f} μmol/m³")
    print(f"Global std: {np.std(data):.4f} μmol/m³")
    print(f"Global min: {np.min(data):.4f} μmol/m³")
    print(f"Global max: {np.max(data):.4f} μmol/m³")
    
    # Check for any NaN or infinite values
    if np.any(np.isnan(data)):
        print("Warning: NaN values found in data!")
    if np.any(np.isinf(data)):
        print("Warning: Infinite values found in data!")

def main():
    """
    Main function to run the analysis
    """
    print("Loading FCO2 data...")
    data = load_FCO2_data()
    
    if data is None:
        return
    
    # Analyze the data
    analyze_data(data)
    
    # Create visualizations
    print("\nCreating visualizations...")
    
    # Time series plot
    plot_time_series(data)
    
    # Spatial distribution plot
    plot_spatial_distribution(data)
    
    # Heatmap
    plot_heatmap(data)
    
    # Statistics
    plot_statistics(data)
    
    print("\nAnalysis complete! Check the generated PNG files for visualizations.")

if __name__ == "__main__":
    main() 