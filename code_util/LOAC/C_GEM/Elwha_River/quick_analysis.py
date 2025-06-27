"""
Quick Interactive Analysis Script for C-GEM Model Output
Simple script to quickly analyze specific variables
"""

import numpy as np
import matplotlib.pyplot as plt
from config import M, DELTI, DELXI
import os

def quick_plot(filename, variable_name="Variable"):
    """
    Quick plot of a single variable
    """
    if not os.path.exists(filename):
        print(f"File {filename} not found!")
        return
    
    # Load data
    data = np.fromfile(filename, dtype=np.float64)
    n_grid_points = M + 1
    n_time_steps = len(data) // n_grid_points
    data_reshaped = data[:n_time_steps * n_grid_points].reshape(n_time_steps, n_grid_points)
    
    # Create distance array
    distance_km = np.arange(M + 1) * DELXI / 1000
    
    # Create time array
    time_days = np.arange(n_time_steps) * DELTI / (24 * 3600)
    
    # Create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # Plot 1: Spatial distribution at different times
    times_to_plot = [0, n_time_steps//4, n_time_steps//2, -1]
    for t in times_to_plot:
        if t < 0:
            t = n_time_steps + t
        if t < n_time_steps:
            ax1.plot(distance_km, data_reshaped[t, :], 
                    label=f'Day {time_days[t]:.1f}', linewidth=2)
    
    ax1.set_xlabel('Distance from Mouth (km)')
    ax1.set_ylabel(f'{variable_name}')
    ax1.set_title(f'{variable_name} Spatial Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Time series at different locations
    locations = [0, M//4, M//2, 3*M//4, M]
    for loc in locations:
        if loc < n_grid_points:
            ax2.plot(time_days, data_reshaped[:, loc], 
                    label=f'{distance_km[loc]:.1f} km', linewidth=2)
    
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel(f'{variable_name}')
    ax2.set_title(f'{variable_name} Time Series')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Mean and std along estuary
    mean_data = np.mean(data_reshaped, axis=0)
    std_data = np.std(data_reshaped, axis=0)
    
    ax3.plot(distance_km, mean_data, 'b-', linewidth=2, label='Mean')
    ax3.fill_between(distance_km, mean_data - std_data, mean_data + std_data, 
                    alpha=0.3, color='blue', label='±1 Std Dev')
    ax3.set_xlabel('Distance from Mouth (km)')
    ax3.set_ylabel(f'{variable_name}')
    ax3.set_title(f'{variable_name} Statistics Along Estuary')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Heatmap (sampled)
    sample_rate = max(1, n_time_steps // 50)
    sampled_data = data_reshaped[::sample_rate, :]
    sampled_time = time_days[::sample_rate]
    
    im = ax4.imshow(sampled_data.T, aspect='auto', 
                   extent=[sampled_time[0], sampled_time[-1], distance_km[0], distance_km[-1]],
                   origin='lower', cmap='viridis')
    plt.colorbar(im, ax=ax4, label=f'{variable_name}')
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Distance from Mouth (km)')
    ax4.set_title(f'{variable_name} Heatmap')
    
    plt.tight_layout()
    plt.savefig(f'{variable_name}_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print statistics
    print(f"\n=== {variable_name} Statistics ===")
    print(f"Global mean: {np.mean(data_reshaped):.6f}")
    print(f"Global std: {np.std(data_reshaped):.6f}")
    print(f"Global min: {np.min(data_reshaped):.6f}")
    print(f"Global max: {np.max(data_reshaped):.6f}")
    print(f"Data shape: {data_reshaped.shape}")

def main():
    """
    Interactive main function
    """
    print("C-GEM Quick Analysis Tool")
    print("=" * 30)
    
    # List available files
    dat_files = [f for f in os.listdir('.') if f.endswith('.dat')]
    print(f"Available data files:")
    for i, filename in enumerate(dat_files):
        print(f"  {i+1}. {filename}")
    
    if not dat_files:
        print("No .dat files found in current directory!")
        return
    
    # Simple analysis for FCO2
    if 'FCO2.dat' in dat_files:
        print("\nAnalyzing FCO2.dat...")
        quick_plot('FCO2.dat', 'FCO2 (μmol/m³)')
    
    # Simple analysis for pH
    if 'pH.dat' in dat_files:
        print("\nAnalyzing pH.dat...")
        quick_plot('pH.dat', 'pH')
    
    # Simple analysis for DIC
    if 'DIC.dat' in dat_files:
        print("\nAnalyzing DIC.dat...")
        quick_plot('DIC.dat', 'DIC (μmol/kg)')
    
    print("\nAnalysis complete! Check the generated PNG files.")

if __name__ == "__main__":
    main() 