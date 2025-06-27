"""
Comprehensive Data Analysis Script for C-GEM Model Output
Analyzes all .dat files and provides various visualization options
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from config import M, MAXT, DELTI, EL, DELXI
import os
import glob

# Define variable information for better plotting
VARIABLE_INFO = {
    'FCO2': {'unit': 'μmol/m³', 'description': 'Fugacity of CO2'},
    'NEM': {'unit': 'mmol C/m³/s', 'description': 'Net Ecosystem Metabolism'},
    'phy_death': {'unit': 'mmol C/m³/s', 'description': 'Phytoplankton Death'},
    'aer_deg': {'unit': 'mmol C/m³/s', 'description': 'Aerobic Degradation'},
    'denit': {'unit': 'mmol C/m³/s', 'description': 'Denitrification'},
    'nit': {'unit': 'mmol N/m³/s', 'description': 'Nitrification'},
    'O2_ex': {'unit': 'mmol O2/m³/s', 'description': 'O2 Exchange'},
    'NPP': {'unit': 'mmol C/m³/s', 'description': 'Net Primary Production'},
    'ALK': {'unit': 'μmol/kg', 'description': 'Alkalinity'},
    'pH': {'unit': '', 'description': 'pH'},
    'DIC': {'unit': 'μmol/kg', 'description': 'Dissolved Inorganic Carbon'},
    'S': {'unit': 'PSU', 'description': 'Salinity'},
    'SPM': {'unit': 'mg/L', 'description': 'Suspended Particulate Matter'},
    'O2': {'unit': 'μmol/kg', 'description': 'Dissolved Oxygen'},
    'PO4': {'unit': 'μmol/kg', 'description': 'Phosphate'},
    'TOC': {'unit': 'μmol/kg', 'description': 'Total Organic Carbon'},
    'NH4': {'unit': 'μmol/kg', 'description': 'Ammonium'},
    'DIA': {'unit': 'μmol/kg', 'description': 'Diatoms'},
    'dSi': {'unit': 'μmol/kg', 'description': 'Dissolved Silica'},
    'NO3': {'unit': 'μmol/kg', 'description': 'Nitrate'},
    'depth': {'unit': 'm', 'description': 'Water Depth'},
    'width': {'unit': 'm', 'description': 'Estuary Width'}
}

def load_data_file(filename):
    """
    Load data from any .dat file
    Returns: data array with shape (time_steps, grid_points)
    """
    if not os.path.exists(filename):
        print(f"Error: {filename} not found!")
        return None
    
    # Read the binary data
    data = np.fromfile(filename, dtype=np.float64)
    
    # Calculate dimensions
    n_grid_points = M + 1
    n_time_steps = len(data) // n_grid_points
    
    print(f"Loaded {filename}:")
    print(f"  Data shape: {len(data)} total values")
    print(f"  Grid points: {n_grid_points}")
    print(f"  Time steps: {n_time_steps}")
    print(f"  Simulation duration: {n_time_steps * DELTI / (24*3600):.1f} days")
    
    # Reshape data
    data_reshaped = data[:n_time_steps * n_grid_points].reshape(n_time_steps, n_grid_points)
    
    return data_reshaped

def create_distance_array():
    """Create distance array along the estuary"""
    return np.arange(M + 1) * DELXI

def plot_variable_comparison(data_dict, variables_to_plot=None, 
                           filename="variable_comparison.png"):
    """
    Plot multiple variables for comparison
    """
    if variables_to_plot is None:
        variables_to_plot = list(data_dict.keys())[:6]  # Plot first 6 variables
    
    n_vars = len(variables_to_plot)
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    distance_km = create_distance_array() / 1000
    
    for i, var_name in enumerate(variables_to_plot):
        if i >= len(axes):
            break
            
        data = data_dict[var_name]
        mean_data = np.mean(data, axis=0)
        
        ax = axes[i]
        ax.plot(distance_km, mean_data, 'b-', linewidth=2)
        ax.set_xlabel('Distance from Mouth (km)')
        
        unit = VARIABLE_INFO.get(var_name, {}).get('unit', '')
        ax.set_ylabel(f'{var_name} ({unit})')
        ax.set_title(f'{var_name}: {VARIABLE_INFO.get(var_name, {}).get("description", "")}')
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplots
    for i in range(n_vars, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def plot_correlation_matrix(data_dict, variables_to_analyze=None,
                          filename="correlation_matrix.png"):
    """
    Create a correlation matrix heatmap for selected variables
    """
    if variables_to_analyze is None:
        variables_to_analyze = list(data_dict.keys())[:8]  # First 8 variables
    
    # Calculate estuary-wide means for each variable
    means = {}
    for var_name in variables_to_analyze:
        data = data_dict[var_name]
        means[var_name] = np.mean(data, axis=1)  # Mean across space for each time step
    
    # Create correlation matrix
    df = pd.DataFrame(means)
    corr_matrix = df.corr()
    
    # Plot correlation matrix
    plt.figure(figsize=(12, 10))
    im = plt.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    plt.colorbar(im, label='Correlation Coefficient')
    
    # Add correlation values as text
    for i in range(len(corr_matrix.columns)):
        for j in range(len(corr_matrix.columns)):
            plt.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}', 
                    ha='center', va='center', fontsize=10)
    
    plt.xticks(range(len(corr_matrix.columns)), corr_matrix.columns, rotation=45)
    plt.yticks(range(len(corr_matrix.columns)), corr_matrix.columns)
    plt.title('Correlation Matrix of Model Variables')
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def plot_time_evolution(data_dict, variable_name, grid_points=[0, M//4, M//2, 3*M//4, M],
                       filename=None):
    """
    Plot time evolution of a specific variable
    """
    if filename is None:
        filename = f"{variable_name}_time_evolution.png"
    
    data = data_dict[variable_name]
    time_days = np.arange(data.shape[0]) * DELTI / (24 * 3600)
    
    plt.figure(figsize=(12, 8))
    
    for grid_point in grid_points:
        if grid_point < data.shape[1]:
            distance_km = grid_point * DELXI / 1000
            plt.plot(time_days, data[:, grid_point], 
                    label=f'Distance: {distance_km:.1f} km', linewidth=2)
    
    plt.xlabel('Time (days)')
    unit = VARIABLE_INFO.get(variable_name, {}).get('unit', '')
    plt.ylabel(f'{variable_name} ({unit})')
    plt.title(f'{variable_name}: {VARIABLE_INFO.get(variable_name, {}).get("description", "")}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def create_summary_report(data_dict, filename="model_summary_report.txt"):
    """
    Create a comprehensive summary report of all variables
    """
    with open(filename, 'w') as f:
        f.write("C-GEM Model Output Summary Report\n")
        f.write("=" * 50 + "\n\n")
        
        f.write(f"Model Configuration:\n")
        f.write(f"  Estuary Length: {EL/1000:.1f} km\n")
        f.write(f"  Grid Spacing: {DELXI} m\n")
        f.write(f"  Number of Grid Points: {M+1}\n")
        f.write(f"  Time Step: {DELTI} s\n")
        f.write(f"  Simulation Duration: {MAXT/(24*3600):.1f} days\n\n")
        
        f.write("Variable Statistics:\n")
        f.write("-" * 30 + "\n")
        
        for var_name, data in data_dict.items():
            f.write(f"\n{var_name}:\n")
            f.write(f"  Description: {VARIABLE_INFO.get(var_name, {}).get('description', 'N/A')}\n")
            f.write(f"  Unit: {VARIABLE_INFO.get(var_name, {}).get('unit', 'N/A')}\n")
            f.write(f"  Global Mean: {np.mean(data):.6f}\n")
            f.write(f"  Global Std: {np.std(data):.6f}\n")
            f.write(f"  Global Min: {np.min(data):.6f}\n")
            f.write(f"  Global Max: {np.max(data):.6f}\n")
            
            # Check for anomalies
            if np.any(np.isnan(data)):
                f.write(f"  WARNING: NaN values found!\n")
            if np.any(np.isinf(data)):
                f.write(f"  WARNING: Infinite values found!\n")

def main():
    """
    Main function to analyze all data files
    """
    print("C-GEM Model Data Analysis")
    print("=" * 30)
    
    # Find all .dat files
    dat_files = glob.glob("*.dat")
    print(f"Found {len(dat_files)} data files: {dat_files}")
    
    # Load all data
    data_dict = {}
    for filename in dat_files:
        var_name = filename.replace('.dat', '')
        print(f"\nLoading {filename}...")
        data = load_data_file(filename)
        if data is not None:
            data_dict[var_name] = data
    
    if not data_dict:
        print("No data files loaded!")
        return
    
    print(f"\nSuccessfully loaded {len(data_dict)} variables")
    
    # Create summary report
    print("\nCreating summary report...")
    create_summary_report(data_dict)
    
    # Create visualizations
    print("\nCreating visualizations...")
    
    # Variable comparison
    plot_variable_comparison(data_dict)
    
    # Correlation matrix
    plot_correlation_matrix(data_dict)
    
    # Time evolution for key variables
    key_variables = ['FCO2', 'pH', 'DIC', 'O2', 'S', 'NPP']
    for var_name in key_variables:
        if var_name in data_dict:
            plot_time_evolution(data_dict, var_name)
    
    print("\nAnalysis complete!")
    print("Generated files:")
    print("- model_summary_report.txt")
    print("- variable_comparison.png")
    print("- correlation_matrix.png")
    print("- [variable_name]_time_evolution.png files")

if __name__ == "__main__":
    main() 