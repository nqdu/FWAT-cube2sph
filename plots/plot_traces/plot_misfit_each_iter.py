import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from fwat.const import SRC_REC, MISFIT
import glob

# Get max iteration from command line argument
def main():
    if len(sys.argv) != 3:
        print("Usage: python plot_noise_misfit.py simu_type max_iteration ")
        sys.exit(1)

    simu_type = sys.argv[1]
    max_iter = int(sys.argv[2])
    path = "../.."
    # Read parameter file to get period bands

    # Get event list
    evts = np.loadtxt(f"{path}/{SRC_REC}/sources.dat.{simu_type}", usecols=0, dtype=str)
    # Find all unique bands across all iterations (turning points detection)
    all_bands_by_iter = {}  # Store bands available at each iteration
    for iter_num in range(0, max_iter + 1):
        mod = f"M{iter_num:02d}"
        pattern = f"../../{MISFIT}/{mod}/{evts[0]}_*_{simu_type}_window_chi"
        files = glob.glob(pattern)
        # Extract band names from filenames
        bands_at_iter = []
        for f in files:
            basename = os.path.basename(f)
            # Parse: eventname_bandname_simu_type_window_chi
            parts = basename.split('_noise_window_chi')[0].split('_')
            band = parts[1] + "_" + parts[2]  # Extract band portion
            bands_at_iter.append(band)
        all_bands_by_iter[iter_num] = bands_at_iter
        # print(f"Iteration {iter_num}: Found {len(bands_at_iter)} bands bands = {bands_at_iter}")
    
    # Detect turning points (where band set changes)
    turning_points = [0]  # Always start at iteration 0
    prev_bands = all_bands_by_iter[0]
    for iter_num in range(1, max_iter + 1):
        curr_bands = all_bands_by_iter[iter_num]
        if curr_bands != prev_bands:
            turning_points.append(iter_num)
            # print(f"Turning point detected at iteration {iter_num}")
            prev_bands = curr_bands
    turning_points.append(max_iter + 1)  # Add end point for range

    print("Detected {} multi-scale stages".format(len(turning_points)-1))

    # Calculate misfit for each iteration and band
    misfit_all = {'iterations': [], 'misfits': []}
    for iter_num in range(0, max_iter + 1):
        mod = f"M{iter_num:02d}"
        
        total_misfit = 0.
        total_n = 0
        # Determine bands to consider based on turning points
        all_bands = all_bands_by_iter[iter_num]
        for ib, band in enumerate(all_bands):
            sumf = 0.
            sumn = 0
            
            for ievt in range(len(evts)):
                filename = f"../../{MISFIT}/{mod}/{evts[ievt]}_{band}_{simu_type}_window_chi"
                if os.path.exists(filename):
                    d = np.loadtxt(filename, usecols=28, ndmin=2)
                    chi = np.sum(d[:, 0])
                    if not np.isnan(chi):
                        sumf += chi
                        sumn += d.shape[0]

            total_misfit += sumf
            total_n += sumn

        if total_n > 0:
            misfit_all['iterations'].append(iter_num)
            misfit_all['misfits'].append(total_misfit)

    # Create the plot
    plt.figure(figsize=(10, 6))

    # generate colors for each range
    colors = plt.get_cmap('viridis')(np.linspace(0, 1, len(turning_points)-1))

    # Plot total misfit for each band
    # normalize misfit by using the turning points
    npts = len(turning_points) - 1
    for ib in range(npts):
        iter_start = turning_points[ib]
        iter_end = turning_points[ib + 1]
        misfits_in_this_range = misfit_all['misfits'][iter_start:iter_end]
        iters_in_this_range = misfit_all['iterations'][iter_start:iter_end]
        bands = all_bands_by_iter[iter_start]

        # only get the first/last band in this range for plotting
        band = bands[0].split('_')[0] + "-" + bands[-1].split('_')[1]

        if ib > 0:
            iters_in_this_range = np.array(iters_in_this_range) - 1 # adjust x-axis to make sure overlap at turning point

        # normalize
        misfits_in_this_range = np.array(misfits_in_this_range)  / misfits_in_this_range[0]
        plt.plot(iters_in_this_range, misfits_in_this_range, 
                 'o-', label=f'Band {band} (iter {iter_start}-{iter_end-1})', linewidth=2, markersize=6)

    plt.xlabel('Iteration', fontsize=12)
    plt.ylabel('Normalized Misfit', fontsize=12)
    plt.title(f'Misfit by Band and Iteration (up to iter {max_iter})', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'noise_misfit.png', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    main()