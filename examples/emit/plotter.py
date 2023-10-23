import numpy as np
import matplotlib.pyplot as plt
#

def low_pass_filter(data, cutoff_frequency):
    # Perform FFT
    transformed = np.fft.fft(data)
    frequency = np.fft.fftfreq(len(data))

    # Set high-frequency components to zero
    transformed[(frequency > cutoff_frequency)] = 0

    # Perform inverse FFT to get filtered data
    filtered_data = np.fft.ifft(transformed)
    
    return np.real(filtered_data)  # Take only the real part of the data after IFFT


def parse_file(filename):
    results = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("ITEM: TIMESTEP"):
                timestep = int(lines[i+1].strip())
                number_of_cells = int(lines[i+3].strip())
                xcs = []
                ycs  = []
                density_values_oii = []
                density_values_w  = []
                for j in range(number_of_cells):
                    cell_data_index = i + 5 + j
                    if cell_data_index >= len(lines):  # Check to avoid going out of bounds
                        break
                    cell_data_line = lines[cell_data_index]
                    cell_data = cell_data_line.split()
                    if len(cell_data) < 4:  # Ensure we have enough data in the line
                        continue
                    try:
                        xc = float(cell_data[1])
                        yc = float(cell_data[2])
                        density_oii = float(cell_data[3])  # 10 O+ -1 W
                        print(xc, yc, density_oii)
                        density_w = float(cell_data[-1])  # 10 O+ -1 W
                        xcs.append(xc)
                        ycs.append(yc)
                        density_values_oii.append(density_oii)
                        density_values_w.append(density_w)
                    except ValueError:
                        # Skip if not a valid float
                        continue
                results[timestep] = (xcs, ycs, density_values_oii, density_values_w)
    return results



#import matplotlib.pyplot as plt
#
data = parse_file("tmp.grid")


last_timestep = max(data.keys())
if last_timestep != 0:  # Ensure it's not the zero timestep you wanted to skip
    print(last_timestep)

    xcs, ycs, densities_oii, density_w = data[last_timestep]
    mask = [i for i, xc_val in enumerate(xcs) if xc_val > 2.6]

    # Apply this mask to all other lists
    xcs0 = [xcs[i] for i in mask]
    ycs0 = [ycs[i] for i in mask]
    densities_oii_ = [densities_oii[i] for i in mask]
    density_w_ = [density_w[i] for i in mask]

#    fig, (ax1, ax2, ax3, ax4) = plt.subplots(2, 2, figsize=(6, 10))
    fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # Sort ycs and densities based on ycs values
    sorted_indices = sorted(range(len(ycs0)), key=lambda k: ycs0[k])
    ycs_sorted = [ycs0[i] for i in sorted_indices]
    densities_oii_sorted = [densities_oii_[i] for i in sorted_indices]
    densities_w_sorted = [density_w_[i] for i in sorted_indices]

#    ax1.plot(ycs_sorted, densities_oii_sorted, 'b')  # 'ko-' means black circles connected with lines


    # Filter out zero densities for oii
    non_zero_indices_oii = [i for i, density in enumerate(densities_oii_sorted) if density != 0]
    ycs_oii_non_zero = [ycs_sorted[i] for i in non_zero_indices_oii]
    densities_oii_non_zero = [densities_oii_sorted[i] for i in non_zero_indices_oii]

    # Filter out zero densities for w
    non_zero_indices_w = [i for i, density in enumerate(densities_w_sorted) if density != 0]
    ycs_w_non_zero = [ycs_sorted[i] for i in non_zero_indices_w]
    densities_w_non_zero = [densities_w_sorted[i] for i in non_zero_indices_w]

    cutoff_frequency = 0.1  # example value; adjust based on your needs

    # Apply the filter to the data
    densities_oii_filtered = low_pass_filter(densities_oii_non_zero, cutoff_frequency)
    densities_w_filtered = low_pass_filter(densities_w_non_zero, cutoff_frequency)


    ax1.plot(ycs_oii_non_zero, densities_oii_filtered, 'k', marker='o')  # black lines connecting non-zero points
#    ax3.plot(ycs_w_non_zero, densities_w_non_zero, 'k', marker='o')  # red lines connecting non-zero points for densities_w (assuming you want to plot it in red)


#    ax1.set_ylim(1e10, max(densities_oii_non_zero)*1.1)  # Uncomment if you need this xlim
#    ax1.set_xlim(-1,0.75)

    ax1.set_title(f'OII: {last_timestep}' )
    ax1.set_xlabel('Z [m]')
    ax1.set_ylabel('Density [1/m$^3$]')
    ax1.grid(alpha=0.3)

    # 2D Plot
    # Convert the data to a structured 2D grid
    x = np.unique(xcs)
    y = np.unique(ycs)
    density_grid_oii = np.zeros((len(y), len(x)))

    for i, (xc, yc, density_oii) in enumerate(zip(xcs, ycs, densities_oii)):
        j = np.where(x == xc)[0]
        k = np.where(y == yc)[0]
        density_grid_oii[k, j] = density_oii

#    pcm = ax2.pcolormesh(x, y, density_grid_oii, shading='auto')  # 'shading' is set to 'auto' to eliminate warning in newer matplotlib versions
    pcm = ax2.pcolormesh(x, y, density_grid_oii, shading='auto', cmap='plasma')  # You can replace 'viridis' with other colormaps like 'plasma', 'inferno', 'magma', etc.
#    ax2.set_xlim([x.min(), x.])
#    ax2.set_ylim([ymin, ymax])
    ax2.tick_params(direction='out', length=6, width=2, colors='black', grid_color='gray', grid_alpha=0.5)


#    fig.colorbar(pcm, ax=ax2, label='Density [1/m$^3$]')
    cbar = fig.colorbar(pcm, ax=ax2, label='Density [1/m$^{-3}$]')
    cbar.ax.tick_params(labelsize=10)  # Change tick font size

    ax2.set_xlabel('R [m]')
    ax2.set_ylabel('Z [m]')


    plt.tight_layout()
    plt.show()
