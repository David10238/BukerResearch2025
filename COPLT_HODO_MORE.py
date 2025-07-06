import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from scipy.interpolate import CubicSpline
import time


from matplotlib.patches import Ellipse
barb_heights = None

# Function to generate the full hodograph using a spline
# Create the Tkinter GUI

# root.title("Hodograph Generator with SRH and Superhelicity")

# # Create sliders for u_start and v_start
# tk.Label(root, text="Starting U (m/s):").grid(row=8, column=0)
# u_start_slider = tk.Scale(root, from_=-20, to=20, resolution=0.5, orient=tk.HORIZONTAL)
# u_start_slider.set(0)  # Default starting U
# u_start_slider.grid(row=8, column=1)

# tk.Label(root, text="Starting V (m/s):").grid(row=8, column=2)
# v_start_slider = tk.Scale(root, from_=-20, to=20, resolution=0.5, orient=tk.HORIZONTAL)
# v_start_slider.set(0)  # Default starting V
# v_start_slider.grid(row=8, column=3)
# Global variable to alternate barb colors
alternate_barb_color = False

# Function to generate the full hodograph using a spline
def generate_hodograph(curvatures, lengths):
    # Get starting U and V from sliders
    u_start = u_start_slider.get()
    v_start = v_start_slider.get()

    # Reduce lengths by 70%
    lengths = [length * 0.5 for length in lengths]
    
    # Generate control points for the spline
    u_points = [u_start]  # Start at the user-defined U
    v_points = [v_start]  # Start at the user-defined V
    for i in range(len(curvatures)):
        angle = np.pi / 2 * curvatures[i]  # Map curvature to angle
        u_points.append(u_points[-1] + lengths[i] * np.cos(angle))
        v_points.append(v_points[-1] + lengths[i] * np.sin(angle))
    
    # Create a cubic spline for smooth interpolation
    heights = np.linspace(0, 5, len(u_points))  # Heights corresponding to each segment
    spline_u = CubicSpline(heights, u_points, bc_type='natural')
    spline_v = CubicSpline(heights, v_points, bc_type='natural')
    
    # Interpolate 100 points per segment
    interpolated_heights = np.linspace(0, 5, 500)
    u_interpolated = spline_u(interpolated_heights)
    v_interpolated = spline_v(interpolated_heights)
    
    return np.column_stack((u_interpolated, v_interpolated))
# Function to plot the hodograph


# Function to calculate Bunkers storm motion
def calculate_bunkers_motion(surface_wind, wind_5km):
    """
    Calculate Bunkers storm motion for right-moving and left-moving supercells.
    
    Parameters:
        surface_wind (tuple): Surface wind vector (u, v) in m/s.
        wind_5km (tuple): 5 km wind vector (u, v) in m/s.
    
    Returns:
        dict: Right-moving and left-moving storm motion vectors.
    """
    # Mean wind
    mean_wind = ((surface_wind[0] + wind_5km[0]) / 2, (surface_wind[1] + wind_5km[1]) / 2)
    
    # Shear vector
    shear_vector = (wind_5km[0] - surface_wind[0], wind_5km[1] - surface_wind[1])
    
    # Perpendicular shear vector (rotate 90 degrees to the right)
    perp_shear_vector = (-shear_vector[1], shear_vector[0])
    
    # Scale the perpendicular shear vector (7.5 m/s is typical)
    scale_factor = 7.5 / np.sqrt(perp_shear_vector[0]**2 + perp_shear_vector[1]**2)
    scaled_perp_shear = (perp_shear_vector[0] * scale_factor, perp_shear_vector[1] * scale_factor)
    
    # Right-moving storm motion
    right_motion = (mean_wind[0] - scaled_perp_shear[0], mean_wind[1] - scaled_perp_shear[1])
    
    # Left-moving storm motion
    left_motion = (mean_wind[0] + scaled_perp_shear[0], mean_wind[1] + scaled_perp_shear[1])
    
    return {"Right Motion": right_motion, "Left Motion": left_motion, "Mean Wind": mean_wind}

## Function to calculate SRH
def calculate_srh(heights, u, v, storm_motion, depth=3000):
    """
    Calculate Storm-Relative Helicity (SRH) over a specified depth.
    
    Parameters:
        heights (np.ndarray): Heights (in meters) above ground level.
        u (np.ndarray): Zonal wind components (m/s).
        v (np.ndarray): Meridional wind components (m/s).
        storm_motion (tuple): Storm motion vector (u_storm, v_storm) in m/s.
        depth (int): Depth (in meters) over which to calculate SRH (default: 3000 m).
    
    Returns:
        float: Storm-Relative Helicity (SRH) in m^2/s^2.
    """
    # Limit the profile to the specified depth
    mask = heights <= depth
    heights = heights[mask]
    u = u[mask]
    v = v[mask]
    
    # Calculate storm-relative wind
    u_relative = u - storm_motion[0]
    v_relative = v - storm_motion[1]
    
    # Calculate vertical shear
    du_dz = np.diff(u) / np.diff(heights)
    dv_dz = np.diff(v) / np.diff(heights)
    
    # Calculate horizontal vorticity
    vorticity_x = -dv_dz  # ∂v/∂z
    vorticity_y = du_dz  # -∂u/∂z
    
    # Calculate the dot product of storm-relative wind and horizontal vorticity
    srh = 0
    for i in range(len(du_dz)):
        dot_product = (u_relative[i] * vorticity_x[i] + v_relative[i] * vorticity_y[i])
        srh += dot_product * (heights[i + 1] - heights[i])  # Approximate integral
    
    return srh

# Function to calculate Superhelicity
def calculate_superhelicity(heights, u, v, depth=1000):
    """
    Calculate Superhelicity (S) over a specified depth.
    
    Parameters:
        heights (np.ndarray): Heights (in meters) above ground level.
        u (np.ndarray): Zonal wind components (m/s).
        v (np.ndarray): Meridional wind components (m/s).
        depth (int): Depth (in meters) over which to calculate Superhelicity (default: 3000 m).
    
    Returns:
        float: Superhelicity (S) in s^-2.
    """
    # Limit the profile to the specified depth
    mask = heights <= depth
    heights = heights[mask]
    u = u[mask]
    v = v[mask]
    
    # Calculate vertical shear (first derivatives)
    du_dz = np.diff(u) / np.diff(heights)
    dv_dz = np.diff(v) / np.diff(heights)
    
    # Calculate shear of shear (second derivatives)
    d2u_dz2 = np.diff(du_dz) / np.diff(heights[:-1])
    d2v_dz2 = np.diff(dv_dz) / np.diff(heights[:-1])
    
    # Calculate superhelicity
    superhelicity = 0
    for i in range(len(d2u_dz2)):
        term1 = dv_dz[i] * d2u_dz2[i]
        term2 = du_dz[i] * d2v_dz2[i]
        superhelicity += (term1 + term2) * (heights[i + 2] - heights[i + 1])  # Approximate integral
    
    return superhelicity

# Function to update the plot based on slider values
def update_plot(*args):
    # Get the current values from the sliders
    curvatures = [curvature_sliders[i].get() for i in range(5)]
    lengths = [length_sliders[i].get() for i in range(5)]
    
    
    # Generate the updated hodograph
    hodograph = generate_hodograph(curvatures, lengths)
    
    # Get the selected cflag value
    cflag = cflag_var.get()
    eflag = eflag_var.get()
    max_radius = max_radius_slider.get()  # Get the max radius
    
    # Re-plot the hodograph with the updated data and selected cflag
    plot_hodograph(hodograph, canvas, ax, srh_label, cflag, eflag, max_radius)


import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import Normalize




def plot_hodograph(hodograph, canvas, ax, srh_label, cflag=1.0,eflag=0,max_radius=40):
    global barb_heights,barb_indices,u_barbs,v_barbs

    # Clear the current axes and reset the plot to avoid accumulation
    ax.clear()

    # Plot the hodograph with smaller markers
    ax.plot(hodograph[:, 0], hodograph[:, 1], marker='o', markersize=2, label="Hodograph", zorder=1)  # Tiny markers with lower z-order

    # Add radii contours (like a polar plot)
    radii_interval = 10  # Interval between radii
    for radius in range(radii_interval, max_radius + 1, radii_interval):
        circle = plt.Circle((0, 0), radius, color='gray', fill=False, linestyle='dotted', linewidth=0.8, zorder=0)
        ax.add_artist(circle)
        ax.text(radius, 0, f"{radius} m/s", color="gray", fontsize=8, ha="left", va="center")  # Annotate the radii

    # Set the plot limits based on max_radius
    ax.set_xlim(-max_radius, max_radius)
    ax.set_ylim(-max_radius, max_radius)

    # Calculate heights and wind components
    heights = np.linspace(0, 5000, len(hodograph))  # Assume heights from 0 to 5 km
    u = hodograph[:, 0]
    v = hodograph[:, 1]





    # Plot the wind barbs
    #ax2.barbs(np.zeros_like(barb_heights), barb_heights, u_barbs, v_barbs, length=6, color=barb_color, zorder=3)


    # Redraw the canvas
    canvas.draw()

    
    # Clear the current axes and reset the plot to avoid accumulation
    ax.clear()
    
    # Plot the hodograph with smaller markers
    ax.plot(hodograph[:, 0], hodograph[:, 1], marker='o', markersize=2, label="Hodograph", zorder=1)  # Tiny markers with lower z-order
    
    # Add radii contours (like a polar plot)
    #max_radius = 40 # Maximum wind speed to display
    radii_interval = 10  # Interval between radii
    for radius in range(radii_interval, max_radius + 1, radii_interval):
        circle = plt.Circle((0, 0), radius, color='gray', fill=False, linestyle='dotted', linewidth=0.8, zorder=0)
        ax.add_artist(circle)
        ax.text(radius, 0, f"{radius} m/s", color="gray", fontsize=8, ha="left", va="center")  # Annotate the radii
    
    # Calculate heights and wind components
    heights = np.linspace(0, 5000, len(hodograph))  # Assume heights from 0 to 5 km
    u = hodograph[:, 0]
    v = hodograph[:, 1]

        # Calculate and plot Bunkers storm motion
    surface_wind = hodograph[0]  # 0 km wind
    wind_5km = hodograph[-1]     # 5 km wind
    storm_motion = calculate_bunkers_motion(surface_wind, wind_5km)
    mean_wind = storm_motion["Mean Wind"]
    right_motion = storm_motion["Right Motion"]
    left_motion = storm_motion["Left Motion"]

    # Calculate storm-relative wind
    u_rel = u - right_motion[0]
    v_rel = v - right_motion[1]
    
    # Annotate specific heights (0 km, 1 km, ..., 5 km)
    num_segments = 5
    segment_indices = np.linspace(0, len(u) - 1, num_segments + 1, dtype=int)  # Indices for 0 km to 5 km
    for i, idx in enumerate(segment_indices):
        ax.text(u[idx], v[idx], f"{i} km", fontsize=6, color="blue", ha="left", va="bottom")  # Annotate the height
    
    # Calculate vertical shear (first derivatives)
    du_dz = np.diff(u) / np.diff(heights)
    dv_dz = np.diff(v) / np.diff(heights)
    
    # Calculate shear of shear (second derivatives)
    d2u_dz2 = np.diff(du_dz) / np.diff(heights[:-1])
    d2v_dz2 = np.diff(dv_dz) / np.diff(heights[:-1])
    
    # Calculate helicity contribution

    helicity_contribution = du_dz[:-1] * v_rel[1:-1] - dv_dz[:-1] * u_rel[1:-1]
    superhelicity_contribution = du_dz[:-1] * d2v_dz2 + dv_dz[:-1] * d2u_dz2
    
    # Calculate superhelicity contribution
    sharpening_contribution = helicity_contribution * np.sqrt(du_dz[:-1]**2 + dv_dz[:-1]**2)
    
    # Calculate local magnitude of vorticity
    vorticity_magnitude = np.sqrt(du_dz**2 + dv_dz**2)
    
    # Determine coloring values based on cflag
    if cflag == 0:
        coloring_values = helicity_contribution
    elif cflag == 1:
        coloring_values = superhelicity_contribution
    elif cflag == 2:
        coloring_values = vorticity_magnitude[:-1]
    elif cflag == 3:
        coloring_values = sharpening_contribution
    else:
        raise ValueError("Invalid cflag value. Must be 0, 1, 2 or 3.")
    
    # Normalize the coloring values for the colormap
    cmap = plt.get_cmap('coolwarm')  # Choose a colormap (e.g., 'coolwarm', 'viridis', etc.)
    norm = Normalize(vmin=np.min(coloring_values), vmax=np.max(coloring_values))  # Normalize coloring values
    
    # Plot shear vectors, colored by the selected contribution
    for i in range(len(hodograph) - 1):  # Plot vectors for each segment
        dx = hodograph[i + 1, 0] - hodograph[i, 0]
        dy = hodograph[i + 1, 1] - hodograph[i, 1]
        shear_magnitude = np.sqrt(dx**2 + dy**2)  # Scale by shear vector magnitude
        
        # Get the color corresponding to the selected contribution
        if i < len(coloring_values):  # Ensure we don't exceed the array bounds
            color = cmap(norm(coloring_values[i]))
        else:
            color = 'gray'  # Default color for the last segment
        
        # Tangent vector (scaled, length doubled)
        ax.arrow(hodograph[i, 0], hodograph[i, 1], dx * 0.1 / shear_magnitude, dy * 0.1 / shear_magnitude,
                 linewidth=4, head_width=0.15, head_length=0.15, fc=color, ec=color, zorder=2)  # Color by contribution
        
            # EYES
        
        # Add cartoon "eyes" to the plot
        if eflag == 1:  # Only plot eyes if eflag is set to 1
            eye_index = -20  # 20th-last point of the spline
            if len(hodograph) > abs(eye_index):  # Ensure the index is valid
                eye_center = hodograph[eye_index]
                eye_offset = 1.5  # Distance between the eyes (in m/s)

                # Left eye (football-shaped)
                left_eye = Ellipse((eye_center[0] - eye_offset, eye_center[1]), width=1.5, height=0.75, zorder=5, facecolor='#99FFDD',edgecolor='black')
                ax.add_artist(left_eye)

                # Right eye (football-shaped)
                right_eye = Ellipse((eye_center[0] + eye_offset, eye_center[1]), width=1.5, height=0.75, zorder=5, facecolor='#99FFDD',edgecolor='black')
                ax.add_artist(right_eye)

                # Left pupil
                left_pupil = plt.Circle((eye_center[0] - eye_offset, eye_center[1]), radius=0.2, color='black', zorder=6)
                ax.add_artist(left_pupil)

                # Right pupil
                right_pupil = plt.Circle((eye_center[0] + eye_offset, eye_center[1]), radius=0.2, color='black', zorder=6)
                ax.add_artist(right_pupil)
        

    
    # Plot mean wind
    ax.arrow(0, 0, mean_wind[0], mean_wind[1], linewidth=2, head_width=0.3, head_length=0.3, fc='green', ec='green', label="Mean Wind", zorder=3)
    
    # Plot right-moving storm motion
    ax.arrow(0, 0, right_motion[0], right_motion[1], linewidth=2, head_width=0.3, head_length=0.3, fc='purple', ec='purple', label="Right Motion", zorder=3)
    
    # Plot left-moving storm motion
    ax.arrow(0, 0, left_motion[0], left_motion[1], linewidth=2, head_width=0.3, head_length=0.3, fc='orange', ec='orange', label="Left Motion", zorder=3)
    
    # Annotate storm motion attributes
    ax.text(right_motion[0], right_motion[1], "Right Motion", color='purple', fontsize=10)
    ax.text(left_motion[0], left_motion[1], "Left Motion", color='orange', fontsize=10)
    ax.text(mean_wind[0], mean_wind[1], "Mean Wind", color='green', fontsize=10)

    superh=calculate_superhelicity(heights, u, v, depth=1000)
    # Update SRH and Superhelicity label
    srh_label.config(text=f"SRH (Right): {calculate_srh(heights, u, v, right_motion, depth=3000):.2f} m²/s²\n"
                          f"SRH (Left): {calculate_srh(heights, u, v, left_motion, depth=3000):.2f} m²/s²\n"
                          f"Superhelicity (0–1 km): {superh * 1.e+5:.2f} s⁻²\n"
                          f"Superhelicity (1–2 km): {(calculate_superhelicity(heights, u, v, depth=2000)-superh) * 1.e+5:.2f} s⁻²")
    
    ax.set_xlabel("U (m/s)")
    ax.set_ylabel("V (m/s)")
    ax.set_title("HODO-PYTHON")
    ax.grid(True)
    #ax.legend()
    ax.set_xlim(-max_radius, max_radius)  # Adjust limits to fit the radii
    ax.set_ylim(-max_radius, max_radius)

    # Create GridSpec for the plot and colorbar layout
    gs = gridspec.GridSpec(2, 1)  # 10x1 ratio: one large and one small axis for colorbar
    ax = plt.subplot(gs[0])  # Main plot axis
    cax = plt.subplot(gs[1])  # Colorbar axis
    cax.axis("off")  # Turn off the main plot's axes
    ax.axis("off")  # Turn off the main plot's axes





    
    # Remove existing colorbar if it exists
    # Remove existing colorbar if it exists
    if hasattr(plot_hodograph, 'colorbar') and plot_hodograph.colorbar:
        #print("Removing existing colorbar")
        plot_hodograph.colorbar.remove()
        plot_hodograph.colorbar = None  # Reset the colorbar attribute
        canvas.draw()  # Force the canvas to redraw without the colorbar

    # Add the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plot_hodograph.colorbar = fig.colorbar(sm, ax=cax, orientation='vertical', pad=.02,shrink=0.3)
    plot_hodograph.colorbar.set_label('Color Scale for Contributions', fontsize=8)


    # Improve font rendering for colorbar tick labels
    for label in plot_hodograph.colorbar.ax.get_yticklabels():
        label.set_fontsize(8)  # Adjust font size
        #label.set_fontname('Arial')  # Use a smoother font


     # Add wind barbs on a vertical axis

    #ax2.barbs(np.zeros_like(barb_heights), barb_heights, u_barbs, v_barbs, length=6, color = "white", zorder=3)
        

    



    # Create a secondary vertical axis for wind barbs
    ax2 = ax.twinx()  # Create a twin axis sharing the same x-axis
    ax2.set_ylim(-10, 1)  # Set the vertical axis range from 0 to 5 km
    ax2.axis("off")  # Turn off the axis decorations
    # if previous_u_barbs is not None and previous_v_barbs is not None and previous_barb_heights is not None:
    #     ax2 = ax.twinx()  # Create a twin axis sharing the same x-axis
    #     ax2.set_ylim(0, 5)  # Set the vertical axis range from 0 to 5 km
    #     ax2.axis("off")  # Turn off the axis decorations
    #     print("white?")
    #     ax2.barbs(np.zeros_like(previous_barb_heights), previous_barb_heights, previous_u_barbs, previous_v_barbs,
    #               length=6, color='white', zorder=3)

    # Plot the new wind barbs
    if barb_heights is not None:
        ax2.barbs(np.ones_like(barb_heights), barb_heights-5, u_barbs, v_barbs, linewidth=3,length=6, color = "white", zorder=3)
    barb_heights = np.arange(0, 5.25, 0.25)  # Heights every 0.25 km
    barb_indices = np.linspace(0, len(u) - 1, len(barb_heights), dtype=int)  # Indices corresponding to barb heights
    u_barbs = u[barb_indices]
    v_barbs = v[barb_indices]
    #ax2.barbs(np.zeros_like(barb_heights), barb_heights, u_barbs, v_barbs, length=6, color = "black", zorder=3)
        # Plot each barb individually, making every 5th barb red
    for i, height in enumerate(barb_heights):
        color = "red" if i % 4 == 0 else "black"  # Every 5th barb is red
        ax2.barbs([1], [height-5], [u_barbs[i]], [v_barbs[i]], length=6, color=color, zorder=3)

 


    # Redraw the canvas
    canvas.draw()


def reset_sliders():
    # Reset curvature sliders to 0
    for slider in curvature_sliders:
        slider.set(0)
    # Reset length sliders to 10
    for slider in length_sliders:
        slider.set(10)
    max_radius_slider.set(40)
    # Update the plot after resetting
    update_plot()








# Create the Tkinter GUI
root = tk.Tk()
root.title("THE HODO-PYTHON: Hodograph Generator with SRH and Superhelicity")
cflag_var = tk.IntVar(value=1)

# Add radio buttons to select cflag
tk.Label(root, text="Select Contribution:").grid(row=10, column=0, sticky="w")  # Label for the radio buttons

tk.Radiobutton(root, text="Helicity", variable=cflag_var, value=0, command=update_plot).grid(row=11, column=0, sticky="w")
tk.Radiobutton(root, text="Superhelicity", variable=cflag_var, value=1, command=update_plot).grid(row=12, column=0, sticky="w")
tk.Radiobutton(root, text="Vorticity Magnitude", variable=cflag_var, value=2, command=update_plot).grid(row=13, column=0, sticky="w")
tk.Radiobutton(root, text="Vorticity Sharpening", variable=cflag_var, value=3, command=update_plot).grid(row=14, column=0, sticky="w")
# Create a Tkinter IntVar to store the selected eflag value
eflag_var = tk.IntVar(value=0)  # Default to 0 (off)

# Add radio buttons to toggle eflag
tk.Label(root, text="Enable Eyes:").grid(row=12, column=1, sticky="w")  # Label for the radio buttons

tk.Radiobutton(root, text="Off", variable=eflag_var, value=0, command=update_plot).grid(row=13, column=1, sticky="w")
tk.Radiobutton(root, text="On", variable=eflag_var, value=1, command=update_plot).grid(row=14, column=1, sticky="w")

#
# Add a reset button to the GUI
reset_button = tk.Button(root, text="Reset Sliders", command=reset_sliders, bg="lightblue", font=("Arial", 12, "bold"))
reset_button.grid(row=15, column=0, columnspan=2, pady=10)  # Adjust the position as needed

tk.Label(root, text="Max Radius:").grid(row=15, column=2, sticky="w")  # Label for the slider
max_radius_slider = tk.Scale(root, from_=10, to=50, resolution=1, orient=tk.HORIZONTAL, command=update_plot)
max_radius_slider.set(40)  # Default max radius
max_radius_slider.grid(row=15, column=2, sticky="w")  # Place the slider next to the reset button

# Create the Matplotlib figure and canvas
fig, ax = plt.subplots(figsize=(8, 8))  # Make the plot taller
canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.grid(row=0, column=0, rowspan=10)  # Place the plot in the first column, spanning multiple rows

# Create sliders for curvature and length
curvature_sliders = []
length_sliders = []

for i in range(5):
    # Curvature slider
    tk.Label(root, text=f"Segment {i+1} Curvature:").grid(row=i * 2, column=1, sticky="w")  # Place in column 1
    curvature_slider = tk.Scale(root, from_=-2, to=2, resolution=0.1, orient=tk.HORIZONTAL, command=update_plot)
    curvature_slider.set(0)  # Default curvature
    curvature_slider.grid(row=i * 2, column=2, sticky="w")  # Place in column 2
    curvature_sliders.append(curvature_slider)
    
    # Length slider
    tk.Label(root, text=f"Segment {i+1} Length:").grid(row=i * 2 + 1, column=1, sticky="w")  # Place in column 1
    length_slider = tk.Scale(root, from_=5, to=20, resolution=0.5, orient=tk.HORIZONTAL, command=update_plot)
    length_slider.set(10)  # Default length
    length_slider.grid(row=i * 2 + 1, column=2, sticky="w")  # Place in column 2
    length_sliders.append(length_slider)

# Create sliders for starting U and V
tk.Label(root, text="Starting U (m/s):").grid(row=10, column=1, sticky="w")  # Place in column 1
u_start_slider = tk.Scale(root, from_=-20, to=20, resolution=0.5, orient=tk.HORIZONTAL, command=update_plot)
u_start_slider.set(0)  # Default starting U
u_start_slider.grid(row=10, column=2, sticky="w")  # Place in column 2

tk.Label(root, text="Starting V (m/s):").grid(row=11, column=1, sticky="w")  # Place in column 1
v_start_slider = tk.Scale(root, from_=-20, to=20, resolution=0.5, orient=tk.HORIZONTAL, command=update_plot)
v_start_slider.set(0)  # Default starting V
v_start_slider.grid(row=11, column=2, sticky="w")  # Place in column 2



# Add SRH and Superhelicity label
srh_label = tk.Label(root, text="SRH (Right): 0.00 m²/s²\n"
                                 "SRH (Left): 0.00 m²/s²\n"
                                 "Superhelicity (0–1 km): 0.00 s⁻²\n"
                                 "Superhelicity (0–3 km): 0.00 s⁻²\n"
                                 "Sum (0–1 km): 0.00 (s⁻³)(m⁻¹)\n"
                                 "Sum (0–3 km): 0.00 (s⁻³)(m⁻¹)")#, font=("Arial", 14, "bold"))
srh_label.grid(row=10, column=0, columnspan=4)


#root.mainloop()




# Create the Tkinter GUI
# root = tk.Tk()
# root.title("Hodograph Generator with SRH and Superhelicity")



# Initialize the plot
initial_curvatures = [slider.get() for slider in curvature_sliders]
initial_lengths = [slider.get() for slider in length_sliders]
initial_hodograph = generate_hodograph(initial_curvatures, initial_lengths)
plot_hodograph(initial_hodograph, canvas, ax, srh_label)



# Run the Tkinter event loop
root.mainloop()