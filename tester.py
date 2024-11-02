import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import math

def naca_4_digit_airfoil(m, p, t, c, num_points=100):
    # Generate x coordinates from 0 to c
    x = np.linspace(0, c, num_points)
    
    # Thickness distribution formula for a symmetrical 4-digit NACA airfoil
    yt = 5 * t * c * (0.2969 * np.sqrt(x/c) - 0.1260 * (x/c) - 0.3516 * (x/c)**2 + 0.2843 * (x/c)**3 - 0.1015 * (x/c)**4)

    # Camber line and its slope are zero for NACA 0012 (m=0, p=0)
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)

    # Angle of the camber line
    theta = np.arctan(dyc_dx)

    # Upper surface coordinates
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)

    # Lower surface coordinates
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)

    return xu, yu, xl, yl

def cutoff_trailing_edge_symmetrical(x_upper, y_upper, x_lower, y_lower, cutoff_ratio, c):
    # Calculate the cutoff x-coordinate
    cutoff_x = c * (1 - cutoff_ratio)
    
    # Interpolate to find the y-coordinates at the cutoff x-coordinate
    y_upper_cutoff = np.interp(cutoff_x, x_upper, y_upper)
    y_lower_cutoff = np.interp(cutoff_x, x_lower, y_lower)
    
    # Find the indices where to cut off
    upper_cutoff_index = np.where(x_upper >= cutoff_x)[0][0]
    lower_cutoff_index = np.where(x_lower >= cutoff_x)[0][0]
    
    # Include the exact cutoff point
    x_upper_cutoff = np.append(x_upper[:upper_cutoff_index], cutoff_x)
    y_upper_cutoff = np.append(y_upper[:upper_cutoff_index], y_upper_cutoff)
    x_lower_cutoff = np.append(x_lower[:lower_cutoff_index], cutoff_x)
    y_lower_cutoff = np.append(y_lower[:lower_cutoff_index], y_lower_cutoff)

    return x_upper_cutoff, y_upper_cutoff, x_lower_cutoff, y_lower_cutoff

def tangent_line_at_cutoff(x_cutoff, y_cutoff):
    delta = 1e-16
    x1 = x_cutoff[-1] - delta
    y1 = np.interp(x1, x_cutoff, y_cutoff)
    x2 = x_cutoff[-1] + delta
    y2 = np.interp(x2, x_cutoff, y_cutoff)
    
    dy_dx = (y2 - y1) / (x2 - x1)
    
    # Coordinates of the cutoff point
    x0 = x1
    y0 = y1
    
    return dy_dx, x0, y0

def perpendicular_line_at_cutoff(dy_dx, x0, y0):
    # Slope of the perpendicular line
    dy_dx_perp = -1 / dy_dx
    
    return dy_dx_perp, x0, y0

def find_intersection(dy_dx1, x0_1, y0_1, dy_dx2, x0_2, y0_2):
    # Solve the system of linear equations to find the intersection point
    A = np.array([[-dy_dx1, 1], [-dy_dx2, 1]])
    B = np.array([y0_1 - dy_dx1 * x0_1, y0_2 - dy_dx2 * x0_2])
    intersection_point = np.linalg.solve(A, B)
    
    return intersection_point

def distance(point1, point2):
    return np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

# Define parameters for NACA 0012
m = 0.0     # Maximum camber
p = 0.0     # Location of maximum camber
t = 0.12    # Maximum thickness as a fraction of the chord
c = 1.0     # Chord length
cutoff_ratio = 0.1

# Generate airfoil coordinates
x_upper, y_upper, x_lower, y_lower = naca_4_digit_airfoil(m, p, t, c)

# Cut off at 0.9 of the chord length symmetrically
x_upper_cutoff, y_upper_cutoff, x_lower_cutoff, y_lower_cutoff = cutoff_trailing_edge_symmetrical(x_upper, y_upper, x_lower, y_lower, cutoff_ratio, c)

# Find the tangent lines at the cutoff points
dy_dx_upper, x0_upper, y0_upper = tangent_line_at_cutoff(x_upper_cutoff, y_upper_cutoff)
dy_dx_lower, x0_lower, y0_lower = tangent_line_at_cutoff(x_lower_cutoff, y_lower_cutoff)

# Find the perpendicular lines at the cutoff points
dy_dx_perp_upper, x0_perp_upper, y0_perp_upper = perpendicular_line_at_cutoff(dy_dx_upper, x0_upper, y0_upper)
dy_dx_perp_lower, x0_perp_lower, y0_perp_lower = perpendicular_line_at_cutoff(dy_dx_lower, x0_lower, y0_lower)

# Find the intersection point of the two perpendicular lines
intersection = find_intersection(dy_dx_perp_upper, x0_perp_upper, y0_perp_upper, dy_dx_perp_lower, x0_perp_lower, y0_perp_lower)

# Print the equations of the tangent and perpendicular lines
tangent_upper_eq = f"y = {dy_dx_upper:.4f} * (x - {[x_upper_cutoff[-1]]}) + {[y_upper_cutoff[-1]]}"
tangent_lower_eq = f"y = {dy_dx_lower:.4f} * (x - {[x_lower_cutoff[-1]]}) + {[y_lower_cutoff[-1]]}"
perpendicular_upper_eq = f"y = {dy_dx_perp_upper:.4f} * (x - {[x_upper_cutoff[-1]]}) + {[y_upper_cutoff[-1]]}"
perpendicular_lower_eq = f"y = {dy_dx_perp_lower:.4f} * (x - {[x_lower_cutoff[-1]]}) + {[y_lower_cutoff[-1]]}"
top_cutoff_point = (x_upper_cutoff[-1], y_upper_cutoff[-1])
distance_to_top_cutoff = distance(intersection, top_cutoff_point)

start_point = (x_upper_cutoff[-1], y_upper_cutoff[-1])
end_point   = (x_lower_cutoff[-1], y_lower_cutoff[-1])

arcsine_radians_top = math.asin(dy_dx_upper)
arcsine_radians_bot = math.asin(dy_dx_lower)
start_angle = np.degrees(np.arctan2(start_point[1] - intersection[1], start_point[0] - intersection[0]))
end_angle   = np.degrees(np.arctan2(end_point[1] - intersection[1], end_point[0] - intersection[0]))
r = distance_to_top_cutoff
arc = Arc((end_point[0], end_point[1]), 2*r, 2*r, angle=0, theta1=start_angle, theta2=end_angle, color='blue')
plt.gca().add_patch(arc)

print("Equation of the tangent line at the upper surface cutoff point:", tangent_upper_eq)
print("Equation of the tangent line at the lower surface cutoff point:", tangent_lower_eq)
print("Equation of the perpendicular line at the upper surface cutoff point:", perpendicular_upper_eq)
print("Equation of the perpendicular line at the lower surface cutoff point:", perpendicular_lower_eq)
print(f"Intersection point of the two perpendicular lines: {intersection}")
print(f"Distance from the intersection point to the top of the cutoff: {distance_to_top_cutoff}")



# Ensure arc starts and ends at the correct points
#arc_x = np.concatenate([[x_upper_cutoff[-1]], arc_x, [x_lower_cutoff[-1]]])
#arc_y = np.concatenate([[y_upper_cutoff[-1]], arc_y, [y_lower_cutoff[-1]]])

# Create the combined coordinates
#x_combined = np.concatenate([x_upper_cutoff, arc_x, x_lower_cutoff[::-1]])
#y_combined = np.concatenate([y_upper_cutoff, arc_y, y_lower_cutoff[::-1]])

# Plot the original and modified airfoil
plt.figure(figsize=(10, 5))
plt.plot(x_upper, y_upper, label='Original NACA 0012 Upper Surface')
plt.plot(x_lower, y_lower, label='Original NACA 0012 Lower Surface')
#plt.plot(x_combined, y_combined, label='Modified NACA 0012 with Arc (10% Cutoff)', linestyle='--')

plt.xlabel('x')
plt.ylabel('y')
plt.title('NACA 0012 Airfoil with 10% Trailing Edge Cutoff and Arc Connection')
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
