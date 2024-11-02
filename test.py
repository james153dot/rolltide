import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Arc

def naca_4_digit_airfoil(m, p, t, c, num_points=100):
    x = np.linspace(0, c, num_points)
    yt = 5 * t * c * (0.2969 * np.sqrt(x/c) - 0.1260 * (x/c) - 0.3516 * (x/c)**2 + 0.2843 * (x/c)**3 - 0.1015 * (x/c)**4)
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)
    theta = np.arctan(dyc_dx)
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    return xu, yu, xl, yl

def cutoff_trailing_edge_symmetrical(x_upper, y_upper, x_lower, y_lower, cutoff_ratio, c):
    cutoff_x = c * (1 - cutoff_ratio)
    y_upper_cutoff = np.interp(cutoff_x, x_upper, y_upper)
    y_lower_cutoff = np.interp(cutoff_x, x_lower, y_lower)
    upper_cutoff_index = np.where(x_upper >= cutoff_x)[0][0]
    lower_cutoff_index = np.where(x_lower >= cutoff_x)[0][0]
    x_upper_cutoff = np.append(x_upper[:upper_cutoff_index], cutoff_x)
    y_upper_cutoff = np.append(y_upper[:upper_cutoff_index], y_upper_cutoff)
    x_lower_cutoff = np.append(x_lower[:lower_cutoff_index], cutoff_x)
    y_lower_cutoff = np.append(y_lower[:lower_cutoff_index], y_lower_cutoff)
    return x_upper_cutoff, y_upper_cutoff, x_lower_cutoff, y_lower_cutoff

def tangent_line_at_cutoff(x_cutoff, y_cutoff):
    delta = 1e-8
    x1 = x_cutoff[-1] - delta
    y1 = np.interp(x1, x_cutoff, y_cutoff)
    x2 = x_cutoff[-1] + delta
    y2 = np.interp(x2, x_cutoff, y_cutoff)
    dy_dx = (y2 - y1) / (x2 - x1)
    x0 = x1
    y0 = y1
    return dy_dx, x0, y0

def perpendicular_line_at_cutoff(dy_dx, x0, y0):
    dy_dx_perp = -1 / dy_dx
    return dy_dx_perp, x0, y0

def find_intersection(dy_dx1, x0_1, y0_1, dy_dx2, x0_2, y0_2):
    A = np.array([[-dy_dx1, 1], [-dy_dx2, 1]])
    B = np.array([y0_1 - dy_dx1 * x0_1, y0_2 - dy_dx2 * x0_2])
    intersection_point = np.linalg.solve(A, B)
    return intersection_point

def distance(point1, point2):
    return np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

m = 0.0
p = 0.0
t = 0.12
c = 1.0
cutoff_ratio = 0.1

x_upper, y_upper, x_lower, y_lower = naca_4_digit_airfoil(m, p, t, c)
x_upper_cutoff, y_upper_cutoff, x_lower_cutoff, y_lower_cutoff = cutoff_trailing_edge_symmetrical(x_upper, y_upper, x_lower, y_lower, cutoff_ratio, c)

dy_dx_upper, x0_upper, y0_upper = tangent_line_at_cutoff(x_upper_cutoff, y_upper_cutoff)
dy_dx_lower, x0_lower, y0_lower = tangent_line_at_cutoff(x_lower_cutoff, y_lower_cutoff)

dy_dx_perp_upper, x0_perp_upper, y0_perp_upper = perpendicular_line_at_cutoff(dy_dx_upper, x0_upper, y0_upper)
dy_dx_perp_lower, x0_perp_lower, y0_perp_lower = perpendicular_line_at_cutoff(dy_dx_lower, x0_lower, y0_lower)

intersection = find_intersection(dy_dx_perp_upper, x0_perp_upper, y0_perp_upper, dy_dx_perp_lower, x0_perp_lower, y0_perp_lower)

tangent_upper_eq = f"y = {dy_dx_upper:.4f} * (x - {[x_upper_cutoff[-1]]}) + {[y_upper_cutoff[-1]]}"
tangent_lower_eq = f"y = {dy_dx_lower:.4f} * (x - {[x_lower_cutoff[-1]]}) + {[y_lower_cutoff[-1]]}"
perpendicular_upper_eq = f"y = {dy_dx_perp_upper:.4f} * (x - {[x_upper_cutoff[-1]]}) + {[y_upper_cutoff[-1]]}"
perpendicular_lower_eq = f"y = {dy_dx_perp_lower:.4f} * (x - {[x_lower_cutoff[-1]]}) + {[y_lower_cutoff[-1]]}"
top_cutoff_point = (x_upper_cutoff[-1], y_upper_cutoff[-1])
distance_to_top_cutoff = distance(intersection, top_cutoff_point)

arcsine_radians_top = math.asin(dy_dx_upper)
arcsine_radians_bot = math.asin(dy_dx_lower)
theta_top = (math.pi/2) + arcsine_radians_top
theta_bot = -(math.pi/2) + arcsine_radians_bot
arcsine_degrees_top2 = math.degrees(theta_top)
arcsine_degrees_bot2 = math.degrees(theta_bot)

print(f"Arcsine of top line in degrees: {arcsine_degrees_top2}")
print(f"Arcsine of bottom line in degrees: {arcsine_degrees_bot2}")
print("Equation of the tangent line at the upper surface cutoff point:", tangent_upper_eq)
print("Equation of the tangent line at the lower surface cutoff point:", tangent_lower_eq)
print("Equation of the perpendicular line at the upper surface cutoff point:", perpendicular_upper_eq)
print("Equation of the perpendicular line at the lower surface cutoff point:", perpendicular_lower_eq)
print(f"Intersection point of the two perpendicular lines: {intersection}")
print(f"Distance from the intersection point to the top of the cutoff: {distance_to_top_cutoff}")

plt.figure(figsize=(10, 5))
plt.plot(np.concatenate([x_upper_cutoff, x_lower_cutoff[::-1]]), np.concatenate([y_upper_cutoff, y_lower_cutoff[::-1]]), label='Modified NACA 0012 Airfoil (10% Cutoff)', linestyle='--')

arc = Arc(intersection, 2*distance_to_top_cutoff, 2*distance_to_top_cutoff, theta2=arcsine_degrees_top2, theta1=arcsine_degrees_bot2, edgecolor='r')
plt.gca().add_patch(arc)

# Plot tangent lines at cutoff points
plt.plot([x_upper_cutoff[-1], x_upper_cutoff[-1] + 0.001], [y_upper_cutoff[-1], y_upper_cutoff[-1] + 0.001 * dy_dx_upper], 'b-', label='Upper Tangent Line')
plt.plot([x_lower_cutoff[-1], x_lower_cutoff[-1] + 0.001], [y_lower_cutoff[-1], y_lower_cutoff[-1] + 0.001 * dy_dx_lower], 'g-', label='Lower Tangent Line')

# Plot intersection point
plt.plot(intersection[0], intersection[1], 'ro', label='Intersection Point')

plt.xlabel('x')
plt.ylabel('y')
plt.title('NACA 0012 Airfoil with 10% Trailing Edge Cutoff and Tangent Lines')
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()
