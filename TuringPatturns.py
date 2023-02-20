# import numba as nb
import math
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
from mpl_toolkits import mplot3d

def ToroidalDomain(r, R, N):
    """
    Creates the surface of a toroidal domain.
    r: The radial distance from the center of the torus to the center of the tube.
    R: The radial distance of the center of the tube of the surface of the torus.
    N: The number of points in each axial domain.
    """
    U = np.linspace(0, 2*np.pi, N)
    V = np.linspace(0, 2*np.pi, N)
    U, V = np.meshgrid(U, V)
    X = (r+R*np.cos(V))*np.cos(U)
    Y = (r+R*np.cos(V))*np.sin(U)
    Z = r*np.sin(V)
    return X, Y, Z


def HelicalDomain(h, R, a, L, N):
    """
    Creates a surface of a helical domain.
    h:
    R:
    a:
    L: The height of the coiled helix.
    N: The number of points in each plotting axes.
    """
    u = np.linspace(0, 2*math.pi, N)
    t = np.linspace(0, L, N)
    u, t = np.meshgrid(u, t)
    X = h*t + R*a*np.sin(u)/np.sqrt(R**2 + h**2)
    Y = (R*np.cos(t) - a*np.cos(t)*np.cos(u)+(h*a*np.sin(t)*np.sin(u))/np.sqrt(R**2 + h**2))
    Z = R*np.sin(t) - a*np.sin(t)*np.cos(u) - (h*a*np.cos(t)*np.sin(u))/np.sqrt(R**2 + h**2)
    return X, Y, Z

if __name__ == "__main__":
    # Change matplotlib plotting font to LaTeX font for consistency in the Capstone document. 
    #plt.rc('text', usetex = True); plt.rc('font', family = 'serif')

    # Toroidal Domain:
    x, y, z = ToroidalDomain(r = 2, R = 1, N = 300)
    fig = plt.figure(figsize = (12, 8)) # Creating figure
    ax = plt.axes(projection = '3d')
    ax.plot_surface(x, y, z, antialiased = True, color = "Gray") # Creating plot
    ax.set_box_aspect([1,1,1])
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$z$")
    plt.show() # show plot

    # Helical Domain:
    x, y, z = HelicalDomain(h = 1, R = 3, a = 0.4, L = 10, N = 300)
    fig = plt.figure(figsize = (12, 8)) # Creating figure
    ax = plt.axes(projection = '3d')
    ax.plot_surface(x, y, z, antialiased = True, color = "Grey") # Creating plot
    ax.set_box_aspect([1,1,1])
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$z$")
    plt.show() # show plot
