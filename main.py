"""
A simple simulation of the solar system that uses astroquery Horizons data.

Author: James Horlock
"""

import numpy as np
import os
from Particle3D import Particle3D
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
from astropy.time import Time

def calculate_gravatational_force(planet_1, planet_2, g, relativistic):
    """
    Calculates the gravatational force on planet_1 due to planet_2
    :param planet_1: Particle3D object of a planet
    :param planet_2: Particle3D object of a planet
    :param g: Gravitational constant converted to units of solar masses, AU and days as floats
    :param relativistic: determines whether or not should consider relativistic component of grav force
    :return: force_vec the gravitational force on planet 1 due to planet 2 as numpy array
    """
    if not relativistic:
        sep = planet_1.position - planet_2.position
        force = - g * (planet_1.mass * planet_2.mass) / ((np.linalg.norm(sep) + 0.01) ** 3)
        force_vec = force * sep
    else:
        raise Exception("Relativistic calculation not yet supported. Incorrect configuration detected.")
    return force_vec

def update_planet_coords(planet_coords, planets, planet_indices):
    """
    Updates the coordinates dictionary for each planet with its current position.
    :param planet_coords: Dictionary mapping planet names to lists of their coordinates
    :param planets: List of Particle3D planet objects
    :param planet_indices: Dictionary mapping planet names to their index in planets list
    :return: Updated planet_coords dictionary
    """
    for key in planet_coords:
        planet_coords[key].append(np.array([planets[planet_indices[key]].position[0], planets[planet_indices[key]].position[1]]))
    return planet_coords

def get_key_from_value(index, planet_indices):
    """
    Retrieves the planet name corresponding to a given index from the planet_indices dictionary.
    :param index: Integer index of the planet
    :param planet_indices: Dictionary mapping planet names to their indices
    :return: Planet name as a string
    :raises ValueError: If no key matches the given index
    """
    for key, value in planet_indices.items():
        if value == index:
            return key
    raise ValueError(f"No key found for value {index} in planet_indices.")

def show_progress(current_step, total_steps):
    """
    Clears the terminal and prints the current progress percentage of the simulation.
    :param current_step: Current simulation step (int)
    :param total_steps: Total number of simulation steps (int)
    """
    if current_step % 1000 != 0:
        return

    if os.name == 'nt':
        os.system('cls')
    else:
        os.system('clear')

    print(str(round(current_step * 100 / total_steps, 1)) + " %")

def plot_planets(planet_coords, planet_name):
    """
    Plots the orbits of planets using their coordinates.
    :param planet_coords: Dictionary mapping planet names to lists of their coordinates
    :param planet_name: Name of the planet to set axis limits for (str)
    """
    orbit_x = list(map((lambda x: x[0]), planet_coords[planet_name]))
    orbit_y = list(map((lambda x: x[1]), planet_coords[planet_name]))
    max_x = max(orbit_x) + 1
    max_y = max(orbit_y) + 1
    min_x = min(orbit_x) - 1
    min_y = min(orbit_y) - 1
    plt.style.use('dark_background')
    plt.xlim(min_x, max_x)
    plt.ylim(min_y, max_y)
    for key in planet_coords:
        x_coords = [coord[0] for coord in planet_coords[key]]
        y_coords = [coord[1] for coord in planet_coords[key]]
        plt.plot(x_coords, y_coords, label=key)
    plt.legend()
    plt.title("Orbits of the Solar System")
    plt.show()

def main():
    # Initialises constants and variables
    g = 2.95 * 10 ** - 4 # Gravitational constant in units of solar masses, AU and days

    planet_indices = { # Dictionary mapping planets order to 0 index
        "Sun": 0,
        "Mercury": 1,
        "Venus": 2,
        "Earth": 3,
        "Mars": 4,
        "Jupiter": 5,
        "Saturn": 6,
        "Uranus": 7,
        "Neptune": 8
    }

    inner_planet_coords = {
        "Mercury" : [],
        "Venus" : [],
        "Earth" : [],
        "Mars" : [],
    }
    outer_planet_coords = {
        "Jupiter": [],
        "Saturn": [],
        "Uranus": [],
        "Neptune": []
    }

    masses = {
        "Sun": 1.0,
        "Mercury": 0.166 * 10 ** -6,
        "Venus": 2.447 * 10 ** -6,
        "Earth": 3.04 * 10 ** -6,
        "Mars": 3.23 * 10 ** -7,
        "Jupiter": 954 * 10 ** -6,
        "Saturn": 285 * 10 ** -6,
        "Uranus": 43.6 * 10 ** -6,
        "Neptune": 51.2 * 10 ** -6
    }

    planets = []
    # Initial conditions for simulation
    planets.append(Particle3D(np.array([0, 0, 0]), np.array([0, 0, 0]), masses["Sun"], 0))
    # Desired number of time steps
    total_steps = int(165 * 365 / 0.01) # One neptune year 
    # Inputs the positions and velocities of the 8 planets on the 1st of January 2019
    for i in range(8):
        planet = Horizons(id=i+1, epochs=Time("2019-01-01").jd).vectors()
        planet_name = get_key_from_value(i+1, planet_indices)
        mass = masses[planet_name]
        x = np.float64(planet['x'])
        y = np.float64(planet['y'])
        z = np.float64(planet['z'])
        vx = np.float64(planet['vx'])
        vy = np.float64(planet['vy'])
        vz = np.float64(planet['vz'])
        planets.append(Particle3D(np.array([x, y, z]), np.array([vx, vy, vz]), mass, i+1))

    inner_planet_coords = update_planet_coords(inner_planet_coords, planets, planet_indices)
    outer_planet_coords = update_planet_coords(outer_planet_coords, planets, planet_indices)

    # Begins simulation loop
    for i in range(total_steps):
        # Creates list for each time steps force calculation
        forces = []
        
        show_progress(i, total_steps)
            
        # Calculates the force on each planet due to the Sun
        force_array = np.zeros([len(planets), len(planets), 3])
        for j in range(len(planets)):
            for k in range(j + 1, len(planets)):
                force = calculate_gravatational_force(planets[j], planets[k], g, False)
                force_array[j, k, :] = force
                force_array[k, j, :] = - force
            # Combines the forces on a single planet into one vector
            forces.append(sum(force_array[j, :, :]))
            # Performs a second-order position and velocity update
            planets[j].leap_pos_second(0.01, forces[j])
            planets[j].leap_velocity(0.01, forces[j])

        if i % 100 == 0:
            inner_planet_coords = update_planet_coords(inner_planet_coords, planets, planet_indices)
            outer_planet_coords = update_planet_coords(outer_planet_coords, planets, planet_indices)
    plot_planets(inner_planet_coords, "Mars") # Uses Mars to set axis limits
    plot_planets(outer_planet_coords, "Neptune") # Uses Neptune to set axis limits

if __name__ == "__main__":
    main()