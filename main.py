"""
A simple simulation of the solar system that uses astroquery Horizons data. This only considers
gravitational forces due to the sun and assumes all other interaction are negligible.
Author: James Horlock
Version: 13/02/2024
"""

import numpy as np
import os
from Particle3D import Particle3D
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
from astropy.time import Time


def force_calc(planet_1, planet_2, g, relativistic):

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
        print("Not figured out yet")
        force_vec = 345678
        exit()
    return force_vec


def main():

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

    planets = []
    g = 2.95 * 10 ** - 4
    # initial conditions for simulation
    planets.append(Particle3D(np.array([0, 0, 0]), np.array([0, 0, 0]), 1, 0))
    masses = [1, 0.166 * 10 ** - 6, 2.447*10**-6, 3.04 * 10 ** -6, 3.23 * 10 ** -7, 954*10**-6,  285*10**-6,
              43.6*10**-6, 51.2*10**-6]
    # desired number of time steps
    numstep = int(165 * 365 / 0.01)
    # inputs the positions and velocities of the 8 planets on the 1st of January 2019
    for i in range(8):
        planet = Horizons(id=i+1, epochs=Time("2019-01-01").jd).vectors()
        x = np.float64(planet['x'])
        y = np.float64(planet['y'])
        z = np.float64(planet['z'])
        vx = np.float64(planet['vx'])
        vy = np.float64(planet['vy'])
        vz = np.float64(planet['vz'])
        planets.append(Particle3D(np.array([x, y, z]), np.array([vx, vy, vz]), masses[i+1], i+1))
    # creates empty lists for use of visualisation

    inner_planet_coords["Mercury"].append(np.array(planets[1].position[0], planets[1].position[1]))
    inner_planet_coords["Venus"].append(np.array(planets[2].position[0], planets[2].position[1]))
    inner_planet_coords["Earth"].append(np.array(planets[3].position[0], planets[3].position[1]))
    inner_planet_coords["Mars"].append(np.array(planets[4].position[0], planets[4].position[1]))
    outer_planet_coords["Jupiter"].append(np.array(planets[5].position[0], planets[5].position[1]))
    outer_planet_coords["Saturn"].append(np.array(planets[6].position[0], planets[6].position[1]))
    outer_planet_coords["Uranus"].append(np.array(planets[6].position[0], planets[6].position[1]))
    outer_planet_coords["Neptune"].append(np.array(planets[8].position[0], planets[8].position[1]))

    # begins simulation loop
    for i in range(numstep):

        # creates list for each time steps force calculation
        forces = []
        # shows progress
        if i % 1000 == 0:
            os.system('cls')
            print(str(round(i*100/numstep, 1)) + " %")
        # calculates the force on each planet due to the Sun
        for j in range(len(planets)):
            for k in range(j+1, len(planets)):
                force = force_calc(planets[j], planets[k], g, False)
                force_array[j, k, :] = force
                force_array[k, j, :] = - force
            # combines the forces on a single planet to one vecto
            forces.append(sum(force_array[j, :, :]))
            # performs a second-order position and velocity update
            planets[j].leap_pos_second(0.01, forces[j])
            planets[j].leap_velocity(0.01, forces[j])

        if i % 100 == 0:
            inner_planet_coords["Mercury"].append(np.array(planets[1].position[0], planets[1].position[1]))
            inner_planet_coords["Venus"].append(np.array(planets[2].position[0], planets[2].position[1]))
            inner_planet_coords["Earth"].append(np.array(planets[3].position[0], planets[3].position[1]))
            inner_planet_coords["Mars"].append(np.array(planets[4].position[0], planets[4].position[1]))
            outer_planet_coords["Jupiter"].append(np.array(planets[5].position[0], planets[5].position[1]))
            outer_planet_coords["Saturn"].append(np.array(planets[6].position[0], planets[6].position[1]))
            outer_planet_coords["Uranus"].append(np.array(planets[6].position[0], planets[6].position[1]))
            outer_planet_coords["Neptune"].append(np.array(planets[8].position[0], planets[8].position[1]))

    plot_planets(inner_planet_coords, "Mars")
    plot_planets(outer_planet_coords, "Neptune")

def plot_planets(planet_coords, largest_orbit):
    largest_orbit_x = map((lambda x: x[0]), planet_coords[largest_orbit])
    largest_orbit_y = map((lambda x: x[1]), planet_coords[largest_orbit])
    max_x = max(largest_orbit_x) + 1
    max_y = max(largest_orbit_y) + 1
    min_x = min(largest_orbit_x) - 1
    min_y = min(largest_orbit_y) - 1
    plt.style.use('dark_background')
    plt.xlim(min_x, max_x)
    plt.ylim(min_y, max_y)
    for key in planet_coords:
        plt.plot(planet_coords[key])
    plt.show()

main()