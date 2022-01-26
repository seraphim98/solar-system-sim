"""
A simple simulation of the solar system that uses astroquery Horizons data. This only considers
gravitational forces due to the sun and assume all other interaction are negligible.
Author: James Horlock
Version: 26/01/2022
"""

import numpy as np
from Particle3D import Particle3D
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
from astropy.time import Time

def force_calc(planet_1, planet_2, G):

    """
    Calculates the gravatational force on planet_1 due to planet_2
    :param planet_1: Particle3D object of a planet
    :param planet_2: Particle3D object of a planet
    :param G: Gravitational constant converted to units of solar masses, AU and days as floats
    :return: force_vec the gravitational force on planet 1 due to planet 2 as numpy array
    """

    sep = planet_1.position - planet_2.position
    force = - G * (planet_1.mass * planet_2.mass) / (np.linalg.norm(sep) ** 3)
    force_vec = force * sep

    return force_vec




def main():

    planets = []
    G = 2.95 * 10 ** - 4
    #initial conditions for simulation
    planets.append(Particle3D(np.array([0, 0, 0]), np.array([0, 0, 0]), 1, 0))
    masses = [1, 0.166 * 10 ** - 6, 2.447*10**-6, 3.04 * 10 ** -6, 3.23 * 10 ** -7, 954*10**-6,  285*10**-6,
              43.6*10**-6, 51.5*10**-6]
    # desired number of time steps
    numstep = int(3 * 365 / 0.002)
    # inputs the positions and velocities of the 8 planets on the 1st of January 2019
    for i in range(8):
        planet = Horizons(id=i+1, epochs=Time("2019-01-01").jd, location="@sun").vectors()
        x = np.float64(planet['x'])
        y = np.float64(planet['y'])
        z = np.float64(planet['z'])
        vx = np.float64(planet['vx'])
        vy = np.float64(planet['vy'])
        vz = np.float64(planet['vz'])
        planets.append(Particle3D(np.array([x, y, z]), np.array([vx, vy, vz]), masses[i+1], i+1))
    # creates empty lists for use of visulisation
    mercury_pos, venus_pos, earth_pos, mars_pos, mars_x, mars_y, sun_pos = [], [], [], [], [], [], []
    mercury_pos_y, venus_pos_y, earth_pos_y, sun_pos_y = [], [], [], []
    # begins simualtion loop
    for i in range(numstep):

        # creates list for each time steps force calculation
        forces = [0]
        # shows progress
        if i % 1000 == 0:
            print(str(round(i*100/numstep, 1)) + " %")
        # calculates the force on each planet due to the Sun
        for j in range(len(planets)):
            if j != 0:
                forces.append(force_calc(planets[j], planets[0], G))
        # performs a second-order position and velocity update
        for j in range(len(planets)):
            planets[j].leap_pos_second(0.002, forces[j])
            planets[j].leap_velocity(0.002, forces[j])
            # stores positions of the sun and first four planets
            if j == 1 and i % 100 == 0:
                mercury_pos.append(planets[j].position[0])
                mercury_pos_y.append(planets[j].position[1])
            if j == 2 and i % 100 == 0:
                venus_pos.append(planets[j].position[0])
                venus_pos_y.append(planets[j].position[1])
            if j == 3 and i % 100 == 0:
                earth_pos.append(planets[j].position[0])
                earth_pos_y.append(planets[j].position[1])
            if j == 4 and i % 100 == 0:
                mars_pos.append(planets[j].position)
                mars_x.append(planets[j].position[0])
                mars_y.append(planets[j].position[1])
                sun_pos.append(0)
                sun_pos_y.append(0)
    #calculates the plot limits from Mars' orbit
    max_x = max(mars_x) + 1
    max_y = max(mars_y) + 1
    min_x = min(mars_x) - 1
    min_y = min(mars_y) - 1
    # plots the orbits of the planets
    plt.style.use('dark_background')
    plt.xlim(min_x, max_x)
    plt.ylim(min_y, max_y)
    plt.plot(sun_pos, sun_pos_y, 'y', linewidth=4)
    plt.plot(mercury_pos, mercury_pos_y, 'brown')
    plt.plot(venus_pos, venus_pos_y, 'm')
    plt.plot(earth_pos,earth_pos_y, 'blue')
    plt.plot(mars_x, mars_y, 'red')
    plt.show()


main()



