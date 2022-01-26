"""
Script for defining a Particle3D object for use in N-body simulations
Author: James Horlock
Version: 26/01/2022
"""

import numpy as np



class Particle3D(object):

    def __init__(self, pos, vel, mass, label):
        """
        Initialise a Particle3D instance

        :param pos : position as numpy array
        :param vel : velocity as numpy array 
        :param mass: mass as float
        
        """
        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label = label
        

    def __str__ (self):
        """
        Define output format
        for p=([1.0,2.0,3.0],[4.0,5.0,6.0],1.0)
        this will print as 
        "position = ([1,2,3]),velocity= ([4,5,6]) m=1.0
        """
        return "position= " + str(self.position) + ",v = " + str(self.velocity) + ",m = "+ str(self.mass)

    def kinetic_energy (self):
        """
        return kinetic energy as 
        1/2*mass*velocity
        """
        return 0.5*self.mass*np.linalgnorm(self.velocity)**2

    def leap_velocity(self, dt, force):
        """
        First-order velocity update 
        v(t+dt)= v(t) +dt*F(t)    
        :param dt : timestep as float
        :param force : force on particle as numpy array

        """

        self.velocity = self.velocity + dt * force / self.mass

    def leap_pos_first(self, dt):
        """
        first-order position update 
        x(t+dt)=x(t) + dt*v(t)
        param dt : timestep as float
            
        """
        self.position = self.position + dt * self.velocity

    def leap_pos_second(self, dt, force):
        """
        second-order position update
        x(t+dt)=x(t) + dt*v(t) + 0.5 * dt^2 * F(t)
        :param dt: timestep as float
        :param force : force on particle as numpy array

        """
        first_order_step = dt * self.velocity
        second_order_step = 0.5 * force * dt**2 / self.mass
        first_plus_second = first_order_step + second_order_step
        self.position = self.position + first_plus_second

                            
           

    



            
                    
                    

                
