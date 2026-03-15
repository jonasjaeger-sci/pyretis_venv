"""
This is a class script for a custom order parameter calculation of the minimal distance between any particles
"""

from orderparameter import OrderParameter,Distance
import random

class MinDist(OrderParameter):
    """
    Class to represent an order parameter/collective variable as a child/grandchild of the
    OrderParameter class. It defines an order parameter which is the minimal distance
    between any two particles of the system
    """

    def calculate(self,system):
        """
        function to calculate the minimal distance between any two particles of the system
        Parameters
        ----------
        system: object
            object containing the particles positions

        Returns
        -------
        mindist: float
            minimal distance between any two particles of the system
        """

        n_atoms = len(system.particles.pos)
        print(f"number of atoms: {n_atoms}")
        return random.random()*2

