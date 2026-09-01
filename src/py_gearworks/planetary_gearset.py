import numpy as np

# planetary base equation (willis)
# w_ring * n_ring = w_carry * (n_ring + n_sun) - w_sun * n_sun


class PlanetaryGearset:
    def __init__(self, n_sun, n_ring, n_planet, num_planets):
        self.n_sun = n_sun
        self.n_ring = n_ring
        self.n_planet = n_planet
        self.num_planets = num_planets

    def get_planet_position_angles(self):
        """Returns all suitable angular positions (in radians) for the planets, given sun and ring gear sizes."""
        # w_ring * n_ring = w_carry * (n_ring + n_sun) - w_sun * n_sun
        # 1 tooth up for ring, sun is stationary
        # w_ring = 2pi=   w_carry * (n_ring + n_sun)
        # w_carry = 2pi / (n_ring + n_sun)

        angles = (
            np.linspace(0, 1, (self.n_sun + self.n_ring), endpoint=False) * 2 * np.pi
        )
        return angles

    def min_sun_ratio_from_planets(self):
        """
        Returns the minimum size of the sun gear relative to the ring gear, given the number of planets.
        This is derived from the fact that planets must fit next to each other. Geometrically modeled with touching circles,
        the minimum sun gear size for 1 and 2 planets is 0, for 3 planets it is 0.071, 4 planets give 0.172, and so on.
        Note that this only holds for pitch circles, actual sun gear size need to be larger and planets slightly smaller to fit.
        """

        R_poly = np.sin(np.pi / self.num_planets)
        R_sun = (1 - R_poly) / (1 + R_poly)
        return R_sun

    def max_num_planets(self):
        """Returns the maximum number of planets that can fit between the sun and ring gear, given their sizes."""
        R_poly = (1 - self.n_sun / self.n_ring) / (1 + self.n_sun / self.n_ring)
        num_planets = np.pi / np.arcsin(R_poly)
        return int(np.floor(num_planets))

    def get_planet_angles_distributed(self):
        """Returns the angles of the planets distributed evenly around the sun gear, given the number of planets."""
        angles_normalized = self.get_planet_position_angles() / np.pi / 2
        pos_ref = np.linspace(0, 1, self.num_planets, endpoint=False)
        increment = angles_normalized[1] - angles_normalized[0]

        angles = np.round(pos_ref / increment) * increment
        return angles * np.pi * 2


if __name__ == "__main__":
    test_gearset = PlanetaryGearset(n_sun=52, n_ring=92, n_planet=10, num_planets=3)

    print("planet positions:", test_gearset.get_planet_position_angles())

    print("min sun ratio:", test_gearset.min_sun_ratio_from_planets())

    print(
        "distributed planet angles:",
        test_gearset.get_planet_angles_distributed() / np.pi / 2,
    )

    print("max num planets:", test_gearset.max_num_planets())
