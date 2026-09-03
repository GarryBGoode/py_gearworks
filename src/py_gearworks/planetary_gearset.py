# Copyright 2026 Gergely Bencsik
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from dataclasses import dataclass
from enum import Enum
import numpy as np

# planetary base equation (willis)
# w_ring * n_ring = w_carry * (n_ring + n_sun) - w_sun * n_sun


@dataclass
class PlanetaryRotations:
    """Angular velocities for all members of a planetary gearset.

    Attributes
    ----------
    ring, carry, sun, planet : float
        Angular velocity of the corresponding gearset member. Units are inherited
        from the input speed used to calculate the result.
    """

    ring: float
    carry: float
    sun: float
    planet: float


class PlanetaryGearType(Enum):
    """Selectable members of a planetary gearset."""

    SUN = 1
    PLANET = 2
    RING = 3
    CARRY = 4


class PlanetaryGearset:
    """
    A helper class for managing planetary gearsets.
    This class handles transmission ratios, possible planet locations, and other basic
    data related to planetary gears.
    This class is independent of gear types or tooth profiles, only considers number of
    teeth.
    """

    def __init__(self, n_sun, n_ring, n_planet, num_planets):
        """Initialize a planetary gearset.

        Parameters
        ----------
        n_sun, n_ring, n_planet : int
            Tooth counts of the sun, ring, and each planet gear.
        num_planets : int
            Number of planet gears in the assembly.
        """
        self.n_sun = n_sun
        self.n_ring = n_ring
        self.n_planet = n_planet
        self.num_planets = num_planets

    def get_planet_position_angles(self):
        """Return tooth-aligned candidate planet positions in radians.

        Returns
        -------
        np.ndarray
            Equally spaced angles for the ``n_sun + n_ring`` valid tooth phases.
        """
        # w_ring * n_ring = w_carry * (n_ring + n_sun) - w_sun * n_sun
        # 1 tooth up for ring, sun is stationary
        # w_ring = 2pi=   w_carry * (n_ring + n_sun)
        # w_carry = 2pi / (n_ring + n_sun)

        angles = (
            np.linspace(0, 1, (self.n_sun + self.n_ring), endpoint=False) * 2 * np.pi
        )
        return angles

    def min_sun_ratio_from_planets(self):
        """Return the minimum sun-to-ring pitch-radius ratio for planet clearance.

        The result models adjacent planet pitch circles as touching. It is a
        geometric lower bound; addenda and tooth clearance require additional
        design margin.

        Returns
        -------
        float
            Minimum allowable ratio of sun pitch radius to ring pitch radius.
        """

        R_poly = np.sin(np.pi / self.num_planets)
        R_sun = (1 - R_poly) / (1 + R_poly)
        return R_sun

    def max_num_planets(self, addendum_ratio=1.0):
        """Return the maximum evenly spaced planets that fit within the gearset.

        Parameters
        ----------
        addendum_ratio : float, default=1.0
            Addendum height relative to the module used for planet clearance.

        Returns
        -------
        int
            Largest whole number of evenly distributed planet gears that fit.
        """
        R_poly = (self.n_planet * 2) / (self.n_ring + self.n_sun)
        R_poly *= (self.n_planet + addendum_ratio * 2) / self.n_planet
        num_planets = np.pi / np.arcsin(R_poly)
        return int(np.floor(num_planets))

    def get_planet_angles_distributed(self):
        """
        Returns the angles of the planets distributed (close to) evenly around the sun
        gear, given the number of planets. If even distribution is not possible, the
        closest possible valid distribution is returned. The angles are in radians.
        """
        angles_normalized = self.get_planet_position_angles() / np.pi / 2
        pos_ref = np.linspace(0, 1, self.num_planets, endpoint=False)
        increment = angles_normalized[1] - angles_normalized[0]

        angles = np.round(pos_ref / increment) * increment
        return angles * np.pi * 2

    def __willis_equation(self, w_ring, w_carry, w_sun):
        """Willis equation rearranged to zero. Returns zero if the given speeds satisfy
        the equation, otherwise returns the difference."""
        return (
            -w_ring * self.n_ring
            + w_carry * (self.n_ring + self.n_sun)
            - w_sun * self.n_sun,
        )

    def get_ratio(
        self, w=1, driven=PlanetaryGearType.SUN, stationary=PlanetaryGearType.CARRY
    ):
        """Calculate member velocities from a driven and stationary gearset member.

        Parameters
        ----------
        w : float, default=1
            Angular velocity assigned to ``driven``.
        driven : PlanetaryGearType, default=PlanetaryGearType.SUN
            Gearset member supplied with the input velocity.
        stationary : PlanetaryGearType, default=PlanetaryGearType.CARRY
            Gearset member held at zero velocity. It must differ from ``driven``.

        Returns
        -------
        PlanetaryRotations
            Angular velocities of the sun, ring, carrier, and planet members.

        Raises
        ------
        ValueError
            If the driven/stationary combination is not supported.
        """
        if driven == PlanetaryGearType.SUN:
            if stationary == PlanetaryGearType.CARRY:
                w_sun = w
                w_carry = 0
                w_ring = (
                    w_carry * (self.n_ring + self.n_sun) - w_sun * self.n_sun
                ) / self.n_ring
            elif stationary == PlanetaryGearType.RING:
                w_sun = w
                w_ring = 0
                w_carry = (w_ring * self.n_ring + w_sun * self.n_sun) / (
                    self.n_ring + self.n_sun
                )
            else:
                raise ValueError("Invalid stationary gear type for driven SUN.")
        elif driven == PlanetaryGearType.RING:
            if stationary == PlanetaryGearType.CARRY:
                w_ring = w
                w_carry = 0
                w_sun = (
                    w_carry * (self.n_ring + self.n_sun) - w_ring * self.n_ring
                ) / (-self.n_sun)
            elif stationary == PlanetaryGearType.SUN:
                w_ring = w
                w_sun = 0
                w_carry = (w_ring * self.n_ring + w_sun * self.n_sun) / (
                    self.n_ring + self.n_sun
                )
            else:
                raise ValueError("Invalid stationary gear type for driven RING.")
        elif driven == PlanetaryGearType.CARRY:
            if stationary == PlanetaryGearType.SUN:
                w_carry = w
                w_sun = 0
                w_ring = (
                    w_carry * (self.n_ring + self.n_sun) - w_sun * self.n_sun
                ) / self.n_ring
            elif stationary == PlanetaryGearType.RING:
                w_carry = w
                w_ring = 0
                w_sun = (
                    w_carry * (self.n_ring + self.n_sun) - w_ring * self.n_ring
                ) / (-self.n_sun)
            else:
                raise ValueError("Invalid stationary gear type for driven CARRY.")
        else:
            raise ValueError("Invalid driven gear type.")

        # Calculate planet speed based on the other speeds and the Willis equation.
        # Using the fact that the planet speed is the average of the sun and ring speeds.
        # This is a simplification and may not hold in all cases, but it is a common assumption in planetary gear analysis.
        w_planet = (w_sun + w_ring) / 2

        return PlanetaryRotations(
            sun=w_sun, ring=w_ring, carry=w_carry, planet=w_planet
        )
