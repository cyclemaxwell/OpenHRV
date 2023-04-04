from typing import Tuple
from scipy.interpolate import CubicSpline
import numpy as np
import time
from PySide6.QtCore import QObject


class Pacer(QObject):
    def __init__(self):
        super().__init__()

        theta = np.linspace(0, 2 * np.pi, 40)
        self.cos_theta = np.cos(theta)
        self.sin_theta = np.sin(theta)
        self._interpolator = self._create_interpolator(
            2, straight_fraction=0.6
        )  # need to connect these to gui later

    def breathing_pattern(self, breathing_rate, time):
        """Returns radius of pacer disk.

        Radius is modulated according to sinusoidal breathing pattern
        and scaled between 0 and 1.
        """
        # return 0.5 + 0.5 * np.sin(2 * np.pi * breathing_rate / 60 * time)
        return self._interpolator(2 * np.pi * breathing_rate / 60 * time)

    def update(self, breathing_rate):
        """Update radius of pacer disc.

        Make current disk radius a function of real time (i.e., don't
        precompute radii with fixed time interval) in order to compensate for
        jitter or delay in QTimer calls.
        """
        radius = self.breathing_pattern(breathing_rate, time.time())
        x = radius * self.cos_theta
        y = radius * self.sin_theta
        return (x, y)

    def _create_points(
        self, asymmetry_parameter: float, straight_fraction: float = 0.8
    ) -> Tuple[list, list]:
        """Returns a list of points to be used to generate an interpolating function for a single cycle of an asymmetric wave

        Args:
            asymmetry_parameter (float): ratio of downward length to upward
            straight_fraction (float, optional): Fraction of wave that is a straight line. Defaults to 0.8.

        Returns:
            Tuple[list,list]: Two lists, the first time and the second a y value for a wave
        """
        non_straight = (1 - straight_fraction) / 2
        t0, y0 = 0, 1
        t6, y6 = 1, 1
        t3, y3 = (asymmetry_parameter / (asymmetry_parameter + 1)) * t6, 0
        dy = 0.25 * non_straight * ((y0 - y3) / (t3 - t0))
        t1, y1 = non_straight * (t3 - t0), y0 - dy
        t2, y2 = (1 - non_straight) * (t3 - t0), y3 + dy
        t4, y4 = non_straight * (t6 - t3) + t3, y3 + dy
        t5, y5 = (1 - non_straight) * (t6 - t3) + t3, y6 - dy
        return [t0, t1, t2, t3, t4, t5, t6], [y0, y1, y2, y3, y4, y5, y6]

    def _create_interpolator(
        self, asymmetry_parameter: float, straight_fraction: float = 0.8
    ) -> CubicSpline:
        """Returns a periodic cubic spline interpolating function with zero slope at
        the ends and dips to zero at with a ratio of asymmetry_parameter on the left side to the right side.
        The period is 2 pi.

        Args:
            asymmetry_parameter (float): ratio of downward length to upward
            straight_fraction (float, optional): The shorter the straight fraction, the flatter the curves are at the extrema. Small values get weird. Defaults to 0.8.

        Returns:
            scipy.interpolate._cubic.CubicSpline: Interpolating spline
        """
        x, y = self._create_points(
            asymmetry_parameter, straight_fraction=straight_fraction
        )
        return lambda x1: CubicSpline(x, y, bc_type="clamped")((x1 / (2 * np.pi)) % 1.0)
