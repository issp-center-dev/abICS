# ab-Initio Configuration Sampling tool kit (abICS)
# Copyright (C) 2019- The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

from __future__ import annotations

from typing import Any

from abc import ABCMeta, abstractmethod

import sys

from abics import __version__

"""Defines base classes for Monte Carlo simulations"""


class Model(metaclass=ABCMeta):
    """This class defines a model whose energy equals 0 no matter the configuration, and the configuration
    never changes.
    This is a base template for building useful models."""

    model_name = ""

    # def __init__(self):

    @abstractmethod
    def energy(self, config) -> float:
        """
        Calculate energy of configuration: input: config

        Parameters
        ----------
        config: config object
            configuration

        Returns
        -------
        energy: float

        """
        ...

    @abstractmethod
    def trialstep(self, config, energy: float) -> tuple[Any, float]:
        """Define a trial step on config

        Returns dconfig, which can contain the minimal information for
        constructing the trial configuration from config to be used in newconfig().
        Make sure that config is unchanged.

        Parameters
        ----------
        config: config object
            current configuration

        energy: float
            current energy

        Returns
        -------
        dconfig: config object
            The minimal information for constructing the trial configuration
            from config to be used in newconfig()
        dE: float
            Energy difference
        """

        # Return only change in configuration dconfig so that
        # you don't have to copy entire configurations,
        # which can sometimes be costly
        ...

    @abstractmethod
    def newconfig(self, config, dconfig):
        """
        Update config by using the trial step, dconfig

        Parameters
        ----------
        config: config object
            Original configuration
            This may be mutated through this function
        dconfig: config object
            Difference of configuration

        Returns
        -------
        config: config object
            updated configuration
        """
        return config

# For backward compatibility
if __version__ < "3":
    model = Model
