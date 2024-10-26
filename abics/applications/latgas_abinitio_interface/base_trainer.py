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

from typing import Sequence, Type
import os

class TrainerBase(object):

    def __init__(
        self,
        structures: Sequence,
        energies: Sequence[float],
        generate_inputdir: os.PathLike,
        train_inputdir: os.PathLike,
        predict_inputdir: os.PathLike,
        generate_exe: str,
        train_exe: str,
        ):
        ...
        
    def prepare(self, latgas_mode = True, st_dir = ""):
        ...

    def generate_run(self, xsfdir="", generate_dir="generate"):
        """ generate training dataset for specific trainer

        Args:
            xsfdir (str, optional): _description_. Defaults to "".
            generate_dir (str, optional): . Defaults to "generate".
        """
        ...

    def generate_wait(self):
        """ wait for generate_run to finish
        """
        ...

    def train(self, train_dir = "train"):
        ...

    def new_baseinput(self, baseinput_dir, train_dir = "train"):
        """generate new baseinput directory/files for prediction

        Args:
            baseinput_dir (str): new baseinput directory
            train_dir (str, optional): directory including training result. Defaults to "train".
        """
        ...


__trainer_table = {}

def register_trainer(trainer_name: str, trainer_class: str, trainer_module: str) -> None:
    """
    Register trainer class.

    Parameters
    ----------
    trainer_name : str
        trainer name (case insensible).
    trainer_class : str
        trainer class, which should be a subclass of trainerBase.
    trainer_module : str
        Module name including the trainer class.
    """

    __trainer_table[trainer_name.lower()] = (trainer_class, trainer_module)


def get_trainer_class(trainer_name) -> Type[TrainerBase]:
    """
    Create trainer instance.

    Parameters
    ----------
    trainer_name : str
        trainer name (case insensible).
    params : ALParams or DFTParams
        Parameters.

    Returns
    -------
    trainer : TrainerBase
        trainer instance.
    """
    sn = trainer_name.lower()
    if sn not in __trainer_table:
        raise ValueError(f"Unknown trainer: {trainer_name}")

    import importlib
    trainer_class_name, trainer_module = __trainer_table[sn]
    mod = importlib.import_module(trainer_module)
    trainer_class = getattr(mod, trainer_class_name)
    if TrainerBase not in trainer_class.mro():
        raise TypeError("trainer_class must be a subclass of TrainerBase")
    
    return trainer_class