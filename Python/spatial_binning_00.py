"""
Copyright (C) 2022 Daniel Utt (utt@mm.tu-darmstadt.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

from ovito.data import *
from ovito.vis import *
import numpy as np


def create(frame: int, data: DataCollection, number_of_points=1000, width=20, seed=123):
    cell_size = 100
    cell_matrix = [
        [cell_size, 0, 0, 0],
        [0, cell_size, 0, 0],
        [0, 0, cell_size, 0],
    ]
    data.create_cell(cell_matrix, pbc=(True, True, True))

    rng = np.random.default_rng(seed)
    positions = np.zeros((number_of_points, 3))
    positions[:, 2] = cell_size * rng.random(size=number_of_points)
    positions[:, :2] = rng.normal(
        loc=cell_size / 2, scale=width, size=(number_of_points, 2)
    )

    particles = data.create_particles(count=len(positions))
    particles.create_property(
        "Position",
        data=positions,
    )

    type_property = particles.create_property("Particle Type")
    if len(type_property.types) == 0:
        type_property.types.append(
            ParticleType(id=1, name="Tealon", color=(1 / 10, 1, 1))
        )
    type_property[...] = [1] * len(positions)
