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


def create(frame: int, data: DataCollection, edge_length=5):
    # edge_length = 5
    vertex_positions = (
        (0, 0, 0),
        (edge_length, 0, 0),
        (0, edge_length, 0),
        (edge_length, edge_length, 0),
        (0, 0, edge_length),
        (edge_length, 0, edge_length),
        (0, edge_length, edge_length),
        (edge_length, edge_length, edge_length),
    )

    particles = data.create_particles(count=len(vertex_positions))
    particles.create_property(
        "Position",
        data=vertex_positions,
    )

    type_property = particles.create_property("Particle Type")
    if len(type_property.types) == 0:
        type_property.types.append(
            ParticleType(id=1, name="Magentanium", color=(1, 0, 1))
        )
    type_property[...] = [1] * len(vertex_positions)

    edge_topology = (
        (0, 1),
        (0, 2),
        (1, 3),
        (2, 3),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
        (4, 5),
        (4, 6),
        (5, 7),
        (6, 7),
    )
    bonds = particles.create_bonds(count=len(edge_topology), vis_params={"width": 0.66})
    bonds.create_property("Topology", data=edge_topology)

    cell_matrix = [
        [edge_length + 2, 0, 0, -1],
        [0, edge_length + 2, 0, -1],
        [0, 0, edge_length + 2, -1],
    ]
    data.create_cell(cell_matrix, pbc=(False, False, False))
