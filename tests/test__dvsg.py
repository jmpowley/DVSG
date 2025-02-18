import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import numpy as np
from dvsg import exclude_above_five_sigma, normalise_velocity_map, denormalise_velocity_map, calculate_DVSG

def test__normalise_velocity_map():
    velocity_map = np.array([[-3, -2, -1], [0, 0, 0], [1, 2, 3]])
    normalised_velocity_map = np.array([[-1, -2/3, -1/3], [0, 0, 0], [1/3, 2/3, 1]])

    assert np.allclose(normalise_velocity_map(velocity_map), normalised_velocity_map)

def test__denormalise_velocity_map():
    velocity_map = np.array([[-1, -2/3, -1/3], [0, 0, 0], [1/3, 2/3, 1]])

    unnormalised_velocity_map = np.array([[-3, -2, -1], [0, 0, 0], [1, 2, 3]])
    max_velocity = 3
    min_velocity = -3

    denormalised_velocity_map = denormalise_velocity_map(normalised_velocity_map=velocity_map, 
                                                         max_velocity=max_velocity, min_velocity=min_velocity)

    assert np.allclose(denormalised_velocity_map, unnormalised_velocity_map)