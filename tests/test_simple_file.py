import simple_file
import numpy as np


def test_add():
    a = np.array([1, 2, 3])
    b = simple_file.add(a, 1)
    assert (b == np.array([2, 3, 4])).all()
