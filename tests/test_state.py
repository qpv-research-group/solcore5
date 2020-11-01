import numpy as np
from solcore.state import State


def test_state():
    array1 = np.array([1, 2, 3])
    array2 = np.array([2, 3, 4])

    state1 = State()
    state1.a = 1
    state1.b = True
    state1.c = array1
    assert state1.__len__() == 3
    assert state1.a == 1
    assert state1.b == True
    assert all([input == output for input, output in zip(array1, state1.c)])

    state2 = State()
    state2.b = False
    state2.d = array2
    state1.update(state2)
    assert state1.__len__() == 4
    assert state1.a == 1
    assert state2.b == False
    assert all([input == output for input, output in zip(array1, state1.c)])
    assert all([input == output for input, output in zip(array2, state1.d)])
