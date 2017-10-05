# This file is part of eddylicious
# (c) Timofey Mukha
# The code is released under the GNU GPL Version 3 licence.
# See LICENCE.txt and the Legal section in the User Guide for more information

from eddylicious.generators.helper_functions import blending_function
from eddylicious.generators.helper_functions import chunks_and_offsets
import numpy as np
import pytest


# Test the value of W for where it is known
def test_blending_function_values_with_default_parameters():

    # 0 at eta = 0
    assert blending_function(np.array([0])) == 0

    # 1 at eta = 1
    assert blending_function(np.array([1])) == 1

    # 0.5 at eta = b
    assert blending_function(np.array([0.2])) == 0.5

    # 1 at eta > 1
    assert blending_function(np.array([2])) == 1


# Check that having a non-positive number of procs raises an error
def test_chunks_and_offsets_nonpositive_processor_number():
    with pytest.raises(Exception):
        res = chunks_and_offsets(0, 5)


# Check that there should be data for every proc
def test_chunks_and_offsets_data_smaller_than_n_procs():
    with pytest.raises(Exception):
        res = chunks_and_offsets(7, 5)


# Check situation when nProcs = size
def test_chunks_and_offsets_n_procs_equals_data():
        [chunks, offsets] = chunks_and_offsets(10, 10)
        assert not np.any(chunks - np.ones(10))
        assert not np.any(offsets - np.arange(10, dtype=np.int64))

