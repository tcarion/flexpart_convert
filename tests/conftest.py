import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
import xarray as xr
import numpy as np

@pytest.fixture
def flexpart_nc_path():
    """Path to the FLEXPART NetCDF test file."""
    return Path("tests/test_data/grid_conc_20210905000000.nc")

@pytest.fixture
def flexpart_dataset(flexpart_nc_path):
    """Load the FLEXPART dataset."""
    return xr.open_dataset(flexpart_nc_path)