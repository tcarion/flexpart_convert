import pytest
from pathlib import Path
from flexpart_convert.converter import NetCDFConverter, NetCDFConversionError

def test_converter_init_flexpart(flexpart_nc_path):
    """Test initialization with FLEXPART NetCDF file."""
    converter = NetCDFConverter(flexpart_nc_path)
    assert converter.input_path == flexpart_nc_path
    assert "spec001_mr" in converter.variables
    assert "WD_spec001" in converter.variables
    assert "DD_spec001" in converter.variables
    assert "time" in converter.dimensions["spec001_mr"]
    converter.close()

def test_get_variable_info_spec001(flexpart_nc_path):
    """Test getting info for spec001_mr variable."""
    converter = NetCDFConverter(flexpart_nc_path)
    info = converter.get_variable_info("spec001_mr")
    
    assert info["dims"] == ("nageclass", "pointspec", "time", "height", "latitude", "longitude")
    assert info["shape"] == (1, 1, 2, 1, 16, 10)
    assert info["attrs"]["units"] == "ng kg-1"
    assert info["attrs"]["long_name"] == "AEROSOL"
    converter.close()

def test_convert_spec001_time_slices(flexpart_nc_path, tmp_path):
    """Test conversion of spec001_mr with time dimension."""
    converter = NetCDFConverter(flexpart_nc_path)
    results = converter.convert(
        variable="spec001_mr",
        output_format="tif",
        output_dir=tmp_path
    )
    
    # Should create 2 files (one per time step)
    assert len(results) == 2
    assert all("time_" in f.name for f in results)
    assert any("time_0" in f.name for f in results)
    assert any("time_1" in f.name for f in results)
    converter.close()

def test_convert_spec001_fixed_time(flexpart_nc_path, tmp_path):
    """Test conversion of spec001_mr with fixed time."""
    converter = NetCDFConverter(flexpart_nc_path)
    results = converter.convert(
        variable="spec001_mr",
        output_format="tif",
        output_dir=tmp_path,
        dim_selections={"time": 0, "nageclass": 0, "pointspec": 0, "height": 0}
    )
    
    # Should create 1 file (all non-spatial dimensions fixed)
    assert len(results) == 1
    assert "time_0" in results[0].name
    converter.close()

def test_convert_wd_spec001(flexpart_nc_path, tmp_path):
    """Test conversion of WD_spec001 variable."""
    converter = NetCDFConverter(flexpart_nc_path)
    results = converter.convert(
        variable="WD_spec001",
        output_format="tif",
        output_dir=tmp_path,
        dim_selections={"nageclass": 0, "pointspec": 0}
    )
    
    # Should create 2 files (one per time step)
    assert len(results) == 2
    assert all("WD_spec001" in f.name for f in results)
    converter.close()

def test_convert_with_height_dimension(flexpart_nc_path, tmp_path):
    """Test conversion with height dimension."""
    converter = NetCDFConverter(flexpart_nc_path)
    results = converter.convert(
        variable="spec001_mr",
        output_format="tif",
        output_dir=tmp_path,
        dim_selections={"time": 0, "nageclass": 0, "pointspec": 0}
    )
    
    # Should create 1 file (height dimension is size 1)
    assert len(results) == 1
    assert "height_0" in results[0].name
    converter.close()

def test_convert_invalid_variable(flexpart_nc_path):
    """Test conversion with invalid variable."""
    converter = NetCDFConverter(flexpart_nc_path)
    with pytest.raises(KeyError):
        converter.convert(
            variable="invalid_var",
            output_format="tif",
            output_dir="."
        )
    converter.close()

def test_convert_oro_variable(flexpart_nc_path, tmp_path):
    """Test conversion of ORO variable (2D field)."""
    converter = NetCDFConverter(flexpart_nc_path)
    results = converter.convert(
        variable="ORO",
        output_format="tif",
        output_dir=tmp_path
    )
    
    # Should create 1 file (no extra dimensions)
    assert len(results) == 1
    assert results[0].name == "grid_conc_20210905000000_ORO.tif"
    converter.close()