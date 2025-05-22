import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
from flexpart_convert.cli import main, _process_dim_selections

def test_process_dim_selections_flexpart():
    """Test processing FLEXPART dimension selections."""
    dim_args = [["time", "0"], ["nageclass", "0"]]
    mock_converter = MagicMock()
    mock_converter.dimensions = {"spec001_mr": ["time", "nageclass", "pointspec", "height", "latitude", "longitude"]}
    mock_converter.dataset = {
        "time": MagicMock(values=[0, 3600]),
        "nageclass": MagicMock(values=[0]),
        "pointspec": MagicMock(values=[0]),
        "height": MagicMock(values=[0])
    }
    
    result = _process_dim_selections(dim_args, "spec001_mr", mock_converter)
    assert result == {"time": 0, "nageclass": 0}

@patch("flexpart_convert.cli.NetCDFConverter")
def test_cli_main_flexpart_spec001(mock_converter, capsys, flexpart_nc_path):
    """Test CLI with FLEXPART spec001_mr variable."""
    mock_converter.return_value.convert.return_value = [
        Path("output1.tif"),
        Path("output2.tif")
    ]
    
    test_args = [
        str(flexpart_nc_path),
        "-v", "spec001_mr",
        "-f", "geotiff",
        "-d", "nageclass", "0",
        "-d", "pointspec", "0",
        "-d", "height", "0"
    ]
    
    with patch("sys.argv", ["flexpart_convert"] + test_args):
        main()
    
    captured = capsys.readouterr()
    assert "Conversion complete" in captured.out
    assert "output1.tif" in captured.out
    assert "output2.tif" in captured.out

@patch("flexpart_convert.cli.NetCDFConverter")
def test_cli_main_flexpart_oro(mock_converter, capsys, flexpart_nc_path):
    """Test CLI with FLEXPART ORO variable."""
    mock_converter.return_value.convert.return_value = [
        Path("oro_output.tif")
    ]
    
    test_args = [
        str(flexpart_nc_path),
        "-v", "ORO",
        "-f", "geotiff"
    ]
    
    with patch("sys.argv", ["flexpart_convert"] + test_args):
        main()
    
    captured = capsys.readouterr()
    assert "Conversion complete" in captured.out
    assert "oro_output.tif" in captured.out