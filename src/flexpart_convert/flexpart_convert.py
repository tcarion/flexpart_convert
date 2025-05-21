import argparse
import sys
import itertools
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import xarray as xr
import rioxarray as rio

class NetCDFConverter:
    """Convert NetCDF files to various geospatial formats."""
    
    def __init__(self, input_path: Union[str, Path]):
        """Initialize the converter with a NetCDF file.
        
        Args:
            input_path: Path to the NetCDF file
        """
        self.input_path = Path(input_path)
        if not self.input_path.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_path}")
            
        self.dataset = xr.open_dataset(self.input_path, decode_timedelta=True)
        # self.dataset = xr.open_dataset(self.input_path, engine="rasterio")
        self.variables = list(self.dataset.data_vars)
        self.dimensions = {var: list(self.dataset[var].dims) for var in self.variables}
        
    def get_variable_info(self, variable: str) -> Dict:
        """Get information about a variable's dimensions.
        
        Args:
            variable: Name of the variable to inspect
            
        Returns:
            Dictionary with dimension information
        """
        if variable not in self.variables:
            raise ValueError(f"Variable '{variable}' not found in dataset")
            
        var_data = self.dataset[variable]
        return {
            "dims": var_data.dims,
            "shape": var_data.shape,
            "attrs": dict(var_data.attrs)
        }
        
    def convert(
        self,
        variable: str = "spec001_mr",
        output_format: str = "tif",
        output_dir: Union[str, Path] = ".",
        dim_selections: Optional[Dict[str, Union[int, float, str]]] = None,
    ) -> List[Path]:
        """Convert the NetCDF variable to the specified format.
        
        Args:
            variable: Variable name to convert
            output_format: One of 'tif', 'geojson', or 'kml'
            output_dir: Directory to save output files
            dim_selections: Dictionary of dimension selections (e.g., {'time': 0})
            
        Returns:
            List of paths to created files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if output_format.lower() not in {"tif", "geojson", "kml"}:
            raise ValueError(f"Unsupported output format: {output_format}")
            
        var_data = self.dataset[variable]
        dim_selections = dim_selections or {}
        
        # Apply dimension selections
        for dim, value in dim_selections.items():
            if dim not in var_data.dims:
                raise ValueError(f"Dimension '{dim}' not found in variable '{variable}'")
            var_data = var_data.sel({dim: value})
            
        # Check remaining dimensions
        remaining_dims = [dim for dim in var_data.dims if dim not in {"lat", "latitude", "lon", "longitude", "x", "y"}]
        
        if len(remaining_dims) > 0:
            # Create a list of all dimension combinations we need to process
            dim_combinations = []
            
            # Get all possible values for each remaining dimension
            dim_values = {dim: var_data[dim].values for dim in remaining_dims}
            
            # Create a mesh grid of all combinations
            if len(remaining_dims) == 1:
                # Simple case - just one remaining dimension
                for val in dim_values[remaining_dims[0]]:
                    dim_combinations.append({remaining_dims[0]: val})
            else:
                # Multiple remaining dimensions - need all combinations
                value_combinations = itertools.product(*[dim_values[dim] for dim in remaining_dims])
                
                for combo in value_combinations:
                    dim_dict = {dim: val for dim, val in zip(remaining_dims, combo)}
                    dim_combinations.append(dim_dict)
            # Create files for each remaining dimension
            created_files = []
            for dim_combo in dim_combinations:
                selected_data = var_data.sel(dim_combo)
                
                # Create filename with dimension info
                dim_str = "_".join(
                    f"{dim}_{str(val).replace('.', '_')}" 
                    for dim, val in dim_combo.items()
                )
                output_path = output_dir / f"{self.input_path.stem}_{variable}_{dim_str}.{output_format.lower()}"
                
                self._export_data(selected_data, output_path, output_format)
                created_files.append(output_path)
            return created_files
        else:
            # Single file output
            output_path = output_dir / f"{self.input_path.stem}_{variable}.{output_format.lower()}"
            self._export_data(var_data, output_path, output_format)
            return [output_path]
            
    def _export_data(
        self,
        data: xr.DataArray,
        output_path: Path,
        output_format: str,
    ) -> None:
        """Internal method to handle the actual export."""
        if output_format.lower() == "tif":
            self._to_geotiff(data, output_path)
        elif output_format.lower() == "geojson":
            self._to_geojson(data, output_path)
        elif output_format.lower() == "kml":
            self._to_kml(data, output_path)
            
    def _to_geotiff(self, data: xr.DataArray, output_path: Path) -> None:
        """Export data to GeoTIFF format."""
        data.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
        data.rio.write_crs("epsg:4326", inplace=True)
        data.rio.to_raster(output_path)
        
    def _to_geojson(self, data: xr.DataArray, output_path: Path) -> None:
        """Export data to GeoJSON format."""
        raise NotImplementedError("GeoJSON export not yet implemented")
        
    def _to_kml(self, data: xr.DataArray, output_path: Path) -> None:
        """Export data to KML format."""
        raise NotImplementedError("KML export not yet implemented")
        
    def close(self) -> None:
        """Close the NetCDF dataset."""
        self.dataset.close()


def interactive_conversion():
    """Run an interactive conversion session."""
    print("NetCDF to Geospatial Converter - Interactive Mode")
    print("-----------------------------------------------")
    
    # Get input file
    while True:
        input_path = input("Enter path to NetCDF file: ").strip()
        if Path(input_path).exists():
            break
        print(f"File not found: {input_path}")
    
    converter = NetCDFConverter(input_path)
    
    try:
        # Show available variables
        print("\nAvailable variables:")
        for i, var in enumerate(converter.variables, 1):
            print(f"{i}. {var}")
        
        # Select variable
        while True:
            var_choice = input("\nEnter variable number or name: ").strip()
            try:
                if var_choice.isdigit():
                    variable = converter.variables[int(var_choice) - 1]
                else:
                    variable = var_choice
                
                var_info = converter.get_variable_info(variable)
                print(f"\nSelected variable: {variable}")
                print(f"Dimensions: {var_info['dims']}")
                print(f"Shape: {var_info['shape']}")
                break
            except (ValueError, IndexError) as e:
                print(f"Invalid selection: {e}")
        
        # Select dimensions if needed
        dim_selections = {}
        for dim in var_info["dims"]:
            if dim not in {"lat", "latitude", "lon", "longitude", "x", "y"}:
                dim_values = converter.dataset[dim].values
                print(f"\nAvailable values for dimension '{dim}':")
                for i, val in enumerate(dim_values, 1):
                    print(f"{i}. {val}")
                
                while True:
                    dim_choice = input(f"Select {dim} (number or value, leave blank for all): ").strip()
                    if not dim_choice:
                        break
                    
                    try:
                        if dim_choice.isdigit():
                            selected = dim_values[int(dim_choice) - 1]
                        else:
                            # Try to convert to the same type as the dimension values
                            if len(dim_values) > 0:
                                val_type = type(dim_values[0])
                                selected = val_type(dim_choice)
                            else:
                                selected = dim_choice
                        
                        dim_selections[dim] = selected
                        break
                    except (ValueError, IndexError) as e:
                        print(f"Invalid selection: {e}")
        
        # Select output format
        formats = ["tif", "geojson", "kml"]
        print("\nAvailable output formats:")
        for i, fmt in enumerate(formats, 1):
            print(f"{i}. {fmt}")
        
        while True:
            fmt_choice = input("\nEnter output format number or name: ").strip()
            try:
                if fmt_choice.isdigit():
                    output_format = formats[int(fmt_choice) - 1]
                else:
                    output_format = fmt_choice.lower()
                
                if output_format in formats:
                    break
                print("Invalid format selected")
            except (ValueError, IndexError) as e:
                print(f"Invalid selection: {e}")
        
        # Select output directory
        output_dir = input("\nEnter output directory (leave blank for current): ").strip()
        output_dir = output_dir if output_dir else "."
        
        # Perform conversion
        print("\nConverting...")
        created_files = converter.convert(
            variable=variable,
            output_format=output_format,
            output_dir=output_dir,
            dim_selections=dim_selections if dim_selections else None,
        )
        
        print("\nConversion complete. Created files:")
        for file in created_files:
            print(f"- {file}")
            
    finally:
        converter.close()


def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description="Convert NetCDF files to geospatial formats (GeoTIFF, GeoJSON, KML)"
    )
    parser.add_argument("input_file", help="Input NetCDF file path")
    parser.add_argument("-v", "--variable", help="Variable to convert", required=True)
    parser.add_argument(
        "-f", "--format", 
        choices=["tif", "geojson", "kml"],
        help="Output format",
        required=True,
    )
    parser.add_argument(
        "-o", "--output-dir", 
        default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "-d", "--dimension", 
        action="append",
        nargs=2,
        metavar=("DIMENSION", "VALUE"),
        help="Dimension selection (can be used multiple times)",
    )
    
    args = parser.parse_args()
    
    try:
        converter = NetCDFConverter(args.input_file)
        
        # Process dimension selections
        dim_selections = {}
        if args.dimension:
            for dim, value in args.dimension:
                # Try to convert to the appropriate type
                var_dims = converter.dimensions[args.variable]
                if dim not in var_dims:
                    raise ValueError(f"Dimension '{dim}' not found in variable '{args.variable}'")
                
                # Get dimension values to determine type
                dim_values = converter.dataset[dim].values
                if len(dim_values) > 0:
                    val_type = type(dim_values[0])
                    try:
                        dim_selections[dim] = val_type(value)
                    except ValueError:
                        # If type conversion fails, keep as string
                        dim_selections[dim] = value
                else:
                    dim_selections[dim] = value
        
        # Perform conversion
        created_files = converter.convert(
            variable=args.variable,
            output_format=args.format,
            output_dir=args.output_dir,
            dim_selections=dim_selections if dim_selections else None,
        )
        
        print("Conversion complete. Created files:")
        for file in created_files:
            print(f"- {file}")
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        if 'converter' in locals():
            converter.close()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main()
    else:
        interactive_conversion()