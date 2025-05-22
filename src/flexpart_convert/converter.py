import itertools
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import geojson
import numpy as np
import xarray as xr
import rioxarray as rio

class NetCDFConversionError(Exception):
    """Base exception for conversion errors."""
    pass

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
        output_format: str = "geotiff",
        output_dir: Union[str, Path] = ".",
        dim_selections: Optional[Dict[str, Union[int, float, str]]] = None,
    ) -> List[Path]:
        """Convert the NetCDF variable to the specified format.
        
        Args:
            variable: Variable name to convert
            output_format: One of 'geotiff', 'geojson', or 'kml'
            output_dir: Directory to save output files
            dim_selections: Dictionary of dimension selections (e.g., {'time': 0})
            
        Returns:
            List of paths to created files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if output_format.lower() not in {"geotiff", "geojson", "kml"}:
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
        if output_format.lower() == "geotiff":
            self._to_geotiff(data, output_path)
        elif output_format.lower() == "geojson":
            self._to_geojson(data, output_path)
        elif output_format.lower() == "kml":
            self._to_kml(data, output_path)
            
    def _to_geotiff(self, data: xr.DataArray, output_path: Path) -> None:
        """Export data to GeoTIFF format."""
        data.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
        data.rio.write_crs("epsg:4326", inplace=True)
        data.rio.to_raster(output_path, driver="GTiff")
        
    def _to_geojson(self, data: xr.DataArray, output_path: Path) -> None:
        """Export data to GeoJSON format."""
        
        # To improve performance, we trim the raster, getting rid of the bounding zero values
        data = data.where(data > 0, drop=True)
        lons = data["longitude"].values
        lats = data["latitude"].values
        dx = (lons[1] - lons[0]); dy = (lats[1] - lats[0])
        geojson_features = []
        
        for x in lons:
            for y in lats:
                value = data.sel(longitude=x, latitude=y)
                if value == 0.0: # we don't write the zero values
                    continue
                dx2 = dx/2
                dy2 = dy/2
                left = x - dx2; right = x + dx2; lower = y - dy2; upper = y + dy2
                cell_coords = (
                    (float(left), float(lower)),
                    (float(left), float(upper)),
                    (float(right), float(upper)),
                    (float(right), float(lower)),
                )
                polygon = geojson.Polygon(cell_coords)
                properties = {"value": float(value)}
                feature = geojson.Feature(geometry=polygon, properties=properties)
                geojson_features.append(feature)
                
        with open(output_path, "w") as f:
            geojson.dump(geojson.FeatureCollection(geojson_features), f)
        
    def _to_kml(self, data: xr.DataArray, output_path: Path) -> None:
        """Export data to KML format."""
        raise NotImplementedError("KML export not yet implemented")
        
    def close(self) -> None:
        """Close the NetCDF dataset."""
        self.dataset.close()

