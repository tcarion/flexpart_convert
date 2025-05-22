import argparse
import sys
from pathlib import Path
from typing import Dict, Optional, Union

from .converter import NetCDFConverter, NetCDFConversionError

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
    parser.add_argument(
        "-v", "--variable",
        default="spec001_mr",
        help="Variable to convert",
        required=True
    )
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
        # Process dimension selections
        dim_selections = _process_dim_selections(args.dimension, args.variable, converter)
        
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



def _process_dim_selections(
    dim_args: Optional[list],
    variable: str,
    converter: NetCDFConverter,
) -> Dict[str, Union[int, float, str]]:
    """Process dimension selections from command line arguments."""
    if not dim_args:
        return {}
    
    dim_selections = {}
    var_dims = converter.dimensions[variable]
    
    for dim, value in dim_args:
        if dim not in var_dims:
            raise NetCDFConversionError(f"Dimension '{dim}' not found in variable '{variable}'")
        
        # Try to convert to the appropriate type
        dim_values = converter.dataset[dim].values
        if len(dim_values) > 0:
            try:
                val_type = type(dim_values[0])
                dim_selections[dim] = val_type(value)
            except ValueError:
                dim_selections[dim] = value
        else:
            dim_selections[dim] = value
    
    return dim_selections

if __name__ == "__main__":
    if len(sys.argv) > 1:
        main()
    else:
        interactive_conversion()