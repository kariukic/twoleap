#!/usr/bin/env python3
"""
get_subtraction_directions.py
Read direction names from sky model for subtraction.
"""

import sys
import os

def read_directions_from_skymodel(sky_model_file, exclude_directions=None):
    """Read direction names from a sky model file."""
    directions = []
    
    try:
        # Try to use lsmtool if available
        import lsmtool
        sky = lsmtool.load(sky_model_file)
        directions = sky.getPatchNames().tolist()
        print(f"Using lsmtool: Found {len(directions)} directions", file=sys.stderr)
    except ImportError:
        # Fallback: parse file manually
        print("lsmtool not available, parsing sky model manually", file=sys.stderr)
        with open(sky_model_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(','):
                    parts = line.split(',')
                    if len(parts) > 1:
                        dir_name = parts[1].strip()
                        if dir_name:
                            directions.append(dir_name)
    
    # Filter out excluded directions
    if exclude_directions:
        exclude_list = [d.strip() for d in exclude_directions.split(',') if d.strip()]
        directions = [d for d in directions if d not in exclude_list]
    
    return directions

def main():
    if len(sys.argv) < 2:
        print("Usage: get_subtraction_directions.py <sky_model_file> [exclude_directions]")
        sys.exit(1)
    
    sky_model_file = sys.argv[1]
    exclude_directions = sys.argv[2] if len(sys.argv) > 2 else ""
    
    try:
        directions = read_directions_from_skymodel(sky_model_file, exclude_directions)

        directions = [[d] for d in directions]

        print(str(directions).replace(" ", ""))

        
        # # Format for DP3: [dir1,dir2,dir3]
        # if directions:
        #     dp3_format = '[' + ','.join(directions) + ']'
        # else:
        #     dp3_format = '[]'
        
        # print(dp3_format)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        print("[]")  # Return empty list on error
        sys.exit(1)

if __name__ == "__main__":
    main()