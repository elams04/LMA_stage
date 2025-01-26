#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 11:01:27 2024

@author: lac
"""

def write_material_file(filename, nmat, vp, vs, rho, npow, A, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz):
    with open(filename, 'w+') as fid:
        # Write the number of materials
        fid.write(f"{nmat}\n")
        
        # First 9 'L' lines
        for _ in range(9):
            fid.write(f"L {vp:.8f} {vs:.8f} {rho:.8f} 0.000000 0.000000\n")
        
        # 'F' line
        fid.write(f"F {vp:.8f} {vs:.8f} {rho:.8f} 0.000000 0.000000\n")
        
        # Next 'L' lines from 11 to 27
        for _ in range(11, 28):
            fid.write(f"L {vp:.8f} {vs:.8f} {rho:.8f} 0.000000 0.000000\n")
        
        # Header for PML properties
        fid.write("# PML properties\n")
        fid.write("# npow,Apow,posX,widthX,posY,widthY,posZ,widthZ,mat\n")
        
        # Write PML properties using the given parameters
        pml_data = [
            (xmin, 0, ymin, 0, zmin, -dz),
            (xmin, -dx, ymin, -dy, zmin, -dz),
            (xmin, 0, ymin, -dy, zmin, -dz),
            (xmax, +dx, ymin, -dy, zmin, -dz),
            (xmin, +dx, ymin, 0, zmin, -dz),
            (xmax, +dx, ymax, +dy, zmin, -dz),
            (xmin, 0, ymax, +dy, zmin, -dz),
            (xmin, -dx, ymax, +dy, zmin, -dz),
            (xmin, -dx, ymin, 0, zmin, -dz),
            (xmin, -dx, ymin, -dy, zmin, 0),
            (xmin, 0, ymin, -dy, zmin, 0),
            (xmax, +dx, ymin, -dy, zmin, 0),
            (xmax, +dx, ymin, 0, zmin, 0),
            (xmax, +dx, ymin, +dy, zmin, 0),
            (xmin, 0, ymax, +dy, zmin, 0),
            (xmin, -dx, ymax, +dy, zmin, 0),
            (xmin, -dx, ymin, 0, zmin, 0),
            (xmin, 0, ymin, 0, zmax, +dz),
            (xmin, -dx, ymin, -dy, zmax, +dz),
            (xmin, 0, ymin, -dy, zmax, +dz),
            (xmax, +dx, ymin, -dy, zmax, +dz),
            (xmin, +dx, ymin, 0, zmax, +dz),
            (xmax, +dx, ymax, +dy, zmax, +dz),
            (xmin, 0, ymax, +dy, zmax, +dz),
            (xmin, -dx, ymax, +dy, zmax, +dz),
            (xmin, -dx, ymin, 0, zmax, +dz)
        ]
        
        # Iterate over each PML entry and write to the file
        for posX, widthX, posY, widthY, posZ, widthZ in pml_data:
            fid.write(f"{npow} {A:.6f} {posX:.8f} {widthX:.8f} {posY:.8f} {widthY:.8f} {posZ:.8f} {widthZ:.8f} 9\n")

# Example usage with variables:
write_material_file(
    'material.input',
    nmat=27,
    vp=50.0,
    vs=0.0,
    rho=41.0,
    npow=2,
    A=10.0,
    xmin=-500.0,
    xmax=500.0,
    ymin=-500.0,
    ymax=500.0,
    zmin=-500.0,
    zmax=500.0,
    dx=13.33333333,
    dy=13.33333333,
    dz=13.33333333
)
