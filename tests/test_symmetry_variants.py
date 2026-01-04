#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test for trigonal and tetragonal symmetry variants

This test verifies that the code correctly distinguishes between
high and low symmetry variants of trigonal and tetragonal crystal systems,
as described in the issue about independent elastic tensor components.

References:
- https://github.com/jochym/Elastic/issues/[issue_number]
- Elasticity measurements on minerals: A review, EJM 21(3), 2009
- https://github.com/libAtoms/matscipy/blob/master/matscipy/elasticity.py
"""

from __future__ import print_function, division

import sys
import os
import numpy as np
from ase.spacegroup import crystal

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import elastic module functions
from elastic.elastic import (
    get_lattice_type, get_symmetry_function, get_cij_order,
    tetragonal, tetragonal_high, tetragonal_low,
    trigonal, trigonal_high, trigonal_low
)


def test_tetragonal_high_symmetry():
    """
    Test tetragonal high symmetry crystals (space groups 89-142).
    These should have 6 independent elastic constants.
    Example: TiO2 (rutile), space group 136 (P4_2/mnm)
    """
    print("\n" + "="*60)
    print("Testing Tetragonal High Symmetry")
    print("="*60)
    
    # TiO2 rutile structure - space group 136
    a = 4.60
    c = 2.96
    cryst = crystal(['Ti', 'O'], [(0, 0, 0), (0.302, 0.302, 0)],
                    spacegroup=136, cellpar=[a, a, c, 90, 90, 90])
    
    lattype, bravais, sg_name, sg_nr = get_lattice_type(cryst)
    print(f"Structure: TiO2 (rutile)")
    print(f"Space group: {sg_nr} ({sg_name})")
    print(f"Bravais lattice: {bravais}")
    
    # Get symmetry function
    symm_func, axes, _ = get_symmetry_function(cryst)
    print(f"Symmetry function: {symm_func.__name__}")
    
    # Test matrix dimensions
    u = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
    matrix = symm_func(u)
    print(f"Matrix shape: {matrix.shape}")
    
    # Get elastic constant order
    cij_order = get_cij_order(cryst)
    print(f"Elastic constants: {cij_order}")
    print(f"Number of constants: {len(cij_order)}")
    
    # Assertions
    assert bravais == "Tetragonal", f"Expected Tetragonal, got {bravais}"
    assert symm_func == tetragonal_high, f"Expected tetragonal_high for SG {sg_nr}"
    assert matrix.shape == (6, 6), f"Expected (6,6), got {matrix.shape}"
    assert len(cij_order) == 6, f"Expected 6 constants, got {len(cij_order)}"
    expected_order = ('C_11', 'C_33', 'C_12', 'C_13', 'C_44', 'C_66')
    assert cij_order == expected_order, f"Expected {expected_order}, got {cij_order}"
    
    print("✓ Test PASSED")
    return True


def test_trigonal_high_symmetry():
    """
    Test trigonal high symmetry crystals (space groups 149-167).
    These should have 6 independent elastic constants.
    Example: Sb, space group 166 (R-3m)
    """
    print("\n" + "="*60)
    print("Testing Trigonal High Symmetry")
    print("="*60)
    
    # Sb structure - space group 166
    a = 4.48
    c = 11.04
    cryst = crystal(['Sb'], [(0, 0, 0.24098)],
                    spacegroup=166, cellpar=[a, a, c, 90, 90, 120])
    
    lattype, bravais, sg_name, sg_nr = get_lattice_type(cryst)
    print(f"Structure: Sb")
    print(f"Space group: {sg_nr} ({sg_name})")
    print(f"Bravais lattice: {bravais}")
    
    # Get symmetry function
    symm_func, axes, _ = get_symmetry_function(cryst)
    print(f"Symmetry function: {symm_func.__name__}")
    
    # Test matrix dimensions
    u = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
    matrix = symm_func(u)
    print(f"Matrix shape: {matrix.shape}")
    
    # Get elastic constant order
    cij_order = get_cij_order(cryst)
    print(f"Elastic constants: {cij_order}")
    print(f"Number of constants: {len(cij_order)}")
    
    # Assertions
    assert bravais == "Trigonal", f"Expected Trigonal, got {bravais}"
    assert symm_func == trigonal_high, f"Expected trigonal_high for SG {sg_nr}"
    assert matrix.shape == (6, 6), f"Expected (6,6), got {matrix.shape}"
    assert len(cij_order) == 6, f"Expected 6 constants, got {len(cij_order)}"
    expected_order = ('C_11', 'C_33', 'C_12', 'C_13', 'C_44', 'C_14')
    assert cij_order == expected_order, f"Expected {expected_order}, got {cij_order}"
    
    print("✓ Test PASSED")
    return True


def test_symmetry_function_matrix_shapes():
    """
    Test that all symmetry functions return matrices with correct dimensions.
    """
    print("\n" + "="*60)
    print("Testing Symmetry Function Matrix Dimensions")
    print("="*60)
    
    u = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
    
    tests = [
        (tetragonal_high, (6, 6), "tetragonal_high"),
        (tetragonal_low, (6, 7), "tetragonal_low"),
        (trigonal_high, (6, 6), "trigonal_high"),
        (trigonal_low, (6, 6), "trigonal_low"),
    ]
    
    for func, expected_shape, name in tests:
        matrix = func(u)
        print(f"{name:20s}: {matrix.shape} (expected: {expected_shape})")
        assert matrix.shape == expected_shape, \
            f"{name}: Expected {expected_shape}, got {matrix.shape}"
    
    print("✓ All matrix dimensions correct")
    return True


def test_backward_compatibility():
    """
    Test that the old tetragonal() and trigonal() functions still work
    and are aliases for the high symmetry variants.
    """
    print("\n" + "="*60)
    print("Testing Backward Compatibility")
    print("="*60)
    
    u = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01])
    
    # Test tetragonal() is an alias for tetragonal_high()
    matrix_old = tetragonal(u)
    matrix_new = tetragonal_high(u)
    assert np.array_equal(matrix_old, matrix_new), \
        "tetragonal() should equal tetragonal_high()"
    print("✓ tetragonal() == tetragonal_high()")
    
    # Test trigonal() is an alias for trigonal_high()
    matrix_old = trigonal(u)
    matrix_new = trigonal_high(u)
    assert np.array_equal(matrix_old, matrix_new), \
        "trigonal() should equal trigonal_high()"
    print("✓ trigonal() == trigonal_high()")
    
    print("✓ Backward compatibility maintained")
    return True


def main():
    """Run all tests"""
    print("\n" + "#"*60)
    print("# Test Suite: Trigonal and Tetragonal Symmetry Variants")
    print("#"*60)
    
    try:
        test_symmetry_function_matrix_shapes()
        test_tetragonal_high_symmetry()
        test_trigonal_high_symmetry()
        test_backward_compatibility()
        
        print("\n" + "#"*60)
        print("# All tests PASSED ✓")
        print("#"*60)
        return 0
        
    except AssertionError as e:
        print(f"\n❌ Test FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
