"""
Библиотека для оценки видимости космических аппаратов
"""

from .core import (
    j2000_to_cirs_precession_only,
    geodetic_to_itrs,
    itrs_to_cirs,
    calculate_elevation,
    check_visibility
)

from .visualization import (
    visualize_satellite_visibility,
    check_visibility_with_plot
)

__version__ = "1.0.0"
__author__ = "Ryazanov Daniil"

__all__ = [
    'j2000_to_cirs_precession_only',
    'geodetic_to_itrs', 
    'itrs_to_cirs',
    'calculate_elevation',
    'check_visibility',
    'visualize_satellite_visibility',
    'check_visibility_with_plot'
]