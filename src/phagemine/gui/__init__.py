"""Optional local graphical interface for PhageMine.

This package deliberately does not import Streamlit so core installations can
import :mod:`phagemine.gui` without installing GUI dependencies.
"""

GUI_VERSION = "0.1"

__all__ = ["GUI_VERSION"]
