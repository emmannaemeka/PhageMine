"""Build a native, windowless PhageMine desktop bundle with PyInstaller."""
from __future__ import annotations

import os
import platform
from pathlib import Path

from PIL import Image, ImageDraw
from PyInstaller.__main__ import run as pyinstaller_run

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "build" / "desktop-assets"


def create_icon() -> Path:
    ASSETS.mkdir(parents=True, exist_ok=True)
    image = Image.new("RGBA", (512, 512), "#082f49")
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle((32, 32, 480, 480), radius=96, fill="#0f766e")
    points_a, points_b = [], []
    for x in range(112, 401, 12):
        phase = (x - 112) / 42
        import math
        points_a.append((x, 256 + int(105 * math.sin(phase))))
        points_b.append((x, 256 - int(105 * math.sin(phase))))
    draw.line(points_a, fill="#ecfeff", width=18, joint="curve")
    draw.line(points_b, fill="#a5f3fc", width=18, joint="curve")
    for index in range(0, len(points_a), 3):
        draw.line((points_a[index], points_b[index]), fill="#f0fdfa", width=9)
    png = ASSETS / "phagemine.png"; image.save(png)
    if platform.system() == "Windows":
        ico = ASSETS / "phagemine.ico"
        image.save(ico, sizes=[(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
        return ico
    # PyInstaller converts PNG to the platform-native macOS icon format using
    # Pillow.  This is more portable than invoking iconutil during CI builds.
    return png


def main() -> int:
    icon = create_icon()
    separator = os.pathsep
    args = [str(ROOT / "src" / "phagemine" / "gui" / "desktop.py"),
            "--name=PhageMine", "--onedir", "--windowed", "--noconfirm", "--clean",
            f"--paths={ROOT / 'src'}", f"--icon={icon}",
            f"--add-data={ROOT / 'src' / 'phagemine' / 'gui' / 'app.py'}{separator}phagemine/gui",
            "--collect-data=streamlit", "--collect-submodules=phagemine",
            "--copy-metadata=phagemine", "--copy-metadata=streamlit",
            # These are optional packages that Streamlit can discover in a
            # developer environment; PhageMine does not import or use them.
            "--exclude-module=torch", "--exclude-module=scipy",
            "--exclude-module=matplotlib", "--exclude-module=sphinx",
            "--exclude-module=pytest", "--exclude-module=IPython",
            "--exclude-module=PyQt5", "--exclude-module=tkinter",
            "--exclude-module=plotly", "--exclude-module=black",
            "--exclude-module=zmq", "--exclude-module=numba",
            "--exclude-module=llvmlite", "--exclude-module=botocore",
            "--exclude-module=openpyxl", "--exclude-module=tables",
            "--exclude-module=sqlalchemy", "--exclude-module=qtpy",
            "--exclude-module=distributed", "--exclude-module=dask",
            "--exclude-module=lxml", "--exclude-module=cloudpickle",
            "--exclude-module=bokeh", "--exclude-module=panel",
            "--exclude-module=xarray", "--exclude-module=selenium",
            "--exclude-module=tensorflow", "--exclude-module=pydantic",
            "--exclude-module=statsmodels", "--exclude-module=patsy",
            "--exclude-module=lz4"]
    if platform.system() == "Darwin":
        args.append("--osx-bundle-identifier=org.phagemine.desktop")
    pyinstaller_run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
