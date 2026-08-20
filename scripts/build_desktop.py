"""Build a native, windowless PhageMine desktop bundle with PyInstaller."""
from __future__ import annotations

import os
import platform
import plistlib
import hashlib
import shutil
import urllib.request
from importlib.metadata import distribution
from pathlib import Path

from PIL import Image, ImageDraw
from PyInstaller.__main__ import run as pyinstaller_run
from phagemine import __version__

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "build" / "desktop-assets"
BACKEND_NAME = "PhageMine-PHANOTATE"
UPSTREAM_SOURCES = {
    "phanotate-1.6.7.tar.gz": (
        "https://files.pythonhosted.org/packages/f2/85/d30456168ba839c9cd97ea241e245fef7d6cbdeaa931504c51447f0dd9a1/phanotate-1.6.7.tar.gz",
        "81889ccaa04bd530a7f2f0b93064c748940b573d736c824f2e97450df5b2033e",
    ),
    "fastpath-1.9.tar.gz": (
        "https://files.pythonhosted.org/packages/75/14/0cfe89b016c8f82ff9db25c84a64c146b0e5ec988a1ba2493fd29d109198/fastpath-1.9.tar.gz",
        "3372d306a3c4e4e764b3995946132333726a229e9002879b9112779dd442b31a",
    ),
    "genbank-0.121.tar.gz": (
        "https://files.pythonhosted.org/packages/f7/b8/51598260fadedfe04ad26195597148bd6c08f8e2f7178e73fc3de16683b5/genbank-0.121.tar.gz",
        "959df344018f8e104ef6d3c9d78f576228cb4997f4fa186050cfa69b231cb5af",
    ),
}


def prepare_corresponding_source() -> Path:
    """Download and verify the exact GPL source shipped with the backend."""
    destination = ASSETS / "third-party-source"
    destination.mkdir(parents=True, exist_ok=True)
    for filename, (url, expected) in UPSTREAM_SOURCES.items():
        path = destination / filename
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            with urllib.request.urlopen(url, timeout=60) as response:
                payload = response.read()
            if hashlib.sha256(payload).hexdigest() != expected:
                raise RuntimeError(f"Checksum mismatch for upstream source: {filename}")
            path.write_bytes(payload)
    return destination


def distribution_file(package: str, predicate) -> Path:
    installed = distribution(package)
    for item in installed.files or []:
        path = Path(installed.locate_file(item))
        if predicate(item, path):
            return path
    raise RuntimeError(f"Required {package} distribution file was not found")


def phanotate_program() -> Path:
    return distribution_file(
        "phanotate", lambda item, path: path.name == "phanotate.py"
    )


def license_file(package: str) -> Path:
    return distribution_file(package, lambda item, path: "license" in path.name.lower())


def windows_version_file() -> Path:
    path = ASSETS / "windows-version.txt"
    path.write_text(
        "VSVersionInfo(ffi=FixedFileInfo(filevers=(1,1,0,1), prodvers=(1,1,0,1), "
        "mask=0x3f, flags=0x0, OS=0x40004, fileType=0x1, subtype=0x0, date=(0,0)), "
        "kids=[StringFileInfo([StringTable('040904B0', ["
        f"StringStruct('FileDescription', 'PhageMine Desktop'), StringStruct('FileVersion', '{__version__}'), "
        f"StringStruct('ProductName', 'PhageMine'), StringStruct('ProductVersion', '{__version__}')])]), "
        "VarFileInfo([VarStruct('Translation', [1033, 1200])])])\n"
    )
    return path


def build_phanotate_backend() -> Path | None:
    """Freeze upstream PHANOTATE as a separate GPL-licensed process."""
    if platform.system() == "Windows":
        return None
    pyinstaller_run([
        str(ROOT / "packaging" / "phanotate_backend.py"),
        f"--name={BACKEND_NAME}", "--onedir", "--console", "--noconfirm", "--clean",
        f"--add-data={phanotate_program()}{os.pathsep}.",
        "--collect-submodules=phanotate_modules", "--hidden-import=fastpathz",
        "--hidden-import=genbank", "--copy-metadata=phanotate",
        "--copy-metadata=fastpath", "--copy-metadata=genbank",
    ])
    suffix = ".exe" if platform.system() == "Windows" else ""
    bundle = ROOT / "dist" / BACKEND_NAME
    backend = bundle / f"{BACKEND_NAME}{suffix}"
    if not backend.is_file():
        raise RuntimeError(f"PHANOTATE backend build missing: {backend}")
    stable_backend = ASSETS / BACKEND_NAME
    if stable_backend.exists():
        shutil.rmtree(stable_backend)
    shutil.copytree(bundle, stable_backend)
    return stable_backend


def apply_macos_version() -> None:
    plist = ROOT / "dist" / "PhageMine.app" / "Contents" / "Info.plist"
    with plist.open("rb") as handle:
        payload = plistlib.load(handle)
    payload["CFBundleShortVersionString"] = __version__
    payload["CFBundleVersion"] = "1.1.0.1"
    with plist.open("wb") as handle:
        plistlib.dump(payload, handle)


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
    backend = build_phanotate_backend()
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
    if backend is not None:
        sources = prepare_corresponding_source()
        args.extend([
            f"--add-data={backend}{separator}phagemine_backend",
            f"--add-data={ROOT / 'docs' / 'third-party' / 'PHANOTATE.md'}{separator}licenses",
            f"--add-data={sources}{separator}third-party-source",
        ])
        for package in ("phanotate", "fastpath", "genbank"):
            args.append(f"--add-data={license_file(package)}{separator}licenses/{package}")
    if platform.system() == "Darwin":
        args.append("--osx-bundle-identifier=org.phagemine.desktop")
    if platform.system() == "Windows":
        args.append(f"--version-file={windows_version_file()}")
    pyinstaller_run(args)
    if platform.system() == "Darwin":
        apply_macos_version()
    if backend is not None:
        shutil.rmtree(backend, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
