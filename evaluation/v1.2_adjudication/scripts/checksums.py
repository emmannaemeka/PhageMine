"""Write checksums for adjudication infrastructure, excluding the manifest itself."""
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
files=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.name not in {'SHA256SUMS'} and '__pycache__' not in p.parts and not str(p).endswith('.pyc'))
(ROOT/'SHA256SUMS').write_text(''.join(f'{sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT.parent.parent)}\n' for p in files))
