"""Refresh manifests after regeneration; paths are relative to repository root."""
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[3]
for directory in ['docs/benchmark_v1.2','evaluation/v1.2_functional']:
    folder=ROOT/directory
    files=sorted(p for p in folder.rglob('*') if p.is_file() and p.name!='SHA256SUMS' and '__pycache__' not in p.parts and p.suffix!='.pyc')
    (folder/'SHA256SUMS').write_text(''.join(f'{sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}\n' for p in files))
