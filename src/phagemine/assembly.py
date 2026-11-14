"""Optional paired-end assembly; never overwrites an existing assembly directory."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess


def assemble(reads1, reads2, output, assembler='shovill', executable=None, threads=1):
    if assembler not in {'shovill','spades'} or threads < 1:
        raise ValueError('Choose shovill or spades and positive threads')
    inputs=[Path(reads1).resolve(),Path(reads2).resolve()]
    if not all(p.is_file() for p in inputs): raise ValueError('Both read files must exist')
    root=Path(output).resolve()
    if root.exists(): raise ValueError('Assembly output must not already exist')
    exe=shutil.which(executable or ('shovill' if assembler=='shovill' else 'spades.py'))
    if not exe: raise ValueError(f'{assembler} executable is not available')
    command=([exe,'--R1',str(inputs[0]),'--R2',str(inputs[1]),'--outdir',str(root),'--depth','100','--cpus',str(threads)]
             if assembler=='shovill' else [exe,'-1',str(inputs[0]),'-2',str(inputs[1]),'-o',str(root),'-t',str(threads)])
    def checksum(path):
        h=hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
        return h.hexdigest()
    hashes={str(p):checksum(p) for p in inputs}
    result=subprocess.run(command,capture_output=True,text=True,check=False)
    root.mkdir(parents=True,exist_ok=True)
    fasta=root/'contigs.fasta'
    report={'command':command,'returncode':result.returncode,'input_sha256':hashes,'stdout':result.stdout,'stderr':result.stderr,'coverage_verified':False,'depth_cap':100 if assembler=='shovill' else None}
    (root/'assembly_provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    if result.returncode or not fasta.is_file(): raise RuntimeError(f'Assembly failed; inspect {root / "assembly_provenance.json"}')
    return fasta
