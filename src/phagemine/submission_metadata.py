"""User-supplied submission metadata; no biological or contact fields inferred."""
import json
from pathlib import Path
from .models import SubmissionMetadata


def load_metadata(path=None, **overrides):
    data={}
    if path:
        p=Path(path)
        if p.suffix.lower() in {'.yaml','.yml'}:
            import yaml
            data=yaml.safe_load(p.read_text()) or {}
        else:data=json.loads(p.read_text())
        if not isinstance(data,dict):raise ValueError('Submission metadata must be a mapping')
    if 'country' in data:data.setdefault('geographic_location',data.pop('country'))
    for key,value in overrides.items():
        if value is not None:data[key]=value
    if data.get('topology') is None: data['topology']='linear'
    if data['topology'] not in {'linear','circular'}:raise ValueError('Topology must be linear or circular')
    template=data.get('submitter_template')
    if template:
        p=Path(template).expanduser()
        if not p.is_absolute() and path:p=Path(path).resolve().parent/p
        if not p.is_file():raise ValueError(f'Submitter template not found: {p}')
        data['submitter_template']=str(p.resolve())
    return SubmissionMetadata.from_dict(data)
