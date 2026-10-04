import hashlib
from unittest.mock import patch

from phagemine.swissprot import SwissProtEvidenceAdapter


def metadata(path, product):
    path.write_text('ID   ENTRY Reviewed; 100 AA.\nAC   P00001;\nDE   RecName: Full='
                    + product + ';\nOS   Bacteriophage example.\n//\n')


def test_repeated_hits_hash_metadata_once(tmp_path):
    source = tmp_path / 'reviewed.dat'
    metadata(source, 'capsid protein')
    adapter = SwissProtEvidenceAdapter(metadata_path=source)
    with patch('phagemine.swissprot.hashlib.sha256', wraps=hashlib.sha256) as digest:
        for _ in range(20):
            assert adapter._metadata_for('P00001')['protein_name'] == 'capsid protein'
        assert digest.call_count == 1


def test_same_size_metadata_change_invalidates_verified_index(tmp_path):
    source = tmp_path / 'reviewed.dat'
    metadata(source, 'capsid protein')
    adapter = SwissProtEvidenceAdapter(metadata_path=source)
    assert adapter._metadata_for('P00001')['protein_name'] == 'capsid protein'
    old_size = source.stat().st_size
    metadata(source, 'portal protein')
    assert source.stat().st_size == old_size
    assert adapter._metadata_for('P00001')['protein_name'] == 'portal protein'


def test_deleted_index_is_rebuilt(tmp_path):
    source = tmp_path / 'reviewed.dat'
    metadata(source, 'capsid protein')
    adapter = SwissProtEvidenceAdapter(metadata_path=source)
    assert adapter._metadata_for('P00001')['protein_name'] == 'capsid protein'
    adapter._metadata_index().unlink()
    assert adapter._metadata_for('P00001')['protein_name'] == 'capsid protein'
