from phagemine.pipeline import _consensus_provider_ids


def test_shared_consensus_provider_set_for_all_profiles():
    expected = ("phanotate", "pyrodigal", "prodigal_gv")
    assert _consensus_provider_ids("standard") == expected
    assert _consensus_provider_ids("extended") == expected


def test_command_line_prodigal_is_not_default_consensus_provider():
    assert "prodigal" not in _consensus_provider_ids()
