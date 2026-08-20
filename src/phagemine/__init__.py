"""PhageMine: annotation followed by cautious discovery mining."""

__version__ = "1.1.0rc1"


def release_tag() -> str:
    """Return the canonical Git tag for the authoritative package version."""
    prefix, candidate = __version__.split("rc", 1)
    return f"v{prefix}-rc.{candidate}"
