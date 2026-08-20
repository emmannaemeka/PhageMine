"""Wait for a local Streamlit health endpoint during packaged smoke tests."""
from __future__ import annotations

import argparse
import time
import urllib.request


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument("--timeout", type=float, default=30)
    args = parser.parse_args()
    deadline = time.monotonic() + args.timeout
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(args.url, timeout=1) as response:
                if response.status == 200:
                    print("READY")
                    return 0
        except OSError:
            time.sleep(0.5)
    print("NOT READY")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
