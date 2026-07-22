import sys
from pathlib import Path

PIPELINE_DIR = Path(__file__).resolve().parent.parent / "pipeline"
SIGNATURES_DIR = Path(__file__).resolve().parent.parent / "signatures"

if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))
