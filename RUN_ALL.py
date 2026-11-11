"""End-to-end controller. GEE extraction/download and analysis are repeatable."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from ENGINE_RUN_ALL import main

if __name__=='__main__': main()
