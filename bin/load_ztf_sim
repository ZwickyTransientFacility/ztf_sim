
import sys
from sqlalchemy import create_engine
import pandas as pd


if len(sys.argv) != 2:
    print('Usage: load_sim.py /path/to/sim.db')

engine = create_engine(f'sqlite:///{sys.argv[1]}')

df = pd.read_sql('Summary', engine)

