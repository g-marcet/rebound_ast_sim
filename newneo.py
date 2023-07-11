import numpy as np
import pandas as pd
from io import StringIO
import math
import os

cat = pd.read_csv('../neopop.cat', skiprows=5, sep='\s+')

with open('asteroids.cat', 'w') as fp:
    for i in range(len(cat["!Name"])):
        fp.write(str(cat["a"][i]) + "\n")
        fp.write(str(cat["e"][i]) + "\n")
        fp.write(str(cat["i"][i]) + "\n")
        fp.write(str(cat["node"][i]) + "\n")
        fp.write(str(cat["arg"][i]) + "\n")
        fp.write(str(cat["anomaly"][i]) + "\n")