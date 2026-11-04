# helpers for the curvilinear NHFLOW tests (G 50, phase P3): runs go to ./runs
# binaries: environment REEF3D_BIN, DIVEMESH_BIN (as ../sflow_curvilinear/p1lib.py)
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sflow_curvilinear"))
from p1lib import *
import p1lib
S3 = os.path.dirname(os.path.abspath(__file__))
