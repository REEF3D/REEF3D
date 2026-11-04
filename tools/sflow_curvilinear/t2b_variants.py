# T2 variants: fully wet on the wavy grid, island and fully wet on the logical (Cartesian) nodes
import sys, numpy as np
from p1lib import *
exec(open(os.path.join(S,"t2_lake_at_rest.py")).read().split("for name, grid in")[0].replace("steps = int(sys.argv[1]) if len(sys.argv) > 1 else 200","steps=300"))
Zs = 0.6*np.exp(-((X-17.0)**2+(Y-17.5)**2)/12.0) + 0.05*rng.standard_normal(X.shape)
Zi = 1.3*np.exp(-((x-17.0)**2+(y-17.5)**2)/12.0) + 0.05*rng.standard_normal(X.shape)
Zis = 0.6*np.exp(-((x-17.0)**2+(y-17.5)**2)/12.0) + 0.05*rng.standard_normal(X.shape)
for name, grid in [("wavy_sub",(X,Y,Zs)),("id_island",(x,y,Zi)),("id_sub",(x,y,Zis))]:
    d = os.path.join(S, "runs/t2_"+name)
    rc = run(d, control, ctrl, grid=grid, steps=steps)
    st = final_state(d)
    wet = st["WL"] > 0.002
    print("%-10s rc=%d wet %d  max |P| %.3e  max |Q| %.3e  max |eta| %.3e" % (name, rc, wet.sum(), np.abs(st["P"]).max(), np.abs(st["Q"]).max(), np.abs(st["eta"][wet]).max()))
