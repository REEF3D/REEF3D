# pvbatch kp505_q_render.py DIR [--q 2000] [--size 1600 900] [--fps 12]
# Renders the frames of kp505_q_export.py: iso-surface Q = --q coloured by the axial velocity, the
# propeller surface; DIR/q3d_####.png, DIR/kp505_q3d.mp4 (ffmpeg) and a GIF of one revolution.
import argparse, glob, os, shutil, subprocess
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

ap = argparse.ArgumentParser()
ap.add_argument("dir")
ap.add_argument("--q", type=float, default=400.0)
ap.add_argument("--size", type=int, nargs=2, default=[1600, 900])
ap.add_argument("--fps", type=float, default=12.0)
ap.add_argument("--umax", type=float, default=3.0)
args = ap.parse_args()

grids = sorted(glob.glob(os.path.join(args.dir, "q_*.vtk")))
bodies = sorted(glob.glob(os.path.join(args.dir, "body_*.vtk")))
times = {}
if os.path.exists(os.path.join(args.dir, "times.txt")):
    for l in open(os.path.join(args.dir, "times.txt")):
        j, t = l.split()
        times[int(j)] = float(t)

DARK = [14/255, 33/255, 49/255]
LIGHT = [245/255, 243/255, 238/255]
v = CreateView("RenderView")
v.ViewSize = args.size
v.OrientationAxesVisibility = 0
v.UseColorPaletteForBackground = 0
v.Background = DARK

grid = LegacyVTKReader(FileNames=[grids[0]])
con = Contour(Input=grid)
con.ContourBy = ["POINTS", "Q"]
con.Isosurfaces = [args.q]
con.ComputeScalars = 0
d = Show(con, v)
ColorBy(d, ("POINTS", "u"))
lut = GetColorTransferFunction("u")
lut.ApplyPreset("Cool to Warm", True)
lut.RescaleTransferFunction(0.0, args.umax)
lut.AutomaticRescaleRangeMode = "Never"
d.Specular = 0.3
sb = GetScalarBar(lut, v)
sb.Title = "axial velocity [m/s]"
sb.ComponentTitle = ""
sb.TitleColor = LIGHT
sb.LabelColor = LIGHT
sb.TitleFontSize = 20
sb.LabelFontSize = 18
sb.WindowLocation = "Any Location"
sb.Orientation = "Horizontal"
sb.Position = [0.33, 0.04]
sb.ScalarBarLength = 0.34
sb.Visibility = 1

body = LegacyVTKReader(FileNames=[bodies[0]])
b = Show(body, v)
b.ColorArrayName = ["POINTS", ""]
b.DiffuseColor = [0.82, 0.78, 0.62]
b.AmbientColor = [0.82, 0.78, 0.62]
b.Specular = 0.5

txt = Text()
td = Show(txt, v)
td.FontSize = 22
td.Color = LIGHT
td.WindowLocation = "Upper Left Corner"

# oblique view from the side and slightly upstream, the wake to the right
v.CameraFocalPoint = [0.07, 0.0, 0.0]
v.CameraPosition = [-0.13, -0.66, 0.27]
v.CameraViewUp = [0.0, 0.0, 1.0]
v.CameraViewAngle = 30

for j, (g, bb) in enumerate(zip(grids, bodies)):
    grid.FileNames = [g]
    body.FileNames = [bb]
    txt.Text = f"KP505, J = 0.6, iso-surface Q = {args.q:g} 1/s²     t = {times.get(j, 0.0):.3f} s"
    Render(v)
    SaveScreenshot(os.path.join(args.dir, f"q3d_{j:04d}.png"), v, ImageResolution=args.size)
    print(f"render {j+1}/{len(grids)}", flush=True)

ff = shutil.which("ffmpeg") or "/opt/homebrew/bin/ffmpeg"
subprocess.run([ff, "-loglevel", "error", "-y", "-framerate", f"{args.fps:g}", "-i", os.path.join(args.dir, "q3d_%04d.png"),
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
                os.path.join(args.dir, "kp505_q3d.mp4")], check=True)
nrev = 20 if len(grids) >= 20 else len(grids)
subprocess.run([ff, "-loglevel", "error", "-y", "-framerate", f"{args.fps:g}", "-i", os.path.join(args.dir, "q3d_%04d.png"),
                "-frames:v", str(nrev), "-vf",
                "scale=960:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];[b][p]paletteuse=dither=none",
                "-loop", "0", os.path.join(args.dir, "kp505_q3d.gif")], check=True)
print("written", os.path.join(args.dir, "kp505_q3d.mp4"))
