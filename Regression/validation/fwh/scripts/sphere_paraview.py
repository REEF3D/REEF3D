# ParaView (pvbatch) pictures and animations of the flow of validation case 15 (sphere, Re = 300),
# from the VTU output (P 10 1, P 35 ts te dt), in the run directory:
#
#   pvbatch sphere_paraview.py [--t0 T] [--t1 T] [--out DIR] [--q Q] [--size W H] [--views q,slice,mean]
#                              [--angle DEG] [--label TEXT] [--dx DX] [--fps N]
#
#   q      hairpin vortices: iso-surface Q = --q (default 0.05 (U/D)^2) of the Q-criterion, coloured by
#          the streamwise velocity, with the sphere and the FW-H box (U 20 in ctrl.txt) as a wire frame
#   slice  vorticity magnitude in the plane of the shedding (it contains x and the direction of the mean
#          lateral force, from the FW-H box momentum balance; --angle overrides, degrees from y towards z)
#   mean   time-averaged streamwise velocity in the same plane with the line u = 0 (recirculation)
#
# Frames DIR/q_####.png and DIR/slice_####.png for all VTU times in [t0, t1] (default: all), made into
# DIR/q.mp4 and DIR/slice.mp4 with ffmpeg when it is found; DIR/mean.png. Only the region
# x in [-2, 10], |y|, |z| < 3 goes into the filters, resampled to a uniform grid of spacing --dx
# (default 0.05 = D/20; 0.025 for D/40; --dx 0 uses the cells of the VTU instead). --fps: frames per
# second of the videos (default 12; ~4 for VTU every 1 s, so that a shedding period lasts ~2 s).
# pvbatch: /Applications/ParaView-5.12.1.app/Contents/bin/pvbatch (macOS); on a cluster the paraview module.
import argparse, glob, math, os, shutil, subprocess
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()     # the cameras below are set explicitly

ap = argparse.ArgumentParser()
ap.add_argument("--t0", type=float, default=-1e30)
ap.add_argument("--t1", type=float, default=1e30)
ap.add_argument("--out", default="paraview")
ap.add_argument("--q", type=float, default=0.05)
ap.add_argument("--size", type=int, nargs=2, default=[1920, 1080])
ap.add_argument("--views", default="q,slice,mean")
ap.add_argument("--angle", type=float, default=None)
ap.add_argument("--label", default="sphere Re = 300")
ap.add_argument("--dx", type=float, default=0.05)
ap.add_argument("--fps", type=float, default=12.0)
args = ap.parse_args()
views = args.views.split(",")
os.makedirs(args.out, exist_ok=True)

DARK = [14/255, 33/255, 49/255]
LIGHT = [245/255, 243/255, 238/255]

# FW-H box and the plane of the shedding
box = None
for l in open("ctrl.txt"):
    s = l.split()
    if len(s) == 8 and s[0] == "U" and s[1] == "20":
        box = [float(x) for x in s[2:8]]
angle = args.angle
fn = "REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat"
if angle is None and os.path.exists(fn):
    rows = []
    for l in open(fn):
        s = l.split()
        try:
            rows.append([float(x) for x in s[:7]])
        except ValueError:
            pass
    rows = rows[len(rows)//2:]
    # mean lateral force of the fluid on the body: -(Lp + Lm); the mean of dM/dt vanishes
    fy = -sum(r[2] + r[5] for r in rows)/len(rows)
    fz = -sum(r[3] + r[6] for r in rows)/len(rows)
    angle = math.degrees(math.atan2(fz, fy))
if angle is None:
    angle = 0.0
a = math.radians(angle)
normal = [0.0, -math.sin(a), math.cos(a)]     # x cross (0, cos a, sin a)
lateral = [0.0, math.cos(a), math.sin(a)]
print(f"shedding plane at {angle:.1f} deg from y towards z, slice normal {normal}")

files = sorted(glob.glob("REEF3D_CFD_VTU/REEF3D-CFD-*.pvtu"))
if not files:
    raise SystemExit("no REEF3D_CFD_VTU/REEF3D-CFD-*.pvtu in this directory")
reader = XMLPartitionedUnstructuredGridReader(FileName=files)
reader.PointArrayStatus = ["velocity", "pressure"]
reader.UpdatePipelineInformation()
times = [t for t in reader.TimestepValues if args.t0 <= t <= args.t1] if reader.TimestepValues else [0.0]
print(f"{len(files)} VTU files, {len(times)} times in the window {times[0]:.2f}..{times[-1]:.2f} s")

B = [-2.0, 10.0, -3.0, 3.0, -3.0, 3.0]
if args.dx > 0.0:
    region = ResampleToImage(Input=reader)
    region.UseInputBounds = 0
    region.SamplingBounds = B
    region.SamplingDimensions = [int(round((B[2*i+1]-B[2*i])/args.dx))+1 for i in range(3)]
else:
    region = ExtractCellsByRegion(Input=reader)
    region.IntersectWith = "Box"
    region.IntersectWith.UseReferenceBounds = 0
    region.IntersectWith.Position = [B[0], B[2], B[4]]
    region.IntersectWith.Length = [B[1]-B[0], B[3]-B[2], B[5]-B[4]]
grad = Gradient(Input=region)
grad.ScalarArray = ["POINTS", "velocity"]
grad.ComputeGradient = 0
grad.ComputeVorticity = 1
grad.ComputeQCriterion = 1
grad.VorticityArrayName = "vorticity"
grad.QCriterionArrayName = "Q"
# Q away from the body: the direct forcing gives a large Q in a shell around the sphere
qcut = Calculator(Input=grad)
qcut.ResultArrayName = "Qflow"
qcut.Function = "Q*max(0,min(1,(mag(coords)-0.55)*20))"

if os.path.exists("floating.stl"):
    body = STLReader(FileNames=["floating.stl"])
else:
    body = Sphere(Radius=0.5, ThetaResolution=64, PhiResolution=64)
fwhbox = None
if box:
    fwhbox = Box(XLength=box[1]-box[0], YLength=box[3]-box[2], ZLength=box[5]-box[4],
                 Center=[0.5*(box[0]+box[1]), 0.5*(box[2]+box[3]), 0.5*(box[4]+box[5])])


def new_view():
    v = CreateView("RenderView")
    v.ViewSize = args.size
    v.OrientationAxesVisibility = 0
    v.UseColorPaletteForBackground = 0
    v.Background = DARK
    return v


def add_body(v, wire_color, cut=False):
    d = Show(body, v)
    d.ColorArrayName = ["POINTS", ""]
    d.DiffuseColor = [0.80, 0.80, 0.80]
    d.AmbientColor = [0.80, 0.80, 0.80]
    if fwhbox:
        # 3D: the box as a wire frame; in a slice view: its section with the plane
        src = fwhbox
        if cut:
            src = Slice(Input=fwhbox)
            src.SliceType = "Plane"
            src.SliceType.Origin = [0.0, 0.0, 0.0]
            src.SliceType.Normal = normal
        b = Show(src, v)
        b.ColorArrayName = ["POINTS", ""]
        b.Representation = "Wireframe"
        b.AmbientColor = wire_color
        b.DiffuseColor = wire_color
        b.LineWidth = 2.5


def add_time(v, color=LIGHT):
    t = AnnotateTime()
    t.Format = args.label + "   t = {time:.1f} s"
    d = Show(t, v)
    d.FontSize = 26
    d.Color = color
    d.WindowLocation = "Upper Left Corner"


def bar(v, lut, title, color=LIGHT):
    sb = GetScalarBar(lut, v)
    sb.Title = title
    sb.ComponentTitle = ""
    sb.TitleColor = color
    sb.LabelColor = color
    sb.TitleFontSize = 22
    sb.LabelFontSize = 20
    sb.WindowLocation = "Any Location"
    sb.Orientation = "Horizontal"
    sb.Position = [0.33, 0.05]
    sb.ScalarBarLength = 0.34
    sb.RangeLabelFormat = "%-#.2g"
    sb.Visibility = 1


def frames(v, name):
    for i, t in enumerate(times):
        v.ViewTime = t
        Render(v)
        SaveScreenshot(os.path.join(args.out, f"{name}_{i:04d}.png"), v, ImageResolution=args.size)
        print(f"{name} {i+1}/{len(times)}  t = {t:.2f}", flush=True)
    ff = shutil.which("ffmpeg") or ("/opt/homebrew/bin/ffmpeg" if os.path.exists("/opt/homebrew/bin/ffmpeg") else None)
    if ff and len(times) > 1:
        subprocess.run([ff, "-loglevel", "error", "-y", "-framerate", f"{args.fps:g}", "-i",
                        os.path.join(args.out, f"{name}_%04d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p",
                        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", os.path.join(args.out, f"{name}.mp4")], check=False)
        print("written", os.path.join(args.out, f"{name}.mp4"))


# 1 Q-criterion
if "q" in views:
    v = new_view()
    con = Contour(Input=qcut)
    con.ContourBy = ["POINTS", "Qflow"]
    con.Isosurfaces = [args.q]
    con.ComputeScalars = 0
    d = Show(con, v)
    ColorBy(d, ("POINTS", "velocity", "X"))
    lut = GetColorTransferFunction("velocity")
    lut.ApplyPreset("Cool to Warm", True)
    lut.VectorMode = "Component"
    lut.VectorComponent = 0
    lut.RescaleTransferFunction(-0.2, 1.2)
    lut.AutomaticRescaleRangeMode = "Never"
    d.Specular = 0.3
    bar(v, lut, "u [m/s]")
    add_body(v, [0.91, 0.54, 0.31])
    add_time(v)
    # view from the side of the shedding plane, slightly from above and from upstream
    c = [3.5, 0.0, 0.0]
    v.CameraFocalPoint = c
    v.CameraPosition = [c[0] - 5.0 + 14.0*normal[0], 14.0*normal[1] + 4.0*lateral[1], 14.0*normal[2] + 4.0*lateral[2]]
    v.CameraViewUp = lateral
    v.CameraViewAngle = 30
    frames(v, "q")
    Delete(v)

# 2 vorticity in the shedding plane
if "slice" in views:
    v = new_view()
    sl = Slice(Input=grad)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.0, 0.0, 0.0]
    sl.SliceType.Normal = normal
    sl.Triangulatetheslice = 0
    calc = Calculator(Input=sl)
    calc.ResultArrayName = "vort_n"
    calc.Function = f"vorticity_X*{normal[0]}+vorticity_Y*{normal[1]}+vorticity_Z*{normal[2]}"
    d = Show(calc, v)
    ColorBy(d, ("POINTS", "vort_n"))
    lut = GetColorTransferFunction("vort_n")
    lut.ApplyPreset("Blue Orange (divergent)", True)
    lut.RescaleTransferFunction(-3.0, 3.0)
    lut.AutomaticRescaleRangeMode = "Never"
    bar(v, lut, "vorticity normal to the plane [1/s]", DARK)
    add_body(v, DARK, cut=True)
    add_time(v, DARK)
    c = [4.0, 0.0, 0.0]
    v.CameraFocalPoint = c
    v.CameraPosition = [c[i] + 14.0*normal[i] for i in range(3)]
    v.CameraViewUp = lateral
    v.CameraParallelProjection = 1
    v.CameraParallelScale = 3.2
    frames(v, "slice")
    Delete(v)

# 3 mean streamwise velocity and the recirculation (u = 0)
if "mean" in views and len(times) > 1:
    v = new_view()
    ext = ExtractTimeSteps(Input=region)
    ext.SelectionMode = "Select Time Range"
    all_t = list(reader.TimestepValues)
    ext.TimeStepRange = [all_t.index(times[0]), all_t.index(times[-1])]
    stat = TemporalStatistics(Input=ext)
    stat.ComputeMinimum = 0
    stat.ComputeMaximum = 0
    stat.ComputeStandardDeviation = 0
    sl = Slice(Input=stat)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.0, 0.0, 0.0]
    sl.SliceType.Normal = normal
    sl.Triangulatetheslice = 0
    calc = Calculator(Input=sl)
    calc.ResultArrayName = "u_mean"
    calc.Function = "velocity_average_X"
    d = Show(calc, v)
    ColorBy(d, ("POINTS", "u_mean"))
    lut = GetColorTransferFunction("u_mean")
    lut.ApplyPreset("Cool to Warm", True)
    lut.RescaleTransferFunction(-0.4, 1.2)
    lut.AutomaticRescaleRangeMode = "Never"
    bar(v, lut, f"mean u [m/s], {times[0]:.0f}-{times[-1]:.0f} s")
    zero = Contour(Input=calc)
    zero.ContourBy = ["POINTS", "u_mean"]
    zero.Isosurfaces = [0.0]
    z = Show(zero, v)
    z.ColorArrayName = ["POINTS", ""]
    z.AmbientColor = [0.0, 0.0, 0.0]
    z.DiffuseColor = [0.0, 0.0, 0.0]
    z.LineWidth = 3.0
    add_body(v, LIGHT, cut=True)
    c = [4.0, 0.0, 0.0]
    v.CameraFocalPoint = c
    v.CameraPosition = [c[i] + 14.0*normal[i] for i in range(3)]
    v.CameraViewUp = lateral
    v.CameraParallelProjection = 1
    v.CameraParallelScale = 3.2
    Render(v)
    SaveScreenshot(os.path.join(args.out, "mean.png"), v, ImageResolution=args.size)
    print("written", os.path.join(args.out, "mean.png"))
    # length of the recirculation behind the sphere: the largest x of u_mean = 0 near the axis
    zero.UpdatePipeline()
    from paraview import servermanager as sm
    data = sm.Fetch(zero)
    xs = [data.GetPoint(i)[0] for i in range(data.GetNumberOfPoints())
          if abs(sum(data.GetPoint(i)[k]*lateral[k] for k in range(3))) < 0.6]
    if xs:
        print(f"recirculation ends at x = {max(xs):.3f} D (length behind the sphere {max(xs)-0.5:.3f} D)")
