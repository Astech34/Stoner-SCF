import os
import numpy as np
import plotly.graph_objects as go

# ---------------------------------------------------------------------------
# t2g orbital wavefunctions (3d hydrogen-like radial part: r^2 exp(-r/3))
# The r^2 from R_{3,2} cancels the r^2 in the denominator of the angular part,
# leaving just the product of two Cartesian coordinates times exp(-r/3).
# ---------------------------------------------------------------------------
def orbital(name, X, Y, Z, center=(0, 0, 0)):
    x = X - center[0]
    y = Y - center[1]
    z = Z - center[2]
    R = np.sqrt(x**2 + y**2 + z**2)
    radial = np.exp(-R / 3)
    if name == 'xy':
        return x * y * radial
    elif name == 'xz':
        return x * z * radial
    elif name == 'yz':
        return y * z * radial
    else:
        raise ValueError(f"Unknown orbital '{name}'. Choose from: xy, xz, yz")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
orbital_1 = 'xy'
orbital_2 = 'xy'

lattice_const = 20.0                    # square lattice constant (in same units as grid)
shift = (0.0, lattice_const, 0.0)      # nearest-neighbour shift along x

site_1 = (0.0, 0.0, 0.0)
site_2 = (site_1[0] + shift[0], site_1[1] + shift[1], site_1[2] + shift[2])

threshold = 1.0

# ---------------------------------------------------------------------------
# Grid and orbital evaluation
# ---------------------------------------------------------------------------
coords = np.linspace(-20, 50, 70)
X, Y, Z = np.meshgrid(coords, coords, coords)

V1 = orbital(orbital_1, X, Y, Z, center=site_1)
V2 = orbital(orbital_2, X, Y, Z, center=site_2)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
colorscale = [
    [0.0, '#4B0082'],   # deep indigo  (negative)
    [0.5, '#F5F5F5'],   # near-white   (zero)
    [1.0, '#FF6F00'],   # deep amber   (positive)
]

def make_isosurface(V, name, show_scale=False):
    return go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=V.flatten(),
        isomin=-threshold,
        isomax=threshold,
        surface_count=2,
        colorscale=colorscale,
        opacity=0.4,
        caps=dict(x_show=False, y_show=False, z_show=False),
        name=name,
        showscale=show_scale,
        colorbar=dict(title='+/−', thickness=15, len=0.5) if show_scale else None
    )

def axis_arrow(end, label):
    """Single arrow from origin to `end` with a text label."""
    x, y, z = end
    return [
        go.Scatter3d(
            x=[0, x], y=[0, y], z=[0, z],
            mode='lines+text',
            line=dict(color='black', width=4),
            text=['', label],
            textposition='top center',
            textfont=dict(size=14, color='black'),
            showlegend=False,
            hoverinfo='skip',
        ),
        # Cone as arrowhead
        go.Cone(
            x=[x], y=[y], z=[z],
            u=[0.05*x], v=[0.05*y], w=[0.05*z],
            colorscale=[[0, 'black'], [1, 'black']],
            showscale=False,
            sizemode='absolute',
            sizeref=3,
            hoverinfo='skip',
        )
    ]

arrow_len = 15
axes_traces = (
    axis_arrow((arrow_len, 0, 0), 'x') +
    axis_arrow((0, arrow_len, 0), 'y') +
    axis_arrow((0, 0, arrow_len), 'z')
)

no_axis = dict(
    showgrid=False, zeroline=False, showline=False,
    showticklabels=False, showaxeslabels=False, visible=False
)

# Compute axis ranges that work for any shift direction.
# Each shifted axis spans from -padding past the near site to +padding past the far site;
# unshifted axes are centred on their midpoint with the same total span.
padding = 15
total_span = lattice_const + 2 * padding

def axis_range(lo, hi):
    span = hi - lo
    if span < total_span:
        c = (lo + hi) / 2
        return [c - total_span / 2, c + total_span / 2]
    return [lo - padding, hi + padding]

x_range = axis_range(min(site_1[0], site_2[0]), max(site_1[0], site_2[0]))
y_range = axis_range(min(site_1[1], site_2[1]), max(site_1[1], site_2[1]))
z_range = axis_range(min(site_1[2], site_2[2]), max(site_1[2], site_2[2]))

fig = go.Figure(data=[
    make_isosurface(V1, f"{orbital_1} at {site_1}", show_scale=True),
    make_isosurface(V2, f"{orbital_2} at {site_2}"),
] + axes_traces)
fig.update_layout(
    title=f"{orbital_1} at {site_1}  |  {orbital_2} at {site_2}",
    scene=dict(
        xaxis=dict(**no_axis, range=x_range),
        yaxis=dict(**no_axis, range=y_range),
        zaxis=dict(**no_axis, range=z_range),
        bgcolor='white',
        aspectmode='cube',
        camera=dict(eye=dict(x=0.7, y=0.7, z=0.7)),
    )
)
shift_label = ''.join(
    f"{ax}{'p' if val > 0 else 'm'}"
    for ax, val in zip('xyz', shift) if val != 0.0
) + 'shift'
filename = f"{shift_label}_{orbital_1}_{orbital_2}.png"

images_dir = os.path.join(os.path.dirname(__file__), 'images')
os.makedirs(images_dir, exist_ok=True)

fig.write_image(os.path.join(images_dir, filename))
print(f"Saved to {os.path.join(images_dir, filename)}")
