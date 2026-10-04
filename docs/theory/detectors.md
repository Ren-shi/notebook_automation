# Detector geometry and response

`physim.nuclear.detectors` turns the detectors of a setup file into geometry and response. Geometry covers where
each face is, which directions from the target hit it, and the solid angle and angles it covers. Response covers
what a particle's energy becomes once it is measured. Lengths are in mm, angles in degrees, solid angles in msr.

## Faces

Each detector is a flat face: a rectangle (`width` × `height`, split into `strips_x` × `strips_y`), a disc
(`radius`), or an annulus (`inner_radius` to `outer_radius`, split into `rings` of equal width and `sectors` of
equal angle). It has a centre $\mathbf c$ and a normal $\mathbf n$ pointing back towards the target, unless
`facing` says otherwise. Its in-plane axes are:

- $\mathbf v$: as close to "up" (+y) as the face allows, or +x for a face looking straight up or down;
- $\mathbf u = \mathbf v \times \mathbf n$, horizontal.

Both are turned by `rotation` about $\mathbf n$. Strips along $u$ are counted from $-u$; sectors are counted from
$+u$ towards $+v$.

A straight track from a point $\mathbf s$ (the target, by default) along $\hat{\mathbf d}$ meets the plane at
$t = (\mathbf c - \mathbf s)\cdot\mathbf n / (\hat{\mathbf d}\cdot\mathbf n)$. It hits the detector when
$t > 0$, the face looks back at it ($\hat{\mathbf d}\cdot\mathbf n < 0$), and the local coordinates
$((\mathbf p - \mathbf c)\cdot\mathbf u, (\mathbf p - \mathbf c)\cdot\mathbf v)$ fall inside the shape. The angle of
incidence is $\arccos(-\hat{\mathbf d}\cdot\mathbf n)$.

## Solid angle

$$
\Omega = \int_\text{face} \frac{\hat{\mathbf r}\cdot(-\mathbf n)}{r^2}\, dA
$$

is integrated with Gauss–Legendre quadrature, 48 × 48 points per face or per segment. Rectangles use $(u, v)$;
discs and annuli use polar coordinates. The integrand is smooth, so this is exact to round-off for every placement,
including off-centre and tilted faces where no closed form exists. Two checks are kept:

- **Closed forms where they exist:** $2\pi(1 - \cos\alpha)$ for a disc on its axis, $4\arctan\!\big(ab /
  (d\sqrt{a^2 + b^2 + d^2})\big)$ for a centred rectangle of half-sides $a, b$, and sums and differences of
  corner rectangles for an off-centre one.
- **An independent Monte Carlo count:** isotropic directions inside the cone that just contains the face.

## Angular coverage

The extreme polar angles of a face lie on its edge, unless the beam axis crosses the face, which gives 0° or 180°.
They are found by sampling the edge densely and refining the best point with a golden-section search, to about
1e-9°. Azimuths are reported unwrapped around the detector's own φ, so a detector at φ = 180° spans, say, 168° to
192°. A face around the beam axis covers all φ. `mean_theta` is the solid-angle-weighted mean polar angle.

## Layout checks

- **Shadowing:** for every detector, rays from the target to its quadrature points are tested against every other
  face. The hidden fraction is weighted by solid angle.
- **Beam blocking:** a face crossed by the beam axis downstream of the target is in the unscattered beam. One
  crossed upstream blocks the incoming beam; an annulus with a central hole does not.

## Response

For a particle of energy $E$ arriving at incidence $\alpha$:

1. It crosses the dead layer, along a path $t_\text{dead}/\cos\alpha$. The energy loss comes from
   {doc}`stopping`.
2. In the active thickness it deposits everything if it stops. If it punches through, it deposits
   $E_\text{in} - E_\text{out}$. The punch-through energy is the one whose range equals
   $(t_\text{dead} + t)/\cos\alpha$.
3. A Gaussian resolution of the given FWHM ($\sigma = \text{FWHM}/2.355$) is added when a random generator is
   passed. Energies below the threshold are dropped (NaN).

## Leaving the target

A reaction at depth $x$ in a target of thickness $t$, whose normal is tilted by $\theta_t$ about $y$, sends a track
$\hat{\mathbf d}$ through the rest of the target:

- forward ($\hat{\mathbf d}\cdot\mathbf n_t > 0$): a path $(t - x)/\hat{\mathbf d}\cdot\mathbf n_t$;
- backward: $x/|\hat{\mathbf d}\cdot\mathbf n_t|$.

Close to the target plane this grows as $1/\cos$, which is the sharp dip LISE++ shows at 90° in its "energy at
detector entrance" plots. `exit_path` caps the path at 1000 thicknesses (a real target is finite) and reports that
it did.
