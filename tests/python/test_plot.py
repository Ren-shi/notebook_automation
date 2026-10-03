import numpy as np
import pytest

import physim as ps

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
plt = pytest.importorskip("matplotlib.pyplot")


@pytest.fixture
def traj():
    w = ps.scenarios.figure_eight()
    return w.run(0.01, 300, record_every=10)


def test_static_plots(traj):
    ax = ps.plot.orbits(traj, labels=["a", "b", "c"])
    assert len(ax.lines) == 6  # a path and a start dot per particle
    assert np.allclose(ax.lines[0].get_xydata(), traj.pos[:, 0, :2])
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["a", "b", "c"]
    assert not ax.spines["top"].get_visible()

    ax = ps.plot.orbits(traj, particles=[2], plane="xz", start=False)
    assert len(ax.lines) == 1 and ax.get_ylabel() == "z"

    ax = ps.plot.energy_error({"yoshida4": traj, "again": traj})
    assert ax.get_yscale() == "log" and len(ax.lines) == 2
    y = ax.lines[0].get_ydata()
    assert np.allclose(y, np.abs(ps.relative_energy_error(traj)) + 1e-17)

    ax = ps.plot.phase_space(traj, particle=1, axis="y", momentum=True, mass=2.0)
    assert np.allclose(ax.lines[0].get_ydata(), 2.0 * traj.vel[:, 1, 1])

    # Rigid trajectories work too.
    s = ps.RigidSystem()
    s.add_body(1.0, [1, 1, 1], velocity=[1, 0, 0])
    ps.plot.orbits(s.run(0.1, 10))
    plt.close("all")


def test_bad_arguments(traj):
    with pytest.raises(ValueError, match="plane"):
        ps.plot.orbits(traj, plane="xx")
    with pytest.raises(ValueError, match="out of range"):
        ps.plot.orbits(traj, particles=[5])
    with pytest.raises(ValueError, match="mass"):
        ps.plot.phase_space(traj, momentum=True)
    no_energy = ps.scenarios.figure_eight().run(0.01, 10, energies=False)
    with pytest.raises(ValueError, match="energies"):
        ps.plot.energy_error(no_energy)
    plt.close("all")


def test_animation(traj, tmp_path):
    anim = ps.plot.animate(traj, trail=5, every=2, labels=["a", "b", "c"])
    frames = list(anim.new_frame_seq())
    assert len(frames) == len(traj.t[::2])
    artists = anim._func(7)
    dots = artists[0]
    assert np.allclose(dots.get_offsets(), traj.pos[14, :, :2])
    tail = artists[1].get_xydata()
    assert np.allclose(tail, traj.pos[2 * 2:2 * 7 + 1:2, 0, :2])
    out = tmp_path / "eight.gif"
    anim.save(out, writer="pillow", fps=20)
    assert out.stat().st_size > 1000
    plt.close("all")


def test_view3d(traj):
    pytest.importorskip("plotly")
    fig = ps.plot.view3d(traj, labels=["a", "b", "c"], every=3)
    assert len(fig.data) == 4 and len(fig.frames) == len(traj.t[::3])
    assert np.allclose(fig.frames[2].data[0].x, traj.pos[6, :, 0])


def test_wavefunction_plots(tmp_path):
    s = ps.Schrodinger([128], 0.2, potential=lambda t, x: 0.5 * x**2)
    s.set_gaussian(center=2.0, width=0.7, momentum=1.0)
    ax = ps.plot.wavefunction(s)
    assert np.allclose(ax.lines[0].get_ydata(), s.density)
    assert len(ax.collections) == 1  # phase-coloured fill
    assert len(ax.figure.axes) == 2  # potential on a twin axis
    ax = ps.plot.wavefunction(s.psi, s.axes[0], phase=False, potential=False)
    assert len(ax.figure.axes) == 1

    s2 = ps.Schrodinger([32, 16], 0.25)
    s2.set_gaussian(width=0.6, momentum=[1.0, 0.0])
    ax = ps.plot.wavefunction(s2)
    img = ax.images[0].get_array()
    assert img.shape == (16, 32, 3)
    dim = ps.plot.wavefunction(s2, gamma=0.5).images[0].get_array()
    assert dim.max() == pytest.approx(img.max()) and dim.mean() > img.mean()
    with pytest.raises(ValueError, match="1D or 2D"):
        ps.plot.wavefunction(np.zeros((4, 4, 4), complex))

    t, frames = s.run(0.05, 40, record_every=4)
    anim = ps.plot.animate_wavefunction(t, frames, s.axes[0], potential=s.potential, every=2)
    assert len(list(anim.new_frame_seq())) == len(t[::2])
    fill, line = anim._func(3)[:2]
    assert np.allclose(line.get_ydata(), np.abs(frames[6]) ** 2)
    out = tmp_path / "packet.gif"
    anim.save(out, writer="pillow", fps=10)
    assert out.stat().st_size > 1000
    with pytest.raises(ValueError, match="shape"):
        ps.plot.animate_wavefunction(t, frames, s.axes[0][:10])
    plt.close("all")
