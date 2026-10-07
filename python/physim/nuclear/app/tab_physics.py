"""The Physics tab: the physics behind the numbers, as a reference: kinematics, energy loss, trajectories, the
excitation."""


def render(ctx) -> None:
    ctx.block("kinematics", "Kinematics")
    ctx.block("energy_loss", "Energy loss")
    ctx.block("trajectories", "Trajectories")
    if ctx.P().experiment.excitation is not None:
        ctx.block("gamma", "Excitation and γ rays")
