"""The Plan tab: what one checks before asking for beam: the key numbers, the rates and beam time, the excitation
and its γ rays, and the energy loss."""


def render(ctx) -> None:
    ctx.readouts()
    ctx.block("rates", "Rates and beam time")
    if ctx.P().experiment.excitation is not None:
        ctx.block("gamma", "Excitation and γ rays")
    ctx.block("energy_loss", "Energy loss")
