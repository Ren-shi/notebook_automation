"""The Analysis tab: from the peak to B(E2), the alignment check, and the multi-step solution."""


def render(ctx) -> None:
    exp = ctx.P().experiment
    if exp.excitation is None or not exp.gamma_detectors:
        ctx.ui.markdown("The analysis needs Coulomb excitation with γ-ray detectors: set the reaction in the "
                        "**Target and reaction** card and add γ-ray detectors.")
    else:
        ctx.analysis_block()
        ctx.alignment_block()
    ctx.multistep_block()
