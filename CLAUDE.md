# Working in this repository

## Which model does what (set by the user, 2026-11-24)

Work that is split between a parent session and the sessions or agents it starts:

- **The parent session plans, on Fable.** The design, the split into waves, the briefs, and the file ownership of
  each piece.
- **The work itself runs on Opus at high effort.** Every session or agent the parent starts to write code, tests
  or documents gets `model: claude-opus-5-5` and `effort: high`. Fable does not do the implementation.
- **The parent reviews what comes back, on Fable.** It reads each pull request, makes small changes itself where
  that is quicker than sending the work back, resolves conflicts between the pieces, and merges.

Why: three or four Fable contexts running in parallel used up the 5-hour limit in one afternoon. With the work on
Opus there is no need to watch the token usage.
