## What changed

<!-- One or two sentences a user would understand. -->

## Docs checklist

Docs ship with the code, in the same pull request:

- [ ] `CHANGELOG.md` has a line for this change (CI fails without one when `src/` changes)
- [ ] README updated where users will see it: a new feature gets a recipe, a changed accuracy updates the accuracy table, a new MATLAB/pymap3d equivalent goes in the lookup table
- [ ] Every new public item has a doc comment with units, output ranges and an example
- [ ] The crate docs' "find what you need" table points at anything new
- [ ] `cargo test --all-features` passes (the README examples are doctests)
