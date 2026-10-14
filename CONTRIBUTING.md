# 🤝 Contributing to geo3d

Thanks for helping make geo3d better! This guide covers everything you need, from your first issue
to your first pull request.

## 💡 Ways to help

- 🐛 **Report a bug.** Open an issue with the input you used, what you expected, and what you got. A
  short code snippet is perfect.
- ✨ **Ask for a feature.** Check the [roadmap](https://github.com/FusionVault/geo3d/issues) first.
  If it's already there, a 👍 or a comment about your use case helps us prioritise.
- 🧑‍💻 **Write code.** Pick an issue and say you're working on it. Issues marked
  [good first issue](https://github.com/FusionVault/geo3d/issues?q=is%3Aopen+label%3A%22good+first+issue%22)
  are small and well described.
- 📖 **Improve the docs.** If something confused you, it will confuse others too. Fixes to wording
  are always welcome.

## 🛠️ Get set up

You need Rust 1.85 or newer.

```bash
git clone https://github.com/FusionVault/geo3d
cd geo3d
cargo test --all-features
```

That runs the unit tests, the cross-checks against other libraries, and every example in the
README.

## ✅ Before you open a pull request

Run the same checks as CI:

```bash
cargo fmt --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo test --no-default-features --features libm
cargo doc --no-deps --all-features
```

CI also checks that you haven't accidentally broken the public API, by comparing it with the last
release on crates.io. When you open the pull request, the template has a short checklist. Fill it
in, and you're done.

## 📏 Our three rules

**1. Accuracy comes first.** If there's a choice between faster and more accurate, pick accurate.
New maths needs a test against a trusted reference, such as GeographicLib, a published value, or
an independent crate. Its doc comment should say how accurate it is.

**2. Docs ship with the code.** Every change adds a line to [CHANGELOG.md](CHANGELOG.md), and CI
checks for it. New features also get a README recipe. Every public item gets a doc comment that
gives its units and has an example.

**3. Keep it light.** No new required dependencies, and no `unsafe` code. Integrations with other
crates go behind an optional feature. The library must also build without the standard library:
use `core::`, never `std::`, and don't allocate.

## 🎨 Style

- The public API uses degrees and metres. Names end in their unit: `_deg`, `_m`, `_mps`.
- Latitude comes before longitude, always.
- Keep doc comments short and plain. Write for someone who has never heard of ECEF.
- `cargo fmt` sets the layout (lines up to 120 characters).

## 📊 The comparison suite

The [`comparisons/`](comparisons/README.md) folder tests geo3d against other crates. If your change
affects accuracy or speed, re-run it and update the charts:

```bash
cd comparisons
cargo run --release --bin accuracy > results/accuracy.md
cargo bench
python3 charts.py
```

## 🚀 Releasing (maintainers)

1. Move the "Unreleased" notes in [CHANGELOG.md](CHANGELOG.md) under the new version and date.
2. Bump `version` in `Cargo.toml`.
3. Merge, then tag the commit: `git tag vX.Y.Z && git push --tags`.
4. Run `cargo publish`.
5. Create a GitHub release from the tag, with the changelog notes.

## 💬 Be kind

Be patient, constructive and respectful. We follow the
[Rust Code of Conduct](https://www.rust-lang.org/policies/code-of-conduct).
