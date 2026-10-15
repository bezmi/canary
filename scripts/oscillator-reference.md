# Oscillator reference generation

`tests/oscillator_reference.hpp` supplies binary64 reference values for the
truncated-oscillator energy test. Regenerate it from the repository root with
Python 3.9 or later; only the standard library is required:

```sh
python3 scripts/generate_oscillator_reference.py --output tests/oscillator_reference.hpp
```

The command replaces the output file after all precision comparisons pass. A
precision disagreement raises `ArithmeticError` without replacing the file.
The generated header is checked into version control, so C++ tests do not need
Python or an extended-precision `long double` implementation.

For vibrational temperature `v`, dissociation temperature `d`, and reciprocal
temperature difference `x`, the generator evaluates the defining expression:

```text
G(x) / k_B = v / (exp(v*x) - 1) - d / (exp(d*x) - 1)
G(0) / k_B = (d - v) / 2
```

Temperatures and `G/k_B` are in kelvin; `x` is in inverse kelvin. The input pairs
come from the N2, O2 and NO chemistry-vibration records in
`include/canary/data/air5_park_exchange.hpp`. Dissociation temperature is the
record's molar dissociation energy divided by `N_A*k_B`, using the same binary64
SI constants and arithmetic as construction. The C++ test checks that these
reference temperatures still match the constructed exchange records.

Each pair retains 45 arguments: zero, both signs of `1e-300`, both signs of
`10^e` for integer `e` from -12 to -2, and both signs of `w/d` for
`w = 0.1, 0.4999, 0.5, 0.5001, 1, 3, 10, 100, 1000, 10000`. This includes the
cancellation limit, values around the evaluator's `|d*x| = 0.5` boundary and
large positive and negative arguments.

Python `Decimal.from_float` preserves each binary64 input exactly. The defining
expression is evaluated independently at 800 and 1000 decimal digits, then
rounded to binary64; both results must agree exactly. The high precision covers
both subtraction in `exp(v*x)-1` and cancellation between terms of order
`1/x`, including `x = 1e-300`. This calculation does not call the Canary evaluator
or use its Bernoulli expansion. Hexadecimal C++ literals preserve the exact input
bits and rounded outputs. The C++ comparison retains its relative tolerance of
`1e-12`; the zero-limit equality check is retained separately.
