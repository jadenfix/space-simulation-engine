# 3D3V midpoint nonlinear-residual conditioning

## Scope

The bounded 3D3V electromagnetic PIC reference advances particles and fields with a fixed-point iteration for the midpoint electric field. The iteration must remain fail-closed: a state is committed only when the declared nonlinear tolerance and every charge, field, work, momentum, Courant, and finiteness gate pass.

The multi-period homogeneous-plasma benchmark exposed a numerical-conditioning problem in the convergence diagnostic near an electric-field zero crossing. The particle and field update remained finite, charge conserving, Gauss consistent, subluminal, and tightly closed in energy and momentum, but the old *relative* fixed-point residual could exceed its tolerance because its denominator approached zero.

## Previous normalization

The former diagnostic used

```text
||E_mid^(k+1) - E_mid^k||_2
---------------------------------
max(||E_mid^(k+1)||_2, ||E_mid^k||_2)
```

This is well conditioned while the midpoint field has appreciable magnitude. It is not well conditioned when the physical oscillation passes through `E_mid = 0`: both midpoint norms can become much smaller than the resolved endpoint field change over the timestep. A tiny absolute iteration defect can then appear as a large relative defect.

Increasing the iteration cap did not change the observed residual floor, confirming that this was a normalization problem rather than an insufficient iteration budget.

## Conditioned normalization

The accepted diagnostic retains the same numerator and nonlinear tolerance, but uses the larger of the two midpoint norms and half of the resolved endpoint field jump:

```text
D = ||E_mid^(k+1) - E_mid^k||_2

S = max(
      ||E_mid^k||_2,
      ||E_mid^(k+1)||_2,
      0.5 ||E^(n+1) - E^n||_2,
      DBL_MIN
    )

nonlinear_relative_residual = D / S
```

The endpoint-jump term is not a tolerance relaxation. It supplies a physically resolved field scale when the implicit midpoint itself crosses zero. Away from zero crossings, one of the midpoint norms remains the dominant scale and the diagnostic reduces to the previous behavior.

## Acceptance properties

The repair preserves the following contracts:

- The configured nonlinear tolerance is unchanged.
- The plasma-mode benchmark's maximum nonlinear-iteration budget is restored to 24; the temporary increase to 32 is not retained.
- A timestep still fails closed when the conditioned residual exceeds the requested tolerance.
- No failed step is committed to the caller-visible state.
- Charge closure, discrete continuity, initial and final Gauss law, magnetic divergence, local-current and spectral-current checks, current means, particle and field work, total energy, total momentum, Courant limits, and subluminal motion remain independent gates.

## Validation boundary

Before publication, the conditioned update was exercised through strict GCC and Clang builds, complete CTest graphs, deterministic receipt replay, standalone Make, the root-engine regression suite, and AddressSanitizer/UndefinedBehaviorSanitizer. The multi-period plasma-mode benchmark was also rerun directly against the published source shape.

These checks establish finite software and numerical behavior for the tested bounded periodic cases. They do not establish that the solver is production PIC, exactly energy conserving for arbitrary states, free of numerical dispersion, or sufficient for chamber or flight claims. Spatial, timestep, particle-count, shape-function, domain-size, open-boundary, collision, emission, ionization, sheath, circuit, thermal, and flexible-structure refinement remain separate fidelity gates.
