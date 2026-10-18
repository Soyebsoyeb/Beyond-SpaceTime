"""Central numerical conventions and physical constants.

Signature (-,+,+,+); Riemann R^a_{bcd} = d_c Gamma^a_{db} - d_d Gamma^a_{cb}
+ Gamma^a_{ce} Gamma^e_{db} - Gamma^a_{de} Gamma^e_{cb}; Ricci R_{bd} = R^a_{bad}.
"""

# SI constants (CODATA 2018 values; exact where defined)
C_SI = 299_792_458.0            # m / s
G_SI = 6.67430e-11              # m^3 kg^-1 s^-2
M_SUN_KG = 1.98841e30           # kg (nominal solar mass)
G_EARTH = 9.80665               # m / s^2

# Numerical tolerances used by tests and experiments
TOL_SYMBOLIC = 0.0              # symbolic identities must simplify to exactly zero
TOL_ANALYTIC = 1e-9             # numerical vs closed-form checks (relative)
TOL_ODE = 1e-10                 # integrator rtol used for geodesics
TOL_QUANTUM = 1e-10             # density-matrix validity checks


def geometric_mass_m(mass_kg: float) -> float:
    """Gravitational radius GM/c^2 in metres."""
    return G_SI * mass_kg / C_SI**2


def tidal_accel_si(M_geom_m: float, r_m: float, length_m: float, factor: float = 2.0) -> float:
    """Radial tidal acceleration (m/s^2) across ``length_m`` at areal radius ``r_m``.

    In geometric units the radial tidal component of a radial geodesic is
    ``factor * M / r^3`` per unit separation (factor = 2); multiply by c^2 for SI.
    """
    return C_SI**2 * factor * M_geom_m / r_m**3 * length_m
