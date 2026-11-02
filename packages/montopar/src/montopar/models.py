"""Data models for montopar."""
from dataclasses import dataclass, field


@dataclass
class Vertex:
    """An eikonal vertex on the Wilson loop contour."""
    line: int       # edge number (1-indexed)
    lorentz: str    # Lorentz index label (e.g. "mu1")
    color: str      # color index label (e.g. "a1")


@dataclass
class Propagator:
    """A gluon propagator connecting two eikonal vertices."""
    loop: str               # loop momentum label (e.g. "l1")
    lorentz: tuple[str, str]  # (lorentz_a, lorentz_b)
    color: tuple[str, str]    # (color_a, color_b)


@dataclass
class Bound:
    """Integration bound for a Schwinger parameter."""
    var: str    # variable name (e.g. "s31")
    low: str    # lower bound (e.g. "0")
    high: str   # upper bound (e.g. "s32" or "1")


@dataclass
class Pairing:
    """A gluon propagator's edge pairing with associated Schwinger variables."""
    edge_a: int     # first edge number
    edge_b: int     # second edge number
    loop: str       # loop momentum label
    s_a: str        # Schwinger variable on edge_a
    s_b: str        # Schwinger variable on edge_b


@dataclass
class EdgeDistribution:
    """Distribution of attachment points across Wilson loop edges."""
    counts: dict[int, int]          # {edge: num_attachment_points}
    s_vars: dict[int, list[str]]    # {edge: [s_var_names]}, ordered by nesting


@dataclass
class TGVLeg:
    """A propagator leg connected to a TGV."""
    edge: int       # which Wilson-line edge (1-4)
    s_var: str      # Schwinger parameter name (e.g. "s31")
    loop: str       # loop momentum label (e.g. "l1")


@dataclass
class FGVLeg:
    """A propagator leg connected to an FGV."""
    edge: int
    s_var: str
    loop: str


@dataclass
class DoubleTGVLeg:
    """A propagator leg connected to a double TGV (4e-3g+3g)."""
    edge: int
    s_var: str
    loop: str
    tgv_id: str = ""  # "A" or "B" to indicate which TGV


@dataclass
class Orbit:
    """Parsed orbit definition ready for C×J computation."""
    orbit_id: str
    emission_type: str
    n_gluons: int
    edge_distribution: EdgeDistribution
    pairings: list[Pairing]
    bounds: list[Bound]
    vertices: list[Vertex]
    propagators: list[Propagator]
    source: dict = field(default_factory=dict)

    # --- TGV / FGV dispatch fields ---
    dispatch: str = "pure"
    tgv_legs: list[TGVLeg] = field(default_factory=list)
    fgv_legs: list[FGVLeg] = field(default_factory=list)
    double_tgv_legs: list[DoubleTGVLeg] = field(default_factory=list)


@dataclass
class HypercubeResult:
    """Result of hypercube transformation."""
    s_to_x: dict[str, str]    # {"s31": "x1*x3", "s32": "x1", ...}
    jacobian: str             # "x1*x2"
    x_vars: list[str]         # ["x1", "x2", "x3", "x4"]


@dataclass
class CFactors:
    """C_n prefactor decomposed into 5 factors."""
    color: str
    coupling: str
    ft_base: str
    v14_scale: str
    twoscale: str
    derivative: str | None = None   # TGV ⑤ derivative operator (∂/∂v14)


@dataclass
class JFactors:
    """J_n parametric integrand components."""
    lorentz_z: str              # dimensionless Lorentz factor (e.g. "z1^2")
    jacobian: str               # hypercube Jacobian (e.g. "x1*x2")
    f_product: str              # product of f_p^(-1+eps) terms
    f_polynomials: dict[str, str]  # {"f1": "...", "f2": "..."}
    lorentz_v14_power: str      # power of v14 from Lorentz scalar (e.g. "2")


@dataclass
class TGVJFactors:
    """J integrand for TGV diagrams."""
    jacobian: str               # hypercube Jacobian (e.g. "x2")
    pnum: str                   # numerator polynomial N (含 α_i, x_i, z)
    alpha_monomial: str         # "(alpha1*alpha2*alpha3)^(-eps)"
    S_polynomial: str           # denominator polynomial S (含 α_i, x_i, z)
    denom_power: str            # "2 - 2*eps"
    x_vars: list[str]           # ["x1", "x2", "x3"]


@dataclass
class TGVIBPTerm:
    """One independently integrable term in a TGV contour-IBP sum."""
    term_id: str
    kind: str
    coefficient: str
    c_ref: str
    integrand: str
    variables: list[str]
    alpha_vars: list[str]
    bounds: dict[str, list[int | str]]
    delta_constraints: dict
    denom_power: str
    source_derivative: str | None = None
    fixed_parameters: dict[str, str] = field(default_factory=dict)
    identified_parameters: dict[str, str] = field(default_factory=dict)
    s_to_x: dict[str, str] = field(default_factory=dict)
    combined_C: str = ""


@dataclass
class TGVIBPResult:
    """Contour-IBP representation of one TGV orbit variant."""
    method: str
    status: str
    case: str
    C_star: CFactors
    terms: list[TGVIBPTerm]
    sum_expression: str
    certificates: dict
    derivative_vector: list[str] = field(default_factory=list)
    reason: str | None = None


@dataclass
class FGVJFactors:
    """J integrand for FGV diagrams."""
    jacobian: str               # hypercube Jacobian
    lorentzZ: str               # pure z polynomial (e.g. "-(1+4*z1+z1^2)/4")
    alpha_monomial: str         # "(alpha1*alpha2*alpha3*alpha4)^(-eps)"
    Q_prime_polynomial: str     # denominator polynomial Q' (含 α_i, x_i, z)
    denom_power: str            # "2 - 3*eps"
    x_vars: list[str]           # ["x1", "x2", "x3", "x4"]
    physical_zero: bool = False # True if max(m_i) ≥ 3


@dataclass
class MixedJFactors:
    """J integrand for mixed (TGV + pure) diagrams."""
    # Pure part
    lorentz_z: str              # dimensionless Lorentz factor (e.g. "z1")
    jacobian: str               # hypercube Jacobian
    f_product: str              # product of f_p^(-1+eps) terms
    f_polynomials: dict[str, str]  # {"f1": "...", "f2": "..."}
    # TGV part
    pnum: str                   # numerator polynomial N
    alpha_monomial: str         # "(alpha1*alpha2*alpha3)^(-eps)"
    S_polynomial: str           # denominator polynomial S
    denom_power: str            # "2 - 2*eps"

    # Shared
    x_vars: list[str]           # all x-variables (3+2n)
    n_pure: int                 # number of pure gluon pairings
    edge_pairs: list[tuple[int, int]] = field(default_factory=list)
                                # pure gluon edge pairs in canonical frame


@dataclass
class DoubleTGVJFactors:
    """J integrand for double TGV (4e-3g+3g) diagrams."""
    jacobian: str               # hypercube Jacobian
    B5_polynomial: str           # denominator polynomial B5 (含 α_i, x_i, z)
    U_polynomial: str            # first Symanzik polynomial U
    pnum1: str                   # Pnum^(1) numerator (含 α_i, x_i, z, ε)
    pnum2: str                   # Pnum^(2) numerator (含 α_i, x_i, z, ε)
    alpha_monomial: str          # "(alpha1*alpha2*alpha3*alpha4*alpha5)^(-eps)"
    denom_power_1: str           # "3 - 3*eps"
    denom_power_2: str           # "2 - 3*eps"
    x_vars: list[str]            # x variable names
    physical_zero: bool = False  # True if color = 0 (U-type)
