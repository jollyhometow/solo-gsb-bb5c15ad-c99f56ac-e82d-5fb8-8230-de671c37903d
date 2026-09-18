"""Isotope-enrichment estimation by abundance grid search.

For a candidate enrichment fraction ``f`` the labelled element's isotope
distribution is a blend of the natural distribution and the pure label
isotope::

    P(iso) = (1 - f) * natural_abundance(iso) + f * [iso == label]

so ``f = 0`` is the natural distribution and ``f = 1`` the fully labelled
element. Each candidate envelope is generated with the existing
polynomial-convolution machinery and normalised to base peak 100; the
observed peaks (likewise normalised) are paired with the nearest predicted
mass within the tolerance, unmatched observed peaks count as predicted
intensity 0, and the abundance minimising the intensity mean squared error
wins. Ties go to the lower abundance.
"""

from __future__ import annotations

from bisect import bisect_left

from .envelope import Peak, isotope_envelope
from .isotopes import ISOTOPES
from .parser import ParsedFormula

ObservedPeak = tuple[float, float]  # (mass Da, intensity)

# A pairing at one candidate abundance:
# (observed_mass, observed_intensity, predicted_mass | None,
#  predicted_intensity, mass_error | None)
Pairing = tuple[float, float, float | None, float, float | None]


class PeakListError(ValueError):
    """Invalid observed-peak list.

    Carries the 0-based ``index`` of the offending peak, or ``None`` when
    the problem concerns the list as a whole (e.g. it is empty).
    """

    def __init__(self, message: str, index: int | None = None) -> None:
        super().__init__(message)
        self.message = message
        self.index = index


def validate_observed(peaks: list[ObservedPeak]) -> None:
    """Validate the observed peaks: non-empty, strictly increasing masses,
    at least one positive intensity. Raises :class:`PeakListError`."""
    if not peaks:
        raise PeakListError("Observed peak list must not be empty")
    prev_mass: float | None = None
    for index, (mass, _) in enumerate(peaks):
        if prev_mass is not None and mass <= prev_mass:
            if mass == prev_mass:
                raise PeakListError(
                    f"Duplicate observed mass {mass} at index {index}; "
                    f"masses must be strictly increasing",
                    index,
                )
            raise PeakListError(
                f"Observed masses must be in increasing order; mass at "
                f"index {index} ({mass}) is below the previous mass "
                f"({prev_mass})",
                index,
            )
        prev_mass = mass
    if max(intensity for _, intensity in peaks) <= 0:
        raise PeakListError(
            "Observed intensities must contain at least one positive value"
        )


def abundance_grid(step: float) -> list[float]:
    """Candidate abundances in percent: 0, step, 2*step, ..., 100.

    Both endpoints are always included, so the final interval may be
    smaller than ``step``.
    """
    points = [i * step for i in range(int(100.0 / step) + 1)]
    points = [p for p in points if p < 100.0]
    points.append(100.0)
    return points


def enriched_isotopes(
    symbol: str, mass_number: int, fraction: float
) -> tuple[tuple[float, float], ...]:
    """Isotope distribution of ``symbol`` enriched to ``fraction`` of the
    pure ``mass_number`` isotope; the remaining ``1 - fraction`` follows
    the natural distribution. Zero-abundance entries (which contribute
    nothing) are dropped."""
    out = []
    for mass, abundance in ISOTOPES[symbol]:
        adjusted = abundance * (1.0 - fraction)
        if round(mass) == mass_number:
            adjusted += fraction
        if adjusted > 0.0:
            out.append((mass, adjusted))
    return tuple(out)


def _normalise(envelope: list[Peak]) -> list[Peak]:
    """Normalise to base peak = 100 (an empty envelope stays empty)."""
    if not envelope:
        return []
    base = max(i for _, i in envelope)
    return [(m, 100.0 * i / base) for m, i in envelope]


def _nearest_index(
    masses: list[float], target: float, tolerance: float
) -> int | None:
    """Index of the nearest mass within ``tolerance`` of ``target``."""
    i = bisect_left(masses, target)
    best_index = None
    best_distance = None
    for j in (i - 1, i):
        if 0 <= j < len(masses):
            distance = abs(masses[j] - target)
            if distance <= tolerance and (
                best_distance is None or distance < best_distance
            ):
                best_index, best_distance = j, distance
    return best_index


def _match(
    observed: list[ObservedPeak], envelope: list[Peak], tolerance: float
) -> tuple[float, list[Pairing]]:
    """Pair each observed peak with the nearest predicted mass within
    ``tolerance`` and return ``(mse, pairs)``.

    Unmatched observed peaks pair with predicted intensity 0. The MSE is
    the mean of squared intensity errors over the observed peaks.
    """
    masses = [m for m, _ in envelope]
    squared_error = 0.0
    pairs: list[Pairing] = []
    for obs_mass, obs_intensity in observed:
        index = _nearest_index(masses, obs_mass, tolerance)
        if index is None:
            pairs.append((obs_mass, obs_intensity, None, 0.0, None))
            squared_error += obs_intensity**2
        else:
            pred_mass, pred_intensity = envelope[index]
            pairs.append(
                (
                    obs_mass,
                    obs_intensity,
                    pred_mass,
                    pred_intensity,
                    pred_mass - obs_mass,
                )
            )
            squared_error += (obs_intensity - pred_intensity) ** 2
    return squared_error / len(observed), pairs


def fit_enrichment(
    parsed: ParsedFormula,
    symbol: str,
    mass_number: int,
    observed: list[ObservedPeak],
    step: float,
    mass_tolerance: float,
    merge_tolerance: float,
    prune: float,
) -> dict:
    """Grid-search the enrichment of the ``mass_number`` isotope of
    ``symbol`` that best reproduces the ``observed`` (mass, intensity)
    peaks.

    ``observed`` must be validated first (see :func:`validate_observed`);
    it is normalised to base peak 100 before comparison, as is every
    candidate envelope. ``step`` is the abundance grid step in percentage
    points. Returns a dict with the best abundance (percent), its
    intensity MSE, the per-peak pairing at the best abundance, and the
    number of searched grid points.
    """
    grid = abundance_grid(step)
    observed_max = max(intensity for _, intensity in observed)
    observed_norm = [(m, 100.0 * i / observed_max) for m, i in observed]

    best_mse: float | None = None
    best_percent = 0.0
    best_pairs: list[Pairing] = []
    for percent in grid:
        table = dict(ISOTOPES)
        table[symbol] = enriched_isotopes(symbol, mass_number, percent / 100.0)
        envelope = _normalise(isotope_envelope(parsed, merge_tolerance, prune, table))
        mse, pairs = _match(observed_norm, envelope, mass_tolerance)
        # Strict "<": on ties the lower abundance (earlier grid point) wins.
        if best_mse is None or mse < best_mse:
            best_mse, best_percent, best_pairs = mse, percent, pairs

    return {
        "best_abundance_percent": best_percent,
        "mse": best_mse,
        "search_points": len(grid),
        "pairs": best_pairs,
    }
