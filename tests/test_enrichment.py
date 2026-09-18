import pytest

from app.enrichment import (
    PeakListError,
    abundance_grid,
    enriched_isotopes,
    fit_enrichment,
    validate_observed,
)
from app.envelope import isotope_envelope
from app.isotopes import ISOTOPES
from app.parser import LabelError, parse_formula, parse_label

TOL = 0.0001
PRUNE = 1e-12


def make_observed(formula, symbol, mass_number, fraction, top=8):
    """Synthesise observed peaks (base peak = 100) at a known enrichment."""
    table = dict(ISOTOPES)
    table[symbol] = enriched_isotopes(symbol, mass_number, fraction)
    env = isotope_envelope(parse_formula(formula), TOL, PRUNE, table)
    base = max(i for _, i in env)
    strongest = sorted(env, key=lambda p: p[1], reverse=True)[:top]
    return [(m, 100.0 * i / base) for m, i in sorted(strongest)]


def test_parse_label_ok():
    assert parse_label("13C") == ("C", 13)
    assert parse_label("15N") == ("N", 15)
    assert parse_label("2H") == ("H", 2)


def test_parse_label_malformed():
    with pytest.raises(LabelError) as excinfo:
        parse_label("C13")
    assert excinfo.value.position == 0
    with pytest.raises(LabelError):
        parse_label("")


def test_parse_label_unknown_element():
    with pytest.raises(LabelError) as excinfo:
        parse_label("13Xx")
    assert excinfo.value.position == 2
    assert "Xx" in excinfo.value.message


def test_parse_label_unknown_isotope():
    with pytest.raises(LabelError) as excinfo:
        parse_label("20C")
    assert "20" in excinfo.value.message


def test_abundance_grid_endpoints():
    grid = abundance_grid(0.5)
    assert len(grid) == 201
    assert grid[0] == 0.0 and grid[-1] == 100.0
    assert abundance_grid(10.0) == [i * 10.0 for i in range(11)]
    grid = abundance_grid(3.0)
    assert grid[-1] == 100.0 and len(grid) == 35


def test_enriched_isotopes_bounds():
    assert enriched_isotopes("C", 13, 0.0) == ISOTOPES["C"]
    # Fully enriched: only the pure label isotope remains.
    assert enriched_isotopes("C", 13, 1.0) == ((13.00335483507, 1.0),)
    mid = enriched_isotopes("C", 13, 0.3)
    assert sum(a for _, a in mid) == pytest.approx(1.0)
    assert mid[1][1] == pytest.approx(0.0107 * 0.7 + 0.3)
    assert mid[0][1] == pytest.approx(0.9893 * 0.7)


def test_validate_observed_empty():
    with pytest.raises(PeakListError) as excinfo:
        validate_observed([])
    assert excinfo.value.index is None


def test_validate_observed_non_increasing():
    with pytest.raises(PeakListError) as excinfo:
        validate_observed([(100.0, 1.0), (99.0, 1.0)])
    assert excinfo.value.index == 1


def test_validate_observed_duplicate():
    with pytest.raises(PeakListError) as excinfo:
        validate_observed([(100.0, 1.0), (100.0, 1.0)])
    assert excinfo.value.index == 1
    assert "Duplicate" in excinfo.value.message


def test_validate_observed_all_zero_intensity():
    with pytest.raises(PeakListError):
        validate_observed([(100.0, 0.0), (101.0, 0.0)])


def test_fit_recovers_enrichment():
    observed = make_observed("C6H12O6", "C", 13, 0.2)
    result = fit_enrichment(
        parse_formula("C6H12O6"), "C", 13, observed, 1.0, 0.01, TOL, PRUNE
    )
    assert result["best_abundance_percent"] == pytest.approx(20.0)
    assert result["mse"] == pytest.approx(0.0, abs=1e-9)
    assert result["search_points"] == 101
    assert len(result["pairs"]) == len(observed)
    assert all(pair[2] is not None for pair in result["pairs"])


def test_fit_tie_picks_lower_abundance():
    # 13C does not occur in H2O: every candidate scores identically.
    observed = [(18.01056468, 100.0)]
    result = fit_enrichment(
        parse_formula("H2O"), "C", 13, observed, 10.0, 0.01, TOL, PRUNE
    )
    assert result["best_abundance_percent"] == 0.0
    assert result["search_points"] == 11


def test_unmatched_observed_peak_counts_zero():
    observed = [(18.01056468, 100.0), (500.0, 50.0)]
    result = fit_enrichment(
        parse_formula("H2O"), "H", 2, observed, 10.0, 0.01, TOL, PRUNE
    )
    _, _, pred_mass, pred_intensity, mass_error = result["pairs"][1]
    assert pred_mass is None
    assert pred_intensity == 0.0
    assert mass_error is None
    assert result["mse"] > 0
