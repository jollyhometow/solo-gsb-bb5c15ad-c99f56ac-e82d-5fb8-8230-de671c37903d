import pytest
from fastapi.testclient import TestClient

from app.config import settings


@pytest.fixture()
def client(tmp_path, monkeypatch):
    # Isolate each test run with its own SQLite file.
    monkeypatch.setattr(settings, "db_path", str(tmp_path / "test.db"))
    from app.main import app

    with TestClient(app) as c:
        yield c


def test_health(client):
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_envelope_glucose(client):
    resp = client.post("/api/v1/envelope", json={"formula": "C6H12O6", "top_n": 10})
    assert resp.status_code == 200
    body = resp.json()
    assert body["monoisotopic_mass"] == pytest.approx(180.063388, abs=1e-5)
    assert body["peak_count"] == 10
    assert len(body["peaks"]) == 10
    assert max(p["intensity"] for p in body["peaks"]) == pytest.approx(100.0)
    assert body["cached"] is False
    masses = [p["mass"] for p in body["peaks"]]
    assert masses == sorted(masses)


def test_envelope_with_label(client):
    resp = client.post("/api/v1/envelope", json={"formula": "C6H6[13C]2", "top_n": 5})
    assert resp.status_code == 200
    body = resp.json()
    assert body["labels"] == {"[13C]": 2}
    assert body["monoisotopic_mass"] == pytest.approx(
        6 * 12.0 + 6 * 1.00782503223 + 2 * 13.00335483507, abs=1e-5
    )


def test_envelope_cached_on_second_call(client):
    payload = {"formula": "C2H6O", "top_n": 5}
    first = client.post("/api/v1/envelope", json=payload).json()
    second = client.post("/api/v1/envelope", json=payload).json()
    assert first["cached"] is False
    assert second["cached"] is True
    assert first["peaks"] == second["peaks"]


def test_invalid_formula_reports_position(client):
    resp = client.post("/api/v1/envelope", json={"formula": "C6Xx2"})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_formula"
    assert detail["position"] == 2
    assert "Xx" in detail["message"]


def test_invalid_label_mass_reports_position(client):
    resp = client.post("/api/v1/envelope", json={"formula": "C6[20C]"})
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["position"] == 2
    assert "20" in detail["message"]


def test_top_n_bounds(client):
    resp = client.post("/api/v1/envelope", json={"formula": "H2O", "top_n": 0})
    assert resp.status_code == 422
    resp = client.post(
        "/api/v1/envelope", json={"formula": "H2O", "top_n": settings.top_n_max + 1}
    )
    assert resp.status_code == 422


def test_missing_formula(client):
    resp = client.post("/api/v1/envelope", json={})
    assert resp.status_code == 422


def test_elements_endpoint(client):
    resp = client.get("/api/v1/elements")
    assert resp.status_code == 200
    body = resp.json()
    assert "C" in body and "O" in body
    assert body["C"][0]["mass"] == pytest.approx(12.0)


def _observed_peaks(formula, symbol, mass_number, fraction, top=8):
    """Synthesise observed peaks at a known enrichment for fit tests."""
    from app.enrichment import enriched_isotopes
    from app.envelope import isotope_envelope
    from app.isotopes import ISOTOPES
    from app.parser import parse_formula

    table = dict(ISOTOPES)
    table[symbol] = enriched_isotopes(symbol, mass_number, fraction)
    env = isotope_envelope(parse_formula(formula), 0.0001, 1e-12, table)
    base = max(i for _, i in env)
    strongest = sorted(env, key=lambda p: p[1], reverse=True)[:top]
    return [{"mass": m, "intensity": 100.0 * i / base} for m, i in sorted(strongest)]


def test_enrichment_fit_recovers_abundance(client):
    peaks = _observed_peaks("C6H12O6", "C", 13, 0.2)
    resp = client.post(
        "/api/v1/enrichment-fit",
        json={
            "formula": "C6H12O6",
            "label": "13C",
            "peaks": peaks,
            "step": 1.0,
            "tolerance": 0.01,
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["formula"] == "C6H12O6"
    assert body["label"] == "13C"
    assert body["best_abundance_percent"] == pytest.approx(20.0)
    assert body["mse"] == pytest.approx(0.0, abs=1e-6)
    assert body["search_points"] == 101
    assert len(body["pairs"]) == len(peaks)
    assert all(p["predicted_mass"] is not None for p in body["pairs"])


def test_enrichment_fit_unmatched_peak(client):
    peaks = [
        {"mass": 18.010565, "intensity": 100.0},
        {"mass": 500.0, "intensity": 80.0},
    ]
    resp = client.post(
        "/api/v1/enrichment-fit",
        json={
            "formula": "H2O",
            "label": "2H",
            "peaks": peaks,
            "step": 10.0,
            "tolerance": 0.01,
        },
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["search_points"] == 11
    unmatched = body["pairs"][1]
    assert unmatched["predicted_mass"] is None
    assert unmatched["predicted_intensity"] == 0.0
    assert body["mse"] > 0


def test_enrichment_fit_tie_prefers_lower_abundance(client):
    # 13C does not occur in H2O, so all candidates score identically.
    resp = client.post(
        "/api/v1/enrichment-fit",
        json={
            "formula": "H2O",
            "label": "13C",
            "peaks": [{"mass": 18.010565, "intensity": 100.0}],
            "step": 10.0,
            "tolerance": 0.01,
        },
    )
    assert resp.status_code == 200
    assert resp.json()["best_abundance_percent"] == 0.0


def test_enrichment_fit_invalid_label(client):
    for label in ("C13", "13Xx", "20C", ""):
        resp = client.post(
            "/api/v1/enrichment-fit",
            json={
                "formula": "H2O",
                "label": label,
                "peaks": [{"mass": 18.0, "intensity": 1.0}],
            },
        )
        assert resp.status_code == 422, label
        detail = resp.json()["detail"]
        assert detail["error"] == "invalid_label"
        assert "position" in detail


def test_enrichment_fit_empty_peaks(client):
    resp = client.post(
        "/api/v1/enrichment-fit",
        json={"formula": "H2O", "label": "13C", "peaks": []},
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "invalid_observed_peaks"


def test_enrichment_fit_non_increasing_masses(client):
    resp = client.post(
        "/api/v1/enrichment-fit",
        json={
            "formula": "H2O",
            "label": "13C",
            "peaks": [
                {"mass": 20.0, "intensity": 5.0},
                {"mass": 18.0, "intensity": 100.0},
            ],
        },
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_observed_peaks"
    assert detail["index"] == 1


def test_enrichment_fit_duplicate_masses(client):
    resp = client.post(
        "/api/v1/enrichment-fit",
        json={
            "formula": "H2O",
            "label": "13C",
            "peaks": [
                {"mass": 18.0, "intensity": 100.0},
                {"mass": 18.0, "intensity": 50.0},
            ],
        },
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert detail["error"] == "invalid_observed_peaks"
    assert detail["index"] == 1


def test_enrichment_fit_too_many_peaks(client):
    peaks = [{"mass": float(i + 1), "intensity": 1.0} for i in range(201)]
    resp = client.post(
        "/api/v1/enrichment-fit",
        json={"formula": "H2O", "label": "13C", "peaks": peaks},
    )
    assert resp.status_code == 422


def test_enrichment_fit_step_and_tolerance_bounds(client):
    base = {
        "formula": "H2O",
        "label": "13C",
        "peaks": [{"mass": 18.0, "intensity": 1.0}],
    }
    for bad in ({"step": 0.1}, {"step": 10.5}, {"tolerance": 0}, {"tolerance": -0.5}):
        resp = client.post("/api/v1/enrichment-fit", json={**base, **bad})
        assert resp.status_code == 422, bad


def test_enrichment_fit_invalid_formula(client):
    resp = client.post(
        "/api/v1/enrichment-fit",
        json={
            "formula": "Xx2",
            "label": "13C",
            "peaks": [{"mass": 1.0, "intensity": 1.0}],
        },
    )
    assert resp.status_code == 422
    assert resp.json()["detail"]["error"] == "invalid_formula"
