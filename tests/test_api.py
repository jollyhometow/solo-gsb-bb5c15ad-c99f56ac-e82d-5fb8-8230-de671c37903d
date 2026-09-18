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
