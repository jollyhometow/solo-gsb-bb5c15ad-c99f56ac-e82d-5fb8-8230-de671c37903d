"""FastAPI application: mass-spectrometry isotope-envelope precompute API."""

from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field

from . import __version__, db
from .config import settings
from .envelope import isotope_envelope, monoisotopic_mass, top_peaks
from .isotopes import ISOTOPES
from .parser import FormulaError, parse_formula


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Dependencies are installed at image build time; the SQLite schema is
    # created here so a fresh container is ready to serve immediately.
    db.init_db(settings.db_path)
    yield
    db.close_db()


app = FastAPI(
    title="MS Precompute API",
    version=__version__,
    description=(
        "Computes natural-abundance isotope envelopes (mass/intensity "
        "polynomial convolution) and monoisotopic masses for molecular "
        "formulas, with optional pure-isotope labels such as [13C]."
    ),
    lifespan=lifespan,
)


class EnvelopeRequest(BaseModel):
    formula: str = Field(
        ...,
        description="Molecular formula, e.g. 'C6H12O6' or 'C6H6[13C]2'",
        examples=["C6H12O6"],
    )
    top_n: int = Field(
        settings.top_n_default,
        ge=1,
        le=settings.top_n_max,
        description="Number of most intense peaks to return",
    )


class PeakOut(BaseModel):
    mass: float
    intensity: float


class EnvelopeResponse(BaseModel):
    formula: str
    composition: dict[str, int]
    labels: dict[str, int]
    monoisotopic_mass: float
    base_peak_intensity: float
    peak_count: int
    peaks: list[PeakOut]
    cached: bool


@app.exception_handler(FormulaError)
async def formula_error_handler(request: Request, exc: FormulaError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "detail": {
                "error": "invalid_formula",
                "message": exc.message,
                "position": exc.position,
            }
        },
    )


@app.get("/health")
def health() -> dict:
    return {"status": "ok", "version": __version__}


@app.get("/api/v1/elements")
def elements() -> dict:
    """List supported element symbols and their natural isotopes."""
    return {
        symbol: [{"mass": m, "abundance": a} for m, a in isotopes]
        for symbol, isotopes in sorted(ISOTOPES.items())
    }


@app.post("/api/v1/envelope", response_model=EnvelopeResponse)
def envelope(req: EnvelopeRequest) -> dict:
    """Compute the isotope envelope and monoisotopic mass of a formula."""
    try:
        parsed = parse_formula(req.formula)
    except FormulaError as exc:
        db.log_request(req.formula, req.top_n, ok=False, message=exc.message)
        raise

    cached = db.get_cached(req.formula, req.top_n)
    if cached is not None:
        db.log_request(req.formula, req.top_n, ok=True, message="cache hit")
        return {**cached, "cached": True}

    env = isotope_envelope(parsed, settings.merge_tolerance, settings.prune_threshold)
    peaks = top_peaks(env, req.top_n)
    result = {
        "formula": req.formula,
        "composition": parsed.composition,
        "labels": {f"[{mass_no}{symbol}]": count for (symbol, mass_no), count in parsed.labels.items()},
        "monoisotopic_mass": round(monoisotopic_mass(parsed), 6),
        "base_peak_intensity": 100.0,
        "peak_count": len(peaks),
        "peaks": [
            {"mass": round(m, 6), "intensity": round(i, 6)} for m, i in peaks
        ],
        "cached": False,
    }
    db.put_cache(req.formula, req.top_n, {k: v for k, v in result.items() if k != "cached"})
    db.log_request(req.formula, req.top_n, ok=True, message="computed")
    return result
