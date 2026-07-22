from math import sqrt

import pytest

from s13_mixed_model import fixed_effect_meta


COHORT_EFFECTS = {
    "GSE308231": {"diff": -0.7274126800109078, "se": 0.6880011057762103},
    "GSE246662": {"diff": -1.9104335451243588, "se": 0.37116697424314604},
    "GSE239676": {"diff": -0.9093539555709704, "se": 0.28889380957749605},
}


@pytest.mark.parametrize(
    ("omitted", "expected"),
    [
        (None, -1.2316571205236306),
        ("GSE308231", -1.2870233105385795),
        ("GSE246662", -0.8820827423026468),
        ("GSE239676", -1.643740708123175),
    ],
)
def test_locked_validation_fixed_effect_and_loco_values(omitted, expected):
    cohorts = [c for c in COHORT_EFFECTS if c != omitted]
    beta, se, p = fixed_effect_meta(COHORT_EFFECTS, cohorts)
    assert beta == pytest.approx(expected, abs=1e-12)
    assert se > 0
    assert 0 <= p <= 1


def test_fixed_effect_standard_error_formula():
    expected = sqrt(1 / sum(1 / row["se"] ** 2 for row in COHORT_EFFECTS.values()))
    _, se, _ = fixed_effect_meta(COHORT_EFFECTS, list(COHORT_EFFECTS))
    assert se == pytest.approx(expected, abs=1e-15)
