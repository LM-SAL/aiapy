import aiapy


def test_citation() -> None:
    assert aiapy.__citation__ not in [None, ""]
