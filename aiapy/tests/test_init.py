import aiapy


def test_citation():
    assert aiapy.__citation__ not in [None, ""]
