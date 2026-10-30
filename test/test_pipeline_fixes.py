import os
import pytest
from wgd.codeml import _find_trees

TREES = """tree length for dN:       0.0697
tree length for dS:       4.2333

dS tree:
(A: 0.15, (B: 0.67, C: 0.78): 0.14);
dN tree:
(A: 0.00, (B: 0.01, C: 0.01): 0.00);

w ratios as node labels:
(A #0.0001 , (B #0.017 , C #0.019 ) #0.03 );
"""


@pytest.mark.parametrize("trailer", ["", "\n", "\nlnL(ntime: 5)\n"])
def test_find_trees_independent_of_trailing_lines(trailer):
    ds, dn, w = _find_trees((TREES + trailer).splitlines(keepends=True))
    assert ds.startswith("(A: 0.15")
    assert dn.startswith("(A: 0.00")
    assert w.startswith("(A: 0.0001")


def test_find_trees_missing():
    with pytest.raises(ValueError):
        _find_trees(["nothing here\n"])


def test_make_diamond_db_concurrent(tmpdir):
    """Many processes building one database must all see a complete file."""
    import shutil
    from multiprocessing import Pool
    if shutil.which("diamond") is None:
        pytest.skip("diamond not available")
    from wgd.core import SequenceData
    fa = os.path.join(str(tmpdir), "sp.fa")
    with open(fa, "w") as f:
        for i in range(50):
            f.write(">g{}\nATG{}TAA\n".format(i, "GCT" * (10 + i)))
    s = SequenceData(fa, out_path=str(tmpdir), tmp_path=str(tmpdir), to_stop=True, cds=True)
    with Pool(8) as pool:
        pool.map(_build, [s] * 8)
    assert os.path.isfile(s.pro_db + ".dmnd")
    assert not [x for x in os.listdir(os.path.dirname(s.pro_db)) if x.endswith(".dmnd") and x != os.path.basename(s.pro_db) + ".dmnd"]


def _build(s):
    s.make_diamond_db()
