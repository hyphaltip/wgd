import os
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from wgd.codeml import Codeml, ncbi_table_to_icode
from wgd.core import SequenceData, GeneFamily, get_gene_families
from wgd.tree import parse_translation_parameters


def _write_cds(path, seq="ATGAAACTGGGTAAATAA"):
    with open(path, "w") as f:
        f.write(">g1\n{}\n".format(seq))


def test_ncbi_table_to_icode():
    assert ncbi_table_to_icode(1) == 0
    assert ncbi_table_to_icode(12) == 8
    assert ncbi_table_to_icode("12") == 8
    with pytest.raises(ValueError):
        ncbi_table_to_icode(99)
    with pytest.raises(ValueError):
        ncbi_table_to_icode("x")


def test_sequencedata_table(tmpdir):
    fa = os.path.join(str(tmpdir), "sp.fa")
    _write_cds(fa)  # ATG AAA CTG GGT AAA TAA
    kw = {"out_path": str(tmpdir), "tmp_path": str(tmpdir), "to_stop": True, "cds": True}
    d1 = SequenceData(fa, **kw)
    d12 = SequenceData(fa, table=12, **kw)
    assert d1.icode == 0 and d12.icode == 8
    assert str(d1.pro_sequence["sp.fa_00000"]) == "MKLGK"
    assert str(d12.pro_sequence["sp.fa_00000"]) == "MKSGK"
    with pytest.raises(ValueError):
        SequenceData(fa, table=15, **kw)


def test_codeml_icode(tmpdir):
    c = Codeml(None, tmp=str(tmpdir))
    assert c.control["icode"] == 0
    c = Codeml(None, tmp=str(tmpdir), icode=8)
    assert c.control["icode"] == 8
    assert "icode = 8" in str(c)


def test_genefamily_icode(tmpdir):
    gf = GeneFamily("GF1", {}, {}, os.path.join(str(tmpdir), "GF1"), icode=8)
    assert gf.icode == 8
    assert GeneFamily("GF2", {}, {}, os.path.join(str(tmpdir), "GF2")).icode == 0


def test_parse_translation_parameters():
    assert parse_translation_parameters("True,False") == ("True", "False", 1)
    assert parse_translation_parameters("True, False, 12") == ("True", "False", 12)


def test_get_gene_families_passes_icode(tmpdir):
    import pandas as pd
    fa = os.path.join(str(tmpdir), "sp.fa")
    with open(fa, "w") as f:
        f.write(">g1\nATGAAACTGGGTAAATAA\n>g2\nATGAAACTGGGTAAATAA\n")
    kw = {"out_path": str(tmpdir), "tmp_path": str(tmpdir)}
    fam = pd.DataFrame({"sp": [["g1", "g2"]]}, index=["GF1"])
    assert get_gene_families(SequenceData(fa, table=12, **kw), fam)[0].icode == 8
    assert get_gene_families(SequenceData(fa, **kw), fam)[0].icode == 0
