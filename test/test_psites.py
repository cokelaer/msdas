from msdas import *
from msdas.psites import PSites
from easydev import gsf
import os


def test_psites():

    p = PSites()

    #assert p.get_psites_yeast("S3 S4") == "S3^S4"

    assert p.isvalid("S3^S4") == True
    assert p.isvalid("S3+S4") == True
    assert p.isvalid("S3^S4+T1") == True
    assert p.isvalid("") == False
    assert p.isvalid("S") == False
    assert p.isvalid("1") == False


    assert p.isvalid("S3S4") == False
    assert p.isvalid("S3 S4") == False
    assert p.isvalid("S3;S4") == False


    assert p.sorted(["S3", "S1", "S4+S5"]) == 'S1^S3^S4+S5'
    try:
        p.sorted(["S", "S1", "S4+S5"])
        assert False
    except:
        assert True

    sequence = "KVS(Phospho)VGSM(Oxidation)GS(Phospho)GK"
    psites = "S956^S962"
    assert p.append_oxidation(sequence, psites) == 'S956^M960^S962'
    try:
        p.append_oxidation("KVSM", "2")
        assert False
    except:
        assert True


    assert p.remove_duplicated("S3^S3") == "S3"
    assert p.get_unique_psites(["S145", "S145^S145", "S145", "S145+S145"]) == ["S145"]
    assert p.get_factorised_psites(['S177^Y182', 'T180^Y182']) == 'S177+T180^Y182'

    assert p.remove_spaces("S3 ; S4") == "S3;S4"
    
    assert p.sort_psites_ors_only("S4+S3") == "S3+S4"
    assert p.sort_psites_ors_only("S4^S3") == "S4^S3"  # not sorted because there is an AND
    
    p.append_oxidation("KVS(Phospho)VGSMGS(Phospho)GK", "S956^S962")

