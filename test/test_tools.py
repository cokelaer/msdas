from msdas import *



def test_get_protein():

    t = SequenceTools()
    name = 'sp|Q8IYB3-2|SRRM1_HUMAN'
    res = t.get_protein_from_fasta_header(name)
    assert res == 'SRRM1-2'

def test_accession():
    t = SequenceTools()
    name = 'sp|Q8IYB3-2|SRRM1_HUMAN'
    res = t.get_proteinId_from_fasta_header(name)
    assert res == 'Q8IYB3-2'


def test_class_requires():

    class A(Requires):
        def __init__(self):
            super(A, self).__init__()
            self.a = []
        def testa(self):
            self._requires("a", "required")
        def testb(self):
            self._requires("b", "required")
        def testalist(self):
            self._requires(["a"], "required")
        def testinteger(self):
            self._requires(1, "required")


    a = A()
    a.testa()
    a.testalist()
    try:
        a.testb()
        assert False
    except:
        assert True

    try:
        a.testinteger()
        assert False
    except:
        assert True



def test_modifications():



    m = Modification()
    m.modification
    assert m.get_individual_psites("[3] S+79.9663") == ["S3"]
    assert m.modif2psites("[5] S+7777| [1] T+5555") == "T1^S5"
    try:
        # modif not provided neither during the object instanciation nor durig
        # the call to the method
        m.get_individual_psites()
        assert False
    except:
        assert True

    m = Modification(m.example)
    assert m.psites == "S3^S5"
    assert m.get_individual_psites() == ["S3", "S5"]
