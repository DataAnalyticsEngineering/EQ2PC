"""
Unit tests for eq2pc repo
"""


import src.eq2pc as e2


class TestEQ2PC:

    def test_root_gen(self):
        sets1 = e2.Ceq_search([12], [5])
        sets2 = e2.Ceq_search([4, 3], [5])
        #sets3 = e2.Ceq_search([2, 3, 4], [5])

        out1 = e2.plot_S(sets1[0][0], show=False)
        out2 = e2.plot_S(sets2[0][0], show=False)
        #out3 = e2.plot_S(sets3[0][0], show=False)

        S1 = sets2[0][0]
        C1s = e2.rCs(S1)
        out4 = e2.plot_C(C1s[0], show=False)

        assert len(sets1) > 0 and len(sets2) > 0 and out1 == 0 and out2 == 0 and out4 == 0

    def test_child_gen(self):
        sets = e2.Ceq_search([4, 3], [5])
        Ss = sets[0]
        Ss = [e2.kbe(S) for S in Ss]
        Ss = [e2.phc(S) for S in Ss]
        assert e2.Ceq(Ss[0], Ss[1])

    def test_save(self):
        folder = e2.Ceq_search_save([12], [5])
        assert len(folder) > 0
