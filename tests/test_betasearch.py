import unittest
import betasearch

class test_betasearch(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_Elbow(self):
        elbow_A = betasearch.Elbow(1, 2)
        elbow_B = betasearch.Elbow(1, 2)
        elbow_C = betasearch.Elbow(2, 1)

        self.assertTrue(elbow_A.same(elbow_B))
        self.assertFalse(elbow_A.same(elbow_C))

if __name__ == "__main__":
    unittest.main()
