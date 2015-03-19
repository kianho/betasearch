import unittest
import betasearch

class test_betasearch(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_Elbow(self):
        from betasearch import Elbow

        elbow_A = Elbow(1, 2)
        elbow_B = Elbow(1, 2)
        elbow_C = Elbow(2, 1)

        self.assertTrue(elbow_A.same(elbow_B))
        self.assertFalse(elbow_A.same(elbow_C))

        self.assertRaises(ValueError, lambda : Elbow(-1, -1))
        self.assertRaises(ValueError, lambda : Elbow(-1, 0))
        self.assertRaises(ValueError, lambda : Elbow(-1, 2))

    def test_Span(self):
        from betasearch import Span

        span_A = Span(1, 2)
        span_B = Span(1, 2)
        span_C = Span(2, 1)

        self.assertTrue(span_A.same_direction(span_B))
        self.assertFalse(span_A.same_direction(span_C))

        self.assertTrue(span_A.overlaps(span_B))
        self.assertTrue(span_A.overlaps(span_C))

        self.assertRaises(ValueError, lambda : Span(0, 0))
        self.assertRaises(ValueError, lambda : Span(-1, 0))
        self.assertRaises(ValueError, lambda : Span(-1, 2))

if __name__ == "__main__":
    unittest.main()
