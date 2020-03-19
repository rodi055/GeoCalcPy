import sys
import unittest


#
# from core import geocalc
# sys.path.insert(0, "./")

class CoordinatesConverterTests(unittest.TestCase):

    def test_lla2ecef(self):
        lat, lon, alt = 30, 31, 100
        from core import geocalc
        ecef = geocalc.lla2ecef(lat, lon, alt)
        x, y, z = ecef[0], ecef[1], ecef[2]
        self.assertAlmostEqual(x, 4738715.05, places=1)
        self.assertAlmostEqual(y, 2847307.26, places=1)
        self.assertAlmostEqual(z, 3170423.73, places=1)


if __name__ == '__main__':
    unittest.main()
