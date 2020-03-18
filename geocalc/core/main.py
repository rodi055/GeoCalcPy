import unittest
import geocalc

class CoordinatesConverterTests(unittest.TestCase):
    def test_lla2ecef(self):
        lat, lon, alt = 30, 31, 100
        ecef = geocalc.lla2ecef(lat, lon, alt)
        x, y , z = ecef[0], ecef[1], ecef[2]
        self.assertAlmostEqual(x, 4738715.05, places=1)
        self.assertAlmostEqual(y, 2847307.26, places=1)
        self.assertAlmostEqual(z, 3170423.73, places=1)
    
    def test_ecef2lla(self):
        x, y, z = 4738715.05, 2847307.26, 3170423.73
        lat, lon, alt = geocalc.ecef2lla([x, y, z])
        self.assertAlmostEqual(lat, 30, places=1)
        self.assertAlmostEqual(lon, 31, places=1)
        #self.assertAlmostEqual(alt, 100, places=1)

    def test_ecef2ned(self):
        lat_ref, lon_ref, alt_ref = 29, 30, 0
        x, y, z = 4738715.05, 2847307.26, 3170423.73
        ned = geocalc.ecef2ned([x, y, z], lat_ref, lon_ref, alt_ref)
        north, east, down = ned[0], ned[1], ned[2]
        self.assertAlmostEqual(north, 93129.29, places=1)
        self.assertAlmostEqual(east, 96482.89, places=1)
        self.assertAlmostEqual(down, -6371513.44, places=1)

    def test_ned2ecef(self):
        lat_ref, lon_ref, alt_ref = 29, 30, 0
        north, east, down = 93129.29, 96482.89, -6371513.44
        ecef = geocalc.ned2ecef([north, east, down], lat_ref, lon_ref, alt_ref)
        x, y, z = ecef[0], ecef[1], ecef[2]
        self.assertAlmostEqual(x, 4738715.05, places=1)
        self.assertAlmostEqual(y, 2847307.26, places=1)
        self.assertAlmostEqual(z, 3170423.73, places=1)

    def test_ned2lla(self):
        lat_ref, lon_ref, alt_ref = 29, 30, 0
        north, east, down = 93129.29, 96482.89, -6371513.44
        lat, lon, alt = geocalc.ned2lla([north, east, down], lat_ref, lon_ref, alt_ref)
        self.assertAlmostEqual(lat, 30, places=1)
        self.assertAlmostEqual(lon, 31, places=1)
        #self.assertAlmostEqual(alt, 100, places=1)

    '''def test_lla2ned(self):
        lat_ref, lon_ref, alt_ref = 29, 30, 0
        x, y, z = 4738715.05, 2847307.26, 3170423.73
        ned = geocalc.ned2ecef([x, y, z], lat_ref, lon_ref, alt_ref)
        north, east , down = ned[0], ned[1], ned[2]
        self.assertAlmostEqual(north, 4738715.05, places=1)
        self.assertAlmostEqual(east, 2847307.26, places=1)
        self.assertAlmostEqual(down, 3170423.73, places=1)'''


if __name__ == '__main__':
    unittest.main()