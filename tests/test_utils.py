""" Collection of unittests for the utils module. """

import unittest

import utils

class TestValidDate(unittest.TestCase):


    def test_long_date_string(self):
        longDate = "2015-01-03 12:00"
        self.assertFalse(utils.valid_date(longDate))

    def test_short_date_string(self):
        shortDate = "2015"
        self.assertFalse(utils.valid_date(shortDate))

    def test_month_in_range(self):
        badMonth = "2015000100"
        self.assertFalse(utils.valid_date(badMonth))

    def test_hour_in_range(self):
        badHour = "2015010124"
        self.assertFalse(utils.valid_date(badHour))

    def test_good_input(self):
        goodDatetime = "2015010100"
        self.assertTrue(utils.valid_date(goodDatetime))


if __name__ == '__main__':
    unittest.main()
