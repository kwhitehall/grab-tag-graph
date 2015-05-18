""" Collection of unittests for the utils module. """

import unittest

import utils

class TestValidDate(unittest.TestCase):


    def test_poor_date_format(self):
        longDate = "2015-01-03 12:00"
        self.assertFalse(utils.valid_date(longDate))

    def test_short_date_string(self):
        shortDate = "20150101"
        self.assertFalse(utils.valid_date(shortDate))

    def test_bad_month(self):
        badMonth = "2015000100"
        self.assertFalse(utils.valid_date(badMonth))

    def test_bad_hour(self):
        badHour = "2015010124"
        self.assertFalse(utils.valid_date(badHour))

    def test_good_input_hours_only(self):
        goodDatetime = "2015010100"
        self.assertTrue(utils.valid_date(goodDatetime))

    def test_good_input_hours_minutes(self):
        goodHoursMinutes = "201501010059"
        self.assertTrue(utils.valid_date(goodHoursMinutes))


if __name__ == '__main__':
    unittest.main()
