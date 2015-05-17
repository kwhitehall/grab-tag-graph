import unittest

import utils

class test_validDate(unittest.TestCase):


	def test_long_date_string(self):
		long_date = "2015-01-03 12:00"
		self.assertRaises(ValueError, utils.validDate, long_date)

	def test_short_date_string(self):
		short_date = "2015"
		self.assertRaises(ValueError, utils.validDate, short_date)

	def test_month_in_range(self):
		bad_month = "2015000100"
		self.assertRaises(ValueError, utils.validDate, bad_month)

	def test_hour_in_range(self):
		bad_hour = "2015010124"
		self.assertRaises(ValueError, utils.validDate, bad_hour)

	def test_good_input(self):
		good_datetime = "2015010100"
		self.assertEqual(utils.validDate(good_datetime), 1)









if __name__ == '__main__':
	unittest.main()