import unittest

import utils

class test_validDate(unittest.TestCase):


	def test_long_date_string(self):
		long_date = "2015-01-03 12:00"
		self.assertEqual(utils.validDate(long_date), 0)

	def test_short_date_string(self):
		short_date = "2015"
		self.assertRaises(ValueError, utils.validDate, short_date)

	def test_month_in_range(self):
		bad_months = ["2015000100", "2015130100"]
		codes = [utils.validDate(x) for x in bad_months]
		self.assertItemsEqual([0, 0], codes)

	def test_hour_in_range(self):
		bad_hour = "2015010124"
		self.assertEqual(utils.validDate(bad_hour), 0)









if __name__ == '__main__':
	unittest.main()