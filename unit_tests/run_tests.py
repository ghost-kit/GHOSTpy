import unittest
import invariants_units


suite = unittest.TestLoader().loadTestsFromModule(invariants_units)
unittest.TextTestRunner().run(suite)


