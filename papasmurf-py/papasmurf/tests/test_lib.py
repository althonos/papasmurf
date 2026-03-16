import gzip
import os
import tempfile
import unittest

import papasmurf



class TestBuilder(unittest.TestCase):

    def test_builder_error(self):
        builder = papasmurf.Builder([ ("TGGCGAA", "CCGTG") ])
        with self.assertRaises(TypeError):
            builder.add(1, 2)
        with self.assertRaises(ValueError):
            builder.add("Bacteroides uniformis", "nonsense non-DNA string")

