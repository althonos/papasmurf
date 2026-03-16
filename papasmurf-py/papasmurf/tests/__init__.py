from . import test_doctest, test_lib

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_lib))
    test_doctest.load_tests(loader, suite, pattern)
    return suite
