import os

# If this package has tests data in the tests/data directory, add them to
# the paths here, see commented example
paths = ['coveragerc',
         os.path.join('data', '*txt'),
         ]


def get_package_data():
    return {'aiapy.tests': paths}
