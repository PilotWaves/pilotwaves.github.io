import numpy

# It's important to remember that dat-files, being loaded by numpy.loadtxt()
# do not always have a consistent shape or type with MATLAB originals.
# For instance, a column vector (a 2d-matrix of shape (x, 1)) or
# a row vector (a 2d-matrix of shape (1, x)) in MATLAB
# will become just a numpy-array (a 1d-vector of shape (x, )) in NumPy.
# Numbers will become a "zero-dimentional" numpy arrays --
# matrices of shape () (they actually have "empty" shape, it's not a typo)).
# 2d-matrices thicker than one in both dimensions are likely to keep the right shape and transposition.
# So it's important to bring all the variables back to a NumPy form
# consistent with MATLAB, like I did below, before creating an autograder archive.

# So, according to MATLAB's namespace, in a correct solution...
# A1 should be a column vector
A1 = numpy.loadtxt("A1.dat").reshape((-1, 1))
# A2 should be a row vector
A2 = numpy.loadtxt("A2.dat").reshape((1, -1))
# A3 should be a number
A3 = float(numpy.loadtxt("A3.dat"))
# A4 should be a column vector
A4 = numpy.loadtxt("A4.dat").reshape((-1, 1))
# A5 should be a row vector
A5 = numpy.loadtxt("A5.dat").reshape((1, -1))
# A6 should be a number
A6 = float(numpy.loadtxt("A6.dat"))
# A7 is a matrix
A7 = numpy.loadtxt("A7.dat")
# A8 should be a number
A8 = float(numpy.loadtxt("A8.dat"))
# A9 should be a number
A9 = float(numpy.loadtxt("A9.dat"))
# A10 should be a number
A10 = float(numpy.loadtxt("A10.dat"))
# A11-A15 are matrices
A11 = numpy.loadtxt("A11.dat")
A12 = numpy.loadtxt("A12.dat")
A13 = numpy.loadtxt("A13.dat")
A14 = numpy.loadtxt("A14.dat")
A15 = numpy.loadtxt("A15.dat")

test_suite = [
    {
        "test_name": "A1",
        "variable_name": "A1",
        "description": "Problem 1a -- Numerical solution",
        "score": 1
    },
    {
        "test_name": "A2",
        "variable_name": "A2",
        "description": "Problem 1a -- Error values.",
        "score": 1
    },
    {
        "test_name": "A3",
        "variable_name": "A3",
        "description": "Problem 1a -- Slope.",
        "score": 1
    },
    {
        "test_name": "A4",
        "variable_name": "A4",
        "description": "Problem 1b -- Numerical solution",
        "score": 1
    },
    {
        "test_name": "A5",
        "variable_name": "A5",
        "description": "Problem 1b -- Error values.",
        "score": 1
    },
    {
        "test_name": "A6",
        "variable_name": "A6",
        "description": "Problem 1b -- Slope.",
        "rtol": 1e-5,
        "atol": 1e-2,
        "score": 1
    },
    {
        "test_name": "A7",
        "variable_name": "A7",
        "hint_wrong_size": "Problem 2a -- VdP Solutions",
        "score": 1
    },
    {
        "test_name": "A8",
        "variable_name": "A8",
        "hint_wrong_size": "Problem 2b -- VdP slope ode45",
        "score": 1
    },
    {
        "test_name": "A9",
        "variable_name": "A9",
        "hint_wrong_size": "Problem 2b -- VdP slope ode23",
        "score": 1
    },
    {
        "test_name": "A10",
        "variable_name": "A10",
        "hint_wrong_size": "Problem 2b -- VdP slope ode113",
        "score": 1
    },
    {
        "test_name": "A11",
        "variable_name": "A11",
        "hint_wrong_size": "Problem 3 -- Neurons (1)",
        "score": 1
    },
    {
        "test_name": "A12",
        "variable_name": "A12",
        "hint_wrong_size": "Problem 3 -- Neurons (2)",
        "score": 1
    },
    {
        "test_name": "A13",
        "variable_name": "A13",
        "hint_wrong_size": "Problem 3 -- Neurons (3)",
        "score": 1
    },
    {
        "test_name": "A14",
        "variable_name": "A14",
        "hint_wrong_size": "Problem 3 -- Neurons (4)",
        "score": 1
    },
    {
        "test_name": "A15",
        "variable_name": "A15",
        "hint_wrong_size": "Problem 3 -- Neurons (5)",
        "score": 1
    },
]

number_of_attempts = 10

matlab_credentials = '~/Storage/repos/gspack_uw_amath_matlab_credentials'
matlab_use_template = True
