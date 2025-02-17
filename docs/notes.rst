Notes
=====

Reproducibility
---------------

The randomness in DNA Chisel is entirely determined by the Numpy random generator.
As a consequence, setting the Numpy seed at the very beginning of a script ensures
that the result will always be the same:

.. code:: python

    import numpy
    numpy.random.seed(123)

    # Start the script there...

If you have reproducibility issues despite setting the seed, we would consider
it a bug, please open an issue on Github. Please note that the seed must be set
*before* calling any part of DNA Chisel (e.g. `DnaOptimizationProblem`).
