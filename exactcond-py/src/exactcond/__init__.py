"""exactcond: exact conditional reference distributions for tables and selected linear models.

Reference distributions exact at every sample size, obtained by conditioning on a statistic
sufficient for the nuisance parameter rather than by an asymptotic approximation: the Fisher-Yates
law of a two-way table given both margins, sampled by random permutation of the column labels; a
two-table form conditioning further on the total score to cancel a common association parameter;
and the uniform law of the residual direction on the sphere, which stays exact when the statistic
is maximised over candidate structures and so remains valid after selection by fit. ctypes binding
to a shared C back end also used by the R package 'exactcond'. Pure standard library.
"""
from ._core import (ec_scores, ec_rc_null, ec_rc_test, ec_rc_two, ec_sphere_maxquad, ec_pvalue,
                    ec_multi_test, ec_diffnet_test)
__all__ = ["ec_scores", "ec_rc_null", "ec_rc_test", "ec_rc_two", "ec_sphere_maxquad", "ec_pvalue",
           "ec_multi_test", "ec_diffnet_test"]
__version__ = "0.2.0"
