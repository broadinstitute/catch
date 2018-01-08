"""Functions for optimizing parameter values across datasets.

We would like to find optimal parameter values, with different values
for each dataset. We can do this by treating the problem as a constrained
nonlinear optimization problem. In particular, we define a loss over
the parameter values and seek to minimize the loss, subject to the
constraint that the total number of probes (using those parameter values)
is less than the maximum number of allowed probes. It enforces this
constraint using a barrier function, and uses scipy's optimize module
to minimize the sum of the loss and barrier function.
"""

import logging

from hybseldesign.pool import interpolate_count as ic

import numpy as np
from scipy import optimize

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def _make_loss_fn(probe_counts, max_total_count, coeffs, weights,
        interp_fn_type='standard'):
    """Generate and return a loss function.

    The function calculates a loss over the parameters and adds onto
    that loss to meet a constraint based on the barrier method. It
    uses a logarithmic barrier function to enforce the constraint that
    the total probe count be <= max_total_count.

    The loss over the parameters is:
        sum_{datasets d} (w_d * (sum_{param j} (c_j * (v_{dj})^2)))
    where v_{dj} is the value of the j'th parameter for dataset d,
    c_j is the coefficient for the j'th parameter, and w_d is the
    weight of dataset d

    Args:
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters
        max_total_count: upper bound on the number of total probes
        coeffs: coefficient in the loss function for each parameter, in
            the same order as the parameters are given
        weights: dict giving weight in the loss function for each dataset
        interp_fn_type: 'standard' (only perform interpolation on mismatches
            and cover_extension parameters) or 'nd' (use scipy's interpolate
            package to interpolate over n-dimensions)

    Returns:
        a function that is the sum of a loss defined over the parameters
        and a value designed to enforce a barrier on the total number
        of probes
    """
    total_probe_count_across_datasets = ic._make_total_probe_count_across_datasets_fn(
        probe_counts, interp_fn_type=interp_fn_type)

    def loss(x, *func_args):
        """
        Compute a loss.

        Let the number of datasets (len(probe_counts)) be N.
        x is a list giving all the parameter values across datasets,
        such that x_i is the (i % N)'th parameter of the (i/N)'th dataset,
        for i=0,1,2,...
        """
        num_datasets = len(probe_counts)
        # The number of parameter values must be a multiple of the number
        # of datasets
        assert len(x) % num_datasets == 0

        num_params = int(len(x) / num_datasets)

        # There must be a coefficient for each parameter
        assert len(coeffs) == num_params

        # First compute a loss over the parameters by taking their L2-norm
        # This is the function we really want to minimize
        opt_val = 0
        for i, dataset in enumerate(sorted(probe_counts.keys())):
            opt_val_dataset = 0
            for j in range(num_params):
                v = x[num_params * i + j]
                opt_val_dataset += coeffs[j] * np.power(v, 2.0)
            opt_val += weights[dataset] * opt_val_dataset

        # We also have the constraint that the total probe count be less than
        # max_total_count
        # We add a barrier function to enforce this constraint and weight the
        # barrier by eps
        eps = func_args[0]
        total_probe_count = total_probe_count_across_datasets(x)
        if np.isnan(total_probe_count):
            # If the interp_fn_type is 'nd' and the parameter values are
            # outside the convex hull of computed points (from probe_counts)
            # for a dataset, scipy's interpolator will be unable to
            # interpolate a probe count and will return nan; here, make
            # the loss (through the barrier) high to reflect this
            logger.warning(("Parameter values being searched are outside "
                "the convex hull of computed points; unable to interpolate "
                "a probe count"))
            barrier_val = 10000000
        elif total_probe_count >= max_total_count:
            # Since the count is beyond the barrier, we should in theory
            # return infinity. But if the optimizer does indeed try parameters
            # that put the probe count here, it would be unable to compute
            # an approximate gradient and may get stuck. So help it out
            # by giving a value such that the negative gradient points toward
            # a direction outside the barrier.
            # Add 1 so that, if total_probe_count == max_total_count, we do
            # not take log(0).
            barrier_val = 9999 + 10000.0 * np.log((total_probe_count -
                                                  max_total_count + 1))
        else:
            # The barrier function is -log(max_total_count - total_probe_count), to
            # enforce the constraint that total_probe_count be less than
            # max_total_count.
            # Add 1 so that, if max_total_count - total_probe_count < 1,
            # the argument to log(..) remains >= 1.
            barrier_val = -1.0 * eps * np.log((max_total_count -
                                               total_probe_count + 1))

        return opt_val + barrier_val

    return loss


def _make_param_bounds_standard(probe_counts, step_size=0.001):
    """Calculate bounds on parameter values for only mismatches and cover_extension.

    For each dataset d, this calculates bounds on the values of the
    mismatches and cover extension parameters based on which values
    have a number of probes calculated for them. Namely, we wish to
    ensure we can find a bounding box around an arbitrary point
    (mismatches, cover_extension); based on the values in probe_counts,
    the bounds should ensure this.

    Args:
        probe_counts: dict giving the number of probes for each dataset
            and choice of parameters
        step_size: small value subtracted from each upper bound so
            that the minimizer does not exceed it

    Returns:
        [m_1, e_1, m_2, e_2, ...] where m_i is a tuple (lo, hi) that
        gives bounds on the number of mismatches for the i'th dataset,
        and likewise for e_i
    """
    bounds = []
    for dataset in sorted(probe_counts.keys()):
        params = probe_counts[dataset].keys()

        # This requires that mismatches and cover_extension be the
        # only two parameters
        for p in params:
            assert len(p) == 2

        # Bound cover_extensions by the lowest and highest value for
        # which we have a probe count result
        cover_extensions = [k[1] for k in params]
        cover_extensions_lo = min(cover_extensions)
        cover_extensions_hi = max(cover_extensions)

        # To ensure we can find a rectangular bounding box around an
        # arbitrary point ((mismatches, cover_extension)), our lower
        # bound on mismatches should have a corresponding cover extension
        # of min(cover_extensions) and of max(cover_extensions);
        # so should our upper bound on mismatches
        mismatches = [k[0] for k in params]
        mismatches_with_valid_cover_extension = \
            [m for m in mismatches if ((m, cover_extensions_lo) in params and
             (m, cover_extensions_hi) in params)]
        mismatches_lo = min(mismatches_with_valid_cover_extension)
        mismatches_hi = max(mismatches_with_valid_cover_extension)

        bounds += [(mismatches_lo, mismatches_hi - step_size)]
        bounds += [(min(cover_extensions), max(cover_extensions) - step_size)]
    return bounds


def _make_param_bounds_nd(probe_counts, step_size=0.001):
    """Calculate bounds on parameter values in n dimensions.

    For each dataset d, this calculates bounds on the values of each
    parameter based only on the min/max of what has been computed
    (i.e., is in probe_counts). Note that a point that satisfies these
    bounds may not necessarily be within the convex hull of all the
    computed points.

    Args:
        probe_counts: dict giving the number of probes for each dataset
            and choice of parameters
        step_size: small value subtracted from each upper bound so that
            the minimizer does not exceed it

    Returns:
        list x giving all the parameter values across datasets,
        such that x_i corresponds to the (i % N)'th parameter of the
        (i/N)'th dataset, for i=0,1,2,... where N is the number of datasets;
        x_i is a tuple (lo, hi) giving bounds on the parameter
    """
    bounds = []
    for dataset in sorted(probe_counts.keys()):
        params = list(probe_counts[dataset].keys())
        num_params = len(params[0])

        for j in range(num_params):
            lo = min(params[i][j] for i in range(len(params)))
            hi = max(params[i][j] for i in range(len(params))) - step_size
            bounds += [(lo, hi)]
    return bounds
            

def _make_initial_guess(probe_counts, bounds, num_params):
    """Make initial guess for optimal parameter values.

    This guesses a value for each parameter separately
    for each dataset, choosing uniformly from within the provided
    bounds. (If bounds is None, it chooses an already computed
    point (a parameter in probe_counts) uniformly at random.) It allows
    a guess that exceeds the constraint on the total number of probes,
    but outputs a warning if this occurs.

    Args:
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters
        bounds: lower/upper bounds on each parameter for each dataset
            (bounds[i] is a tuple (lo, hi) giving the bounds for
            the (i % N)'th parameter of the (i/N)'th dataset for
            i=0,1,2... where N is the number of datasets; if None,
            rather than picking a guess uniformly at random from
            within the bounds, this picks an already computed
            point (a paramter value in probe_counts) uniformly
            at random, which ensures the guess is within the convex
            hull of the computed points
        num_params: number of parameters

    Returns:
        list x giving all the initial guesses across datasets,
        such that x_i corresponds to the (i % N)'th parameter of the
        (i/N)'th dataset, for i=0,1,2,... where N is the number of datasets
    """
    num_datasets = len(probe_counts)
    if bounds is not None:
        # The number of bounds should be a multiple of the number of datasets
        assert len(bounds) % num_datasets == 0
        assert num_params == int(len(bounds) / num_datasets)

    x0 = np.zeros(num_datasets * num_params)
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        if bounds is not None:
            # For each parameter, pick a value uniformly at random
            # from within the bounds
            for j in range(num_params):
                lo, hi = bounds[num_params * i + j]
                x0[num_params * i + j] = np.random.uniform(lo, hi)
        else:
            # Pick one of the already computed points
            param_vals = list(probe_counts[dataset])
            guess = param_vals[np.random.randint(len(param_vals))]
            for j in range(num_params):
                x0[num_params * i + j] = guess[j]

    return x0


def _optimize_loss(probe_counts, loss_fn, bounds, x0,
                   initial_eps=10.0, step_size=0.001,
                   interp_fn_type='standard'):
    """Optimize loss function with barrier.

    This uses scipy's optimize.fmin_tnc to minimize the loss function
    in which the barrier is weighted by eps. It repeatedly minimizes
    the loss while decreasing eps so that, by the last iteration, the
    weight on the barrier is very small. On each iteration, it starts
    the initial guess/position at the solution to the previous iteration.

    Args:
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters
        loss_fn: the loss function provided by _make_loss_fn
        bounds: bounds on the parameter values provided by _make_param_bounds_*
        x0: the initial guess of parameter values (i.e., starting position)
        initial_eps: weight of the barrier on the first iteration
        step_size: epsilon value provided to optimize.fmin_tnc
        interp_fn_type: 'standard' (only perform interpolation on mismatches
            and cover_extension parameters) or 'nd' (use scipy's interpolate
            package to interpolate over n-dimensions)

    Returns:
        list of length (number of datasets)*(number of parameters) where
        x_i is the (i % N)'th parameter of the (i/N)'th dataset,
        for i=0,1,2,... where N=(number of datasets)
    """
    eps = initial_eps
    while eps >= 0.01:
        x0_probe_count = ic._make_total_probe_count_across_datasets_fn(
            probe_counts, interp_fn_type=interp_fn_type)(x0)
        logger.info(("Starting an iteration with eps=%f, with x0 yielding %f "
                "probes"), eps, x0_probe_count)

        sol, nfeval, rc = optimize.fmin_tnc(loss_fn, x0, bounds=bounds,
                                            args=(eps,),
                                            approx_grad=True,
                                            epsilon=step_size, disp=1, maxfun=2500)

        if rc in [0, 1, 2]:
            # rc == 0 indicates reaching the local minimum, and rc == 1 or
            # rc == 2 indicates the function value converged
            logger.info("  Iteration was successful")
        else:
            logger.info("  Iteration failed to converge!")

        x0 = sol
        eps = 0.1 * eps

    return sol


def _total_probe_count_without_interp(params, probe_counts):
    """Calculate a total probe count without interpolation.

    This assumes that params are keys in the datasets of probe_counts.

    The result of ic._make_total_probe_count_across_datasets_fn should give
    the same count as this function (if params are keys in the datasets
    of probe_counts). But this uses probe_counts directly and can be
    used as a sanity check -- i.e., it does not do any interpolation.

    Args:
        params: parameter values to use when determining probe counts;
            params[i] is the (i % N)'th parameter of the (i/N)'th dataset,
            where N is the number of datasets
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters

    Returns:
        total number of probes across all datasets, according to the
        given values of params
    """
    num_datasets = len(probe_counts)
    # The total number of parameters must be a multiple of the number
    # of datasets
    assert len(params) % num_datasets == 0

    num_params = int(len(params) / num_datasets)

    s = 0
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        p = tuple(params[num_params * i + j] for j in range(num_params))
        s += probe_counts[dataset][p]
    return s


def _round_params(params, probe_counts, max_total_count, loss_coeffs, weights,
        mismatches_eps=0.01, cover_extension_eps=0.1,
        mismatches_round=1, cover_extension_round=1):
    """Round parameter values while satisfying the constraint on total count.

    This is only applied to the mismatches and cover_extension parameters.

    Parameter values found by the search are floats. We want the mismatches
    and cover_extension parameters to be integers, or to fit on a specified
    grid.

    The floats, as given in params, should satisfy the constraint (i.e.,
    the interpolated total number of probes is less than max_total_count).
    Thus, we can round them up, because (generally) increasing the parameter
    values will decrease the number of probes; therefore, after rounding up
    they should still satisfy the constraint.

    But we also check if the parameter values are within eps of their
    rounded-down value. The loss optimizer has a tendency to make this happen
    for some parameters (e.g., finding an optimal mismatches parameter value
    of 1.00001). The reason likely has to do with the fact that, because
    we are linearly interpolating total probe counts, the gradient of the
    barrier function changes greatly around certain values (i.e., around
    the actual data values). That is, the barrier function is not at all
    smooth around actual data values. This may cause the optimizer to
    yield parameter values that are very close to parameter values for
    which probe counts have actually been computed.

    After rounding up, some parameters are decreased; we repeatedly
    choose to decrease the parameter whose reduction yields the smallest
    loss while still yielding a number of probes that is less than
    max_total_count.

    Args:
        params: parameter values to use when determining probe counts;
            params[2*i] is the number of mismatches of the i'th dataset
            and params[2*i+1] is the cover extension of the i'th dataset
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters
        max_total_count: upper bound on the number of total probes
        loss_coeffs: tuple (m, e) giving the coefficients for the
            mismatches (m) and cover_extension (e) parameters in the
            loss function
        weights: dict giving weight in the loss function for each dataset
        mismatches_eps/cover_extension_eps: eps as defined above for
            mismatches and cover_extension
        mismatches_round/cover_extension_round: round mismatches and
            cover_extension to the nearest multiple of this

    Returns:
        list in which index i corresponds to the parameter given in
        params[i], but rounded
    """
    # This requires that the only two parameters be mismatches and
    # cover_extension
    num_datasets = len(probe_counts)
    assert len(params) == 2*num_datasets
    assert len(loss_coeffs) == 2

    params_rounded = []
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]

        if mismatches - ic._round_down(mismatches, mismatches_round) < mismatches_eps:
            # Round mismatches down
            mismatches = ic._round_down(mismatches, mismatches_round)
        else:
            # Round mismatches up
            mismatches = ic._round_up(mismatches, mismatches_round)

        if cover_extension - ic._round_down(cover_extension, cover_extension_round) < cover_extension_eps:
            # Round cover_extension down
            cover_extension = ic._round_down(cover_extension, cover_extension_round)
        else:
            # Round cover_extension up
            cover_extension = ic._round_up(cover_extension, cover_extension_round)

        params_rounded += [mismatches, cover_extension]

    total_probe_count = ic._make_total_probe_count_across_datasets_fn(
        probe_counts, interp_fn_type='standard')
    # Verify that the probe count satisfies the constraint
    # Note that this assertion may fail if we are dealing with datasets
    # for which few actual probe counts have been computed; in these
    # cases, the interpolation may severely underestimate the number
    # of probes at a particular parameter choice
    assert total_probe_count(params_rounded) < max_total_count

    # Keep decreasing parameters while satisfying the constraint.
    # In particular, choose to decrease the parameter whose reduction
    # yields the smallest loss while still satisfying the constraint.
    loss_fn = _make_loss_fn(probe_counts, max_total_count, loss_coeffs,
        weights, interp_fn_type='standard')
    while True:
        curr_loss = loss_fn(params_rounded, 0)
        # Find a parameter to decrease
        min_loss, min_loss_new_params = curr_loss, None
        for i in range(len(params_rounded)):
            params_tmp = list(params_rounded)
            if params_tmp[i] == 0:
                # This cannot be decreased
                continue
            if i % 2 == 0:
                # This is a mismatch; decrease by the rounding multiple
                params_tmp[i] -= mismatches_round
            else:
                # This is a cover_extension; decrease by the rounding multiple
                params_tmp[i] -= cover_extension_round
            if total_probe_count(params_tmp) >= max_total_count:
                # This change yields too many probes, so skip it
                continue
            new_loss = loss_fn(params_tmp, 0)
            if new_loss < min_loss:
                min_loss = new_loss
                min_loss_new_params = params_tmp

        if min_loss_new_params != None:
            # There was a change that led to a better loss, so
            # update params_rounded
            params_rounded = min_loss_new_params
        else:
            # No parameter change satisfies the constraint and
            # yields an improvement in the loss
            break

    return params_rounded


def _log_params_by_dataset(params, probe_counts, type="float"):
    """Log optimal parameter values for each dataset.

    This only logs mismatches and cover_extension parameters.

    Args:
        params: parameter values to log; params[2*i] is the number of
            mismatches of the i'th dataset and params[2*i+1] is the
            cover extension of the i'th dataset
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters
        type: 'int' or 'float' specifying whether to output parameter
            values as an integer or float
    """
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]
        if type == "float":
            logger.info("%s: (%f, %f)", dataset, mismatches, cover_extension)
        elif type == "int":
            logger.info("%s: (%d, %d)", dataset, mismatches, cover_extension)
        else:
            raise ValueError("Unknown type %s", type)


def standard_search(probe_counts, max_total_count,
        verify_without_interp=False, round_params=None,
        loss_coeffs=None, dataset_weights=None):
    """Search over mismatches and cover extension only.

    This performs the standard search, which finds optimal values of
    the mismatches and cover extension parameters subject to the
    constraint on the total number of probes. It performs linear
    interpolation between calculated values of these two parameters
    for each dataset. It rounds these parameters to integers or, if
    desired, to the nearest value on an evenly spaced grid.

    Args:
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters
        max_total_count: upper bound on the number of total probes
        verify_without_interp: if True, check that the total probe count
            calculated without interpolation is the same as that calculated
            after rounding parameter values
        round_params: tuple (m, e); round mismatches to the nearest
            multiple of m and cover_extension to the nearest multiple
            of e (both m and e are int). If not set, default is
            (m, e) = (1, 1).
        loss_coeffs: tuple (m, e) giving the coefficients for the
            mismatches (m) and cover_extension (e) parameters in the
            loss function; if not set, default is (m, e) = (1, 1/100)
        dataset_weights: dict giving weight in the loss function for each
            dataset; if not set, default is a weight of 1 for each dataset

    Returns:
        tuple (x, y, z) where:
            x is a dict {dataset: p} where p is a tuple
                (mismatches, cover_extension) giving optimal values for
                those parameters
            y is the total number of probes required with the parameters in x
            z is the loss for the parameter values in x
    """
    # Set default values for arguments provided as None
    if loss_coeffs:
        # There should be a coefficient for each of the 2 parameters
        assert len(loss_coeffs) == 2
        loss_coeffs = tuple(loss_coeffs)
    else:
        loss_coeffs = (1.0, 1.0/100.0)
    if dataset_weights:
        # There should be a weight for each dataset
        for d in probe_counts.keys():
            assert d in dataset_weights
    else:
        dataset_weights = {d: 1.0 for d in probe_counts.keys()}
    if round_params:
        mismatches_round, cover_extension_round = round_params
    else:
        mismatches_round, cover_extension_round = 1, 1

    # Setup the loss function, parameter bounds, and make an initial guess
    loss_fn = _make_loss_fn(probe_counts, max_total_count, loss_coeffs,
        dataset_weights, interp_fn_type='standard')
    bounds = _make_param_bounds_standard(probe_counts)
    x0 = _make_initial_guess(probe_counts, bounds, 2)

    # Find the optimal parameter values, interpolating probe counts
    # for parameter values between what have been explicitly calculated
    x_sol = _optimize_loss(probe_counts, loss_fn, bounds, x0,
        interp_fn_type='standard')

    # Log the parameter values for each dataset, and the total probe count
    logger.info("##############################")
    logger.info("Continuous parameter values:")
    _log_params_by_dataset(x_sol, probe_counts, "float")
    x_sol_count = ic._make_total_probe_count_across_datasets_fn(
        probe_counts, interp_fn_type='standard')(x_sol)
    logger.info("TOTAL INTERPOLATED PROBE COUNT: %f", x_sol_count)
    logger.info("##############################")

    # Round the interpolated parameter values
    opt_params = _round_params(x_sol, probe_counts, max_total_count,
        loss_coeffs, dataset_weights,
        mismatches_round=mismatches_round,
        cover_extension_round=cover_extension_round)

    # Log the rounded parameter values, the total probe count, and the
    # loss on the rounded values
    logger.info("##############################")
    logger.info("Rounded parameter values:")
    _log_params_by_dataset(opt_params, probe_counts, "int")
    opt_params_count = ic._make_total_probe_count_across_datasets_fn(
        probe_counts, interp_fn_type='standard')(opt_params)
    opt_params_loss = loss_fn(opt_params, 0)
    logger.info("TOTAL PROBE COUNT: %d", opt_params_count)
    logger.info("TOTAL PARAMS LOSS: %f", opt_params_loss)
    logger.info("##############################")

    # Log the total probe count on the rounded parameter values without
    # using interpolation
    if verify_without_interp:
        logger.info("##############################")
        opt_params_count_no_interp = _total_probe_count_without_interp(opt_params, probe_counts)
        logger.info("TOTAL PROBE COUNT WITHOUT INTERP: %d", opt_params_count_no_interp)
        logger.info("##############################")
        # As a sanity check, verify that we get the same total probe count
        # without using interpolation
        assert opt_params_count == opt_params_count_no_interp

    opt_params_dict = {}
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches = opt_params[2 * i]
        cover_extension = opt_params[2 * i + 1]
        opt_params_dict[dataset] = (mismatches, cover_extension)

    return (opt_params_dict, opt_params_count, opt_params_loss)


def higher_dimensional_search(param_names, probe_counts, max_total_count,
        loss_coeffs=None, dataset_weights=None):
    """Search over multiple arbitrary parameters.

    Unlike the standard search, this can search over any number of
    provided parameters. It interpolates linearly using a function
    from scipy for unstructured data of an arbitrary dimension. Unlike
    the standard search, this does not round values (e.g., they may
    remain fractional even if the parameter is integral).

    Args:
        param_names: tuple giving names of parameters
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters (the tuple specifying
            a choice of parameter values is ordered such that parameters
            correspond to those in param_names)
        max_total_count: upper bound on the number of total probes
        loss_coeffs: coefficient to use for each parameter in the loss
            function. If not set, default is 1 for each parameter
        dataset_weights: dict giving weight in the loss function for each
            dataset; if not set, default is a weight of 1 for each dataset

    Returns:
        tuple (x, y, z) where:
            x is a dict {dataset: p} where p is a tuple giving optimal
                values for parameters, corresponding to the order in
                param_names
            y is the total number of probes required with the parameters in x
            z is the loss for the parameter values in x
    """
    num_params = len(param_names)

    # Set default values for arguments provided as None
    if loss_coeffs is None:
        # The default coefficient is 1 for each parameter
        loss_coeffs = tuple(1.0 for _ in range(num_params))
    else:
        # There must be a coefficient for each parameter
        assert len(loss_coeffs) == num_params
        loss_coeffs = tuple(loss_coeffs)
    if dataset_weights:
        # There should be a weight for each dataset
        for d in probe_counts.keys():
            assert d in dataset_weights
    else:
        dataset_weights = {d: 1.0 for d in probe_counts.keys()}

    # Setup the loss function, parameter bounds, and make an initial guess
    loss_fn = _make_loss_fn(probe_counts, max_total_count, loss_coeffs,
        dataset_weights, interp_fn_type='nd')
    x0 = _make_initial_guess(probe_counts, None, num_params)

    # Find the optimal parameter values, interpolating probe counts
    # for parameter values between what have been explicitly calculated
    bounds = _make_param_bounds_nd(probe_counts)
    x_sol = _optimize_loss(probe_counts, loss_fn, bounds, x0,
        interp_fn_type='nd')

    x_sol_dict = {}
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        x_sol_dict[dataset] = tuple(x_sol[num_params * i + j]
            for j in range(num_params))

    x_sol_count = ic._make_total_probe_count_across_datasets_fn(
        probe_counts, interp_fn_type='nd')(x_sol)
    x_sol_loss = loss_fn(x_sol, 0)

    return (x_sol_dict, x_sol_count, x_sol_loss)
