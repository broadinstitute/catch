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

from collections import defaultdict
import logging
import math

import numpy as np
from scipy import optimize

__author__ = 'Hayden Metsky <hayden@mit.edu>'

logger = logging.getLogger(__name__)


def _round_up(x, b):
    """Round float x up to the nearest multiple of int b
    """
    return int(math.ceil(float(x) / b)) * b
def _round_down(x, b):
    """Round float x down to the nearest multiple of int b
    """
    return int(math.floor(float(x) / b)) * b


def _make_interp_probe_count_for_dataset_fn(probe_counts):
    """Generate and return a function that interpolates probe count for a dataset.

    Args:
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters

    Returns:
        function whose input is a dataset, value for the mismatches parameter,
        and value for the cover extension parameter. The function linearly
        interpolates the number of probes required in that dataset for
        those parameter values, based on the values (which were explicitly
        calculated) in probe_counts.
    """
    memoized_bounding_boxes = {dataset: {} for dataset in probe_counts.keys()}
    def immediate_bounding_box(mismatches, cover_extension):
        """Return the box around the values rounded down/up to the nearest int
        to use when memoizing bounding boxes.
        """
        return (_round_down(mismatches, 1),
                _round_up(mismatches, 1),
                _round_down(cover_extension, 1),
                _round_up(cover_extension, 1))

    def find_bounding_box_around_point(dataset, mismatches, cover_extension):
        """Return a rectangular bounding box around given parameters.

        Based on probe_counts[dataset], choose parameter values that form the
        smallest box around (mismatches, cover_extension) so that we can
        use actual computed probe counts (based on the box's parameter values)
        to interpolate the probe count for (mismatches, cover_extension).
        """
        # 'mismatches' may be a float between the min/max mismatches for
        # dataset. 'cover_extension' may be a float between the min/max
        # cover_extension for dataset. Find the mismatches and cover_extension
        # parameters for dataset that form the smallest rectangle encompassing
        # 'mismatches' and 'cover_extension', so that we can use the points of
        # this rectangle to perform bilinear interpolation. We brute force over
        # all possible rectangles encompassing 'mismatches' and
        # 'cover_extension' (with some heuristics to not be completely
        # exhaustive), which takes O(n^3) time; there are faster methods, but
        # this should suffice here.
        # Consider mismatches on the x-axis (left to right) and cover_extension
        # on the y-axis (bottom to top)
        points = set(probe_counts[dataset].keys())
        points_topleft = set()
        points_topright = set()
        points_bottomleft = set()
        points_bottomright = set()
        for p in points:
            m, ce = p
            if m == mismatches:
                if ce == cover_extension:
                    points_topleft.add(p)
                    points_topright.add(p)
                    points_bottomleft.add(p)
                    points_bottomright.add(p)
                elif ce > cover_extension:
                    points_topleft.add(p)
                    points_topright.add(p)
                else:
                    points_bottomleft.add(p)
                    points_bottomright.add(p)
            elif m > mismatches:
                if ce == cover_extension:
                    points_topright.add(p)
                    points_bottomright.add(p)
                elif ce > cover_extension:
                    points_topright.add(p)
                else:
                    points_bottomright.add(p)
            else:
                if ce == cover_extension:
                    points_topleft.add(p)
                    points_bottomleft.add(p)
                elif ce > cover_extension:
                    points_topleft.add(p)
                else:
                    points_bottomleft.add(p)

        points_topright_by_y = defaultdict(set)
        for p in points_topright:
            m, ce = p
            points_topright_by_y[ce].add(p)
        points_bottomleft_by_x = defaultdict(set)
        for p in points_bottomleft:
            m, ce = p
            points_bottomleft_by_x[m].add(p)

        min_rectangle, min_area = None, float('inf')
        for p_topleft in points_topleft:
            p_topleft_m, p_topleft_ce = p_topleft
            # find all points in the top-right with the same cover_extension
            # (y-value) as p_topleft
            for p_topright in points_topright_by_y[p_topleft_ce]:
                p_topright_m, p_topright_ce = p_topright
                # find all points in the bottom-left with the same mismatches
                # (x-value) as p_topleft
                for p_bottomleft in points_bottomleft_by_x[p_topleft_m]:
                    p_bottomleft_m, p_bottomleft_ce = p_bottomleft
                    # to form a rectangle, we need the point (p_topright_m,
                    # p_bottomleft_ce); this point should be in the
                    # bottom-right, so check if it exists
                    p_bottomright = (p_topright_m, p_bottomleft_ce)
                    if p_bottomright in points_bottomright:
                        # we found a valid rectangle; now compute its 'area'
                        width = p_topright_m - p_topleft_m
                        height = (p_topright_ce - p_bottomleft_ce) / 10.0
                        # add pseudocounts to width and height because if
                        # width or height is 0, we still want the other
                        # dimension to be accounted for
                        area = (width + 0.001) * (height + 0.001)
                        if area < min_area:
                            min_rectangle = (p_topleft, p_bottomright)
                            min_area = area
        return min_rectangle

    def interp_probe_count_for_dataset(dataset, mismatches, cover_extension):
        """
        Using the given probe counts at particular parameter values, interpolate
        the number of probes for 'dataset' and 'mismatches' mismatches and
        at a cover extension of 'cover_extension', where each of these may be
        floats
        """
        immediate_bb = immediate_bounding_box(mismatches, cover_extension)
        if immediate_bb in memoized_bounding_boxes[dataset]:
            # The bounding box for (mismatches, cover_extension) has been
            # memoized
            min_rectangle = memoized_bounding_boxes[dataset][immediate_bb]
        else:
            # Compute and memoize the bounding box for (mismatches,
            # cover_extension)
            min_rectangle = find_bounding_box_around_point(dataset,
                                                           mismatches,
                                                           cover_extension)
            if min_rectangle is None:
                raise ValueError(("Unable to find rectangular bounding box around "
                                  "(mismatches, cover_extension)=(%f, %f) for "
                                  "dataset %s") % (mismatches, cover_extension,
                                  dataset))
            memoized_bounding_boxes[dataset][immediate_bb] = min_rectangle

        rect_topleft, rect_bottomright = min_rectangle
        mismatches_floor, cover_extension_ceil = rect_topleft
        mismatches_ceil, cover_extension_floor = rect_bottomright

        # At cover_extension_floor and at cover_extension_ceil, interpolate the
        # number of probes at 'mismatches' mismatches (interpolate linearly)
        for ce in [cover_extension_floor, cover_extension_ceil]:
            count_left = probe_counts[dataset][(mismatches_floor, ce)]
            count_right = probe_counts[dataset][(mismatches_ceil, ce)]
            mismatches_diff = mismatches_ceil - mismatches_floor
            if mismatches_diff == 0:
                # count_left should equal count_right
                assert count_left == count_right
                count = count_left
            elif count_left <= count_right:
                count_diff = count_right - count_left
                f = float(mismatches - mismatches_floor) / mismatches_diff
                count = f * count_diff + count_left
            else:
                count_diff = count_left - count_right
                f = float(mismatches - mismatches_floor) / mismatches_diff
                count = count_left - f * count_diff
            if ce == cover_extension_floor:
                count_floor = count
            if ce == cover_extension_ceil:
                count_ceil = count

        # Interpolate between cover_extension_floor and cover_extension_ceil
        # using count_floor and count_ceil (interpolate linearly)
        cover_extension_diff = cover_extension_ceil - cover_extension_floor
        if cover_extension_diff == 0:
            # count_floor should equal count_ceil
            assert count_floor == count_ceil
            final_interp = count_floor
        elif count_floor <= count_ceil:
            count_diff = count_ceil - count_floor
            f = float(cover_extension - cover_extension_floor) / cover_extension_diff
            final_interp = f * count_diff + count_floor
        else:
            count_diff = count_floor - count_ceil
            f = float(cover_extension - cover_extension_floor) / cover_extension_diff
            final_interp = count_floor - f * count_diff

        return final_interp

    return interp_probe_count_for_dataset


def _make_total_probe_count_across_datasets_fn(probe_counts):
    """Generate and return a function that interpolates probe count.

    Args:
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters

    Returns:
        function whose input is choices of mismatches and cover extension
        parameter values for all datasets. The function interpolates the
        number of probes required for each dataset, and sums across all
        datasets to determine a total number of probes.
    """
    interp_probe_count_for_dataset = _make_interp_probe_count_for_dataset_fn(
        probe_counts)

    def total_probe_count_across_datasets(x):
        """
        Sum the (interpolated) probe counts across datasets.

        x is a list giving all the parameter values across datasets,
        such that for even i, x_i gives the number of mismatches for the
        (i/2)'th dataset and x_{i+1} gives the cover extension for the
        (i/2)'th dataset
        """
        s = 0
        for i, dataset in enumerate(sorted(probe_counts.keys())):
            mismatches, cover_extension = x[2 * i], x[2 * i + 1]
            s += interp_probe_count_for_dataset(dataset, mismatches,
                                                cover_extension)
        return s

    return total_probe_count_across_datasets


def _make_loss_fn(probe_counts, max_total_count):
    """Generate and return a loss function.

    The function calculates a loss over the parameters and adds onto
    that loss to meet a constraint based on the barrier method. It
    uses a logarithmic barrier function to enforce the constraint that
    the total probe count be <= max_total_count.

    The loss over the parameters is the sum, over datasets d, of
        (m_d)^2 + (e_d / 10.0)^2
    where m_d is the value of the mismatches parameter for d and e_d
    is the value of the cover extension parameter for d.

    Args:
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters
        max_total_count: upper bound on the number of total probes

    Returns:
        a function that is the sum of a loss defined over the parameters
        and a value designed to enforce a barrier on the total number
        of probes
    """
    total_probe_count_across_datasets = _make_total_probe_count_across_datasets_fn(
        probe_counts)

    def loss(x, *func_args):
        """
        Compute a loss.

        x is a list giving all the parameter values across datasets,
        such that for even i, x_i gives the number of mismatches for the
        (i/2)'th dataset and x_{i+1} gives the cover extension for the
        (i/2)'th dataset
        """
        # First compute a loss over the parameters by taking their L2-norm (and
        # down-weighting cover_extension by a factor of 10.0)
        # This is the function we really want to minimize
        opt_val = 0
        for i, dataset in enumerate(sorted(probe_counts.keys())):
            mismatches, cover_extension = x[2 * i], x[2 * i + 1]
            opt_val += np.power(mismatches, 2.0) + np.power(cover_extension / 10.0, 2.0)

        # We also have the constraint that the total probe count be less than
        # max_total_count
        # We add a barrier function to enforce this constraint and weight the
        # barrier by eps
        eps = func_args[0]
        total_probe_count = total_probe_count_across_datasets(x)
        if total_probe_count >= max_total_count:
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


def _make_param_bounds(probe_counts, step_size=0.001):
    """Calculate bounds on parameter values.

    For each dataset d, this calculates bounds on the values of the
    mismatches and cover extension parameters based on which values
    have a number of probes calculated for them. Namely, we wish to
    ensure we can find a bounding box around an arbitrary point
    (mismatches, cover_extension); based on the values in probe_counts,
    the bounds should ensure this.

    Args:
        probe_counts: dict giving the number of probes for each dataset
            and choice of parameters
        step_size: small value subtracted from each upper bound

    Returns:
        [m_1, e_1, m_2, e_2, ...] where m_i is a tuple (lo, hi) that
        gives bounds on the number of mismatches for the i'th dataset,
        and likewise for e_i
    """
    bounds = []
    for dataset in sorted(probe_counts.keys()):
        params = probe_counts[dataset].keys()

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


def _make_initial_guess(probe_counts, bounds, max_total_count):
    """Make initial guess for optimal parameter values.

    This guesses a value for mismatches and cover extension separately
    for each dataset, choosing uniformly from within the provided
    bounds. It allows a guess that exceeds the constraint on the
    total number of probes, but outputs a warning if this occurs.

    Args:
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters
        bounds: lower/upper bounds on mismatches and cover extension
            for each dataset
        max_total_count: upper bound on the number of total probes

    Returns:
        list of length 2*(number of datasets) in which the value at
        index 2*i is a guess for the mismatches of the i'th dataset
        and the value at index 2*i+1 is a guess for the cover
        extension of the i'th dataset
    """
    # Guess uniformly within bounds for each dataset
    x0 = np.zeros(2 * len(probe_counts))
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches_lo, mismatches_hi = bounds[2 * i]
        x0[2 * i] = np.random.uniform(mismatches_lo, mismatches_hi)
        cover_extension_lo, cover_extension_hi = bounds[2 * i + 1]
        x0[2 * i + 1] = np.random.uniform(cover_extension_lo,
                                          cover_extension_hi)

    # Verify that this yields fewer probes than the maximum allowed
    # (i.e., is not beyond the barrier)
    guess_probe_count = _make_total_probe_count_across_datasets_fn(probe_counts)(x0)
    if guess_probe_count >= max_total_count:
        logger.warning(("WARNING: Initial guess is beyond the probe barrier "
                "(%d, but the max is %d)"), guess_probe_count, max_total_count)
        logger.warning("         ...continuing anyway")

    return x0


def _optimize_loss(probe_counts, loss_fn, bounds, x0,
                   initial_eps=10.0, step_size=0.001):
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
        bounds: bounds on the parameter values provided by _make_param_bounds
        x0: the initial guess of parameter values (i.e., starting position)
        initial_eps: weight of the barrier on the first iteration
        step_size: epsilon value provided to optimize.fmin_tnc

    Returns:
        list of length 2*(number of datasets) where the value at index
        2*i is the mismatches parameter for the i'th dataset and the
        value at 2*i+1 is the cover extension parameter for the i'th
        dataset
    """
    eps = initial_eps
    while eps >= 0.01:
        x0_probe_count = _make_total_probe_count_across_datasets_fn(probe_counts)(x0)
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

    The result of _make_total_probe_count_across_datasets_fn should give
    the same count as this function (if params are keys in the datasets
    of probe_counts). But this uses probe_counts directly and can be
    used as a sanity check -- i.e., it does not do any interpolation.

    Args:
        params: parameter values to use when determining probe counts;
            params[2*i] is the number of mismatches of the i'th dataset
            and params[2*i+1] is the cover extension of the i'th dataset
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters

    Returns:
        total number of probes across all datasets, according to the
        given values of params
    """
    s = 0
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]
        s += probe_counts[dataset][(mismatches, cover_extension)]
    return s


def _round_params(params, probe_counts, max_total_count,
        mismatches_eps=0.01, cover_extension_eps=0.1):
    """Round parameter values while satisfying the constraint on total count.

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
        mismatches_epis/cover_extension_eps: eps as defined above for
            mismatches and cover_extension

    Returns:
        list in which index i corresponds to the parameter given in
        params[i], but rounded
    """

    params_rounded = []
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]

        if mismatches - _round_down(mismatches, 1) < mismatches_eps:
            # Round mismatches down
            mismatches = _round_down(mismatches, 1)
        else:
            # Round mismatches up
            mismatches = _round_up(mismatches, 1)

        if cover_extension - _round_down(cover_extension, 10) < cover_extension_eps:
            # Round cover_extension down
            cover_extension = _round_down(cover_extension, 10)
        else:
            # Round cover_extension up
            cover_extension = _round_up(cover_extension, 10)

        params_rounded += [mismatches, cover_extension]

    total_probe_count = _make_total_probe_count_across_datasets_fn(probe_counts)
    # Verify that the probe count satisfies the constraint
    # Note that this assertion may fail if we are dealing with datasets
    # for which few actual probe counts have been computed; in these
    # cases, the interpolation may severely underestimate the number
    # of probes at a particular parameter choice
    assert total_probe_count(params_rounded) < max_total_count

    # Keep decreasing parameters while satisfying the constraint.
    # In particular, choose to decrease the parameter whose reduction
    # yields the smallest loss while still satisfying the constraint.
    loss_fn = _make_loss_fn(probe_counts, max_total_count)
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
                # This is a mismatch, so decrease by 1
                params_tmp[i] -= 1
            else:
                # This is a cover_extension, so decrease by 10
                params_tmp[i] -= 10
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


def standard_search(probe_counts, max_total_count, verify_without_interp=False):
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

    Returns:
        tuple (x, y, z) where:
            x is a dict {dataset: p} where p is a tuple
                (mismatches, cover_extension) giving optimal values for
                those parameters
            y is the total number of probes required with the parameters in x
            z is the loss for the parameter values in x
    """
    # Setup the loss function, parameter bounds, and make an initial guess
    loss_fn = _make_loss_fn(probe_counts, max_total_count)
    bounds = _make_param_bounds(probe_counts)
    x0 = _make_initial_guess(probe_counts, bounds, max_total_count)

    # Find the optimal parameter values, interpolating probe counts
    # for parameter values between what have been explicitly calculated
    x_sol = _optimize_loss(probe_counts, loss_fn, bounds, x0)

    # Log the parameter values for each dataset, and the total probe count
    logger.info("##############################")
    logger.info("Continuous parameter values:")
    _log_params_by_dataset(x_sol, probe_counts, "float")
    x_sol_count = _make_total_probe_count_across_datasets_fn(probe_counts)(x_sol)
    logger.info("TOTAL INTERPOLATED PROBE COUNT: %f", x_sol_count)
    logger.info("##############################")

    # Round the interpolated parameter values
    opt_params = _round_params(x_sol, probe_counts, max_total_count)

    # Log the rounded parameter values, the total probe count, and the
    # loss on the rounded values
    logger.info("##############################")
    logger.info("Rounded parameter values:")
    _log_params_by_dataset(opt_params, probe_counts, "int")
    opt_params_count = _make_total_probe_count_across_datasets_fn(probe_counts)(opt_params)
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
