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
from scipy import interpolate
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


def _make_interp_probe_count_for_dataset_standard_fn(probe_counts):
    """Generate and return a function that interpolates probe count for a dataset.

    This operates only on the mismatches and cover_extension parameters.

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

    def interp_probe_count_for_dataset(dataset, param_vals):
        """
        Using the given probe counts at particular parameter values, interpolate
        the number of probes for 'dataset' and given mismatches (param_vals[0])
        and cover_extension (param_vals[1]), where each of these may be floats
        """
        mismatches, cover_extension = param_vals

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


_interp_nd_fn_memoized = {}
def _make_interp_probe_count_for_dataset_nd_fn(probe_counts):
    """Generate and return a function that interpolates probe count for a dataset.

    This uses a function from scipy's interpolate package to operate on
    an arbitrary number of parameters.

    Args:
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters

    Returns:
        function whose input is a dataset and values for arbitrary
        parameters. The function linearly interpolates the number of
        probes required in that dataset for those parameter values, based
        on the values (which were explicitly calculated) in probe_counts.
    """
    # Reset the memoized dict for a new call (this is useful for unit tests,
    # which may call this function multiple times with different inputs --
    # i.e., values in probe_counts)
    _interp_nd_fn_memoized = {}

    def interp_probe_count_for_dataset(dataset, param_vals):
        """
        Using the given probe counts at particular parameter values, interpolate
        the number of probes for 'dataset' and the parameter values given
        in 'param_vals', where each of these may be floats.
        """
        if dataset in _interp_nd_fn_memoized:
            nd_fn = _interp_nd_fn_memoized[dataset]
        else:
            points = []
            values = []
            for p in probe_counts[dataset].keys():
                points += [p]
                values += [probe_counts[dataset][p]]
            points = np.array(points)
            values = np.array(values)

            nd_fn = interpolate.LinearNDInterpolator(points, values,
                rescale=True)
            _interp_nd_fn_memoized[dataset] = nd_fn

        try:
            return nd_fn(np.array(param_vals))[0]
        except ValueError:
            raise ValueError(param_vals, dataset, probe_counts[dataset])

    return interp_probe_count_for_dataset


def _make_total_probe_count_across_datasets_fn(probe_counts,
        interp_fn_type='standard'):
    """Generate and return a function that interpolates probe count.

    Args:
        probe_counts: dict giving number of probes for each dataset and
            choice of parameters
        interp_fn_type: 'standard' (only perform interpolation on mismatches
            and cover_extension parameters) or 'nd' (use scipy's interpolate
            package to interpolate over n-dimensions)

    Returns:
        function whose input is choices of parameter values for all datasets.
        The function interpolates the number of probes required for each
        dataset, and sums across all datasets to determine a total number
        of probes.
    """
    assert interp_fn_type in ['standard', 'nd']
    if interp_fn_type == 'standard':
        interp_fn = _make_interp_probe_count_for_dataset_standard_fn
    elif interp_fn_type == 'nd':
        interp_fn = _make_interp_probe_count_for_dataset_nd_fn
    interp_probe_count_for_dataset = interp_fn(probe_counts)

    def total_probe_count_across_datasets(x):
        """
        Sum the (interpolated) probe counts across datasets.

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

        s = 0
        for i, dataset in enumerate(sorted(probe_counts.keys())):
            param_vals = [x[num_params * i + j] for j in range(num_params)]
            s += interp_probe_count_for_dataset(dataset, param_vals)
        return s

    return total_probe_count_across_datasets


def _make_loss_fn(probe_counts, max_total_count, interp_fn_type='standard'):
    """Generate and return a loss function.

    The function calculates a loss over the parameters and adds onto
    that loss to meet a constraint based on the barrier method. It
    uses a logarithmic barrier function to enforce the constraint that
    the total probe count be <= max_total_count.

    The loss over the parameters is the sum, over datasets d, of the
    sum of (v_{di})^2 for each parameter i in dataset d.

    Args:
        probe_counts: dict giving number of probes required for each
            dataset and choice of parameters
        max_total_count: upper bound on the number of total probes
        interp_fn_type: 'standard' (only perform interpolation on mismatches
            and cover_extension parameters) or 'nd' (use scipy's interpolate
            package to interpolate over n-dimensions)

    Returns:
        a function that is the sum of a loss defined over the parameters
        and a value designed to enforce a barrier on the total number
        of probes
    """
    total_probe_count_across_datasets = _make_total_probe_count_across_datasets_fn(
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

        # First compute a loss over the parameters by taking their L2-norm (and
        # down-weighting cover_extension by a factor of 10.0)
        # This is the function we really want to minimize
        opt_val = 0
        for i, dataset in enumerate(sorted(probe_counts.keys())):
            param_vals = [x[num_params * i + j] for j in range(num_params)]
            # Assume (perhaps incorrectly) that if there are just 2
            # parameters, these are mismatches and cover_extension; this
            # will be fixed later by allowing custom coefficients
            if num_params == 2:
                mismatches, cover_extension = x[2 * i], x[2 * i + 1]
                opt_val += np.power(mismatches, 2.0) + np.power(cover_extension / 10.0, 2.0)
            else:
                for v in param_vals:
                    opt_val += np.power(v, 2.0)

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
        x0_probe_count = _make_total_probe_count_across_datasets_fn(
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

    The result of _make_total_probe_count_across_datasets_fn should give
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


def _round_params(params, probe_counts, max_total_count,
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

    params_rounded = []
    for i, dataset in enumerate(sorted(probe_counts.keys())):
        mismatches, cover_extension = params[2 * i], params[2 * i + 1]

        if mismatches - _round_down(mismatches, mismatches_round) < mismatches_eps:
            # Round mismatches down
            mismatches = _round_down(mismatches, mismatches_round)
        else:
            # Round mismatches up
            mismatches = _round_up(mismatches, mismatches_round)

        if cover_extension - _round_down(cover_extension, cover_extension_round) < cover_extension_eps:
            # Round cover_extension down
            cover_extension = _round_down(cover_extension, cover_extension_round)
        else:
            # Round cover_extension up
            cover_extension = _round_up(cover_extension, cover_extension_round)

        params_rounded += [mismatches, cover_extension]

    total_probe_count = _make_total_probe_count_across_datasets_fn(
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
    loss_fn = _make_loss_fn(probe_counts, max_total_count,
        interp_fn_type='standard')
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
        verify_without_interp=False, mismatches_round=1,
        cover_extension_round=1):
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
        mismatches_round/cover_extension_round: round mismatches and
            cover_extension parameter to the nearest multiple of this
            int

    Returns:
        tuple (x, y, z) where:
            x is a dict {dataset: p} where p is a tuple
                (mismatches, cover_extension) giving optimal values for
                those parameters
            y is the total number of probes required with the parameters in x
            z is the loss for the parameter values in x
    """
    # Setup the loss function, parameter bounds, and make an initial guess
    loss_fn = _make_loss_fn(probe_counts, max_total_count,
        interp_fn_type='standard')
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
    x_sol_count = _make_total_probe_count_across_datasets_fn(
        probe_counts, interp_fn_type='standard')(x_sol)
    logger.info("TOTAL INTERPOLATED PROBE COUNT: %f", x_sol_count)
    logger.info("##############################")

    # Round the interpolated parameter values
    opt_params = _round_params(x_sol, probe_counts, max_total_count,
        mismatches_round=mismatches_round,
        cover_extension_round=cover_extension_round)

    # Log the rounded parameter values, the total probe count, and the
    # loss on the rounded values
    logger.info("##############################")
    logger.info("Rounded parameter values:")
    _log_params_by_dataset(opt_params, probe_counts, "int")
    opt_params_count = _make_total_probe_count_across_datasets_fn(
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


def higher_dimensional_search(param_names, probe_counts, max_total_count):
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

    Returns:
        tuple (x, y, z) where:
            x is a dict {dataset: p} where p is a tuple giving optimal
                values for parameters, corresponding to the order in
                param_names
            y is the total number of probes required with the parameters in x
            z is the loss for the parameter values in x
    """
    num_params = len(param_names)

    # Setup the loss function, parameter bounds, and make an initial guess
    loss_fn = _make_loss_fn(probe_counts, max_total_count,
        interp_fn_type='nd')
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

    x_sol_count = _make_total_probe_count_across_datasets_fn(
        probe_counts, interp_fn_type='nd')(x_sol)
    x_sol_loss = loss_fn(x_sol, 0)

    return (x_sol_dict, x_sol_count, x_sol_loss)
