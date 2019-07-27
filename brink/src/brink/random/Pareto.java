package brink.random;

import brink.random.Uniform;

/**
 * Generates pareto-distributed random numbers.
 */
public class Pareto {

    /**
     * Generate a random number following the pareto distribution.
     */
    public static double pareto(double shape, double scale) {
        // https://en.wikipedia.org/wiki/Pareto_distribution#Random_sample_generation
        // U in (0,1]
        final double U = Uniform.uniform_exclude_lower(0, 1);
        return scale / Math.pow(U, 1/shape);
    }

    /**
     * Generate a random number following the pareto distribution, where the resultant
     * number will be between `min` and `max`.
     */
    public static double truncated_pareto(double shape, double min, double max) {
        // Unfortunately, the truncated pareto distribution equation can't deal with
        // a lower bound of 0, so we just shift both bounds up by one, and then reverse
        // it at the end of the function
        final boolean adjust_min = min == 0;
        if(adjust_min) {
            min++;
            max++;
        }
        // https://en.wikipedia.org/wiki/Pareto_distribution#Bounded_Pareto_distribution
        // U in (0,1)
        final double U = Uniform.uniform_exclude_both(0, 1);
        final double shaped_min = Math.pow(min, shape);
        final double shaped_max = Math.pow(max, shape);
        final double neg_inv_shape = -(1/shape);

        return Math.pow(-(U*shaped_max - U*shaped_min - shaped_max)/(shaped_min*shaped_max), neg_inv_shape) - (adjust_min ? 1 : 0);
    }
}
