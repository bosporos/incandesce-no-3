package brink.random;

import brink.Brink;

/**
 * For generating uniformly-distributed random numbers.
 */
public class Uniform {

    /**
     * Generate a random integer between `min` and `max`.
     */
    public static int uniform(int min, int max) {
        return (int)Brink.applet_instance.random(min, max);
    }

    /**
     * Generate a random double between `min` and `max`, inclusive.
     * I.E. this function will generate a random number R where R in [a,b]
     */
    public static double uniform(double min, double max) {
        return (double)Brink.applet_instance.random((float)min, (float)max);
    }

    /**
     * Generate a random number between `min` (exclusive) and `max` (inclusive).
     * I.E. this function will generate a random number R where R in (a,b]
     */
    public static double uniform_exclude_lower(double min, double max) {
        double ret;
        do {
            ret = Uniform.uniform(min, max);
        } while(ret == min);
        return ret;
    }

    /**
     * Generate a random number between `min` (inclusive) and `max` (exclusive).
     * I.E. this function will generate a random number R where R in [a,b)
     */
    public static double uniform_exclude_upper(double min, double max) {
        double ret;
        do {
            ret = Uniform.uniform(min, max);
        } while(ret == max);
        return ret;
    }

    /**
     * Generate a random numebr between `min` and `max`, exclusive
     * I.E. this function will generate a random number R where R in (a,b)
     */
    public static double uniform_exclude_both(double min, double max) {
        double ret;
        do {
            ret = Uniform.uniform(min, max);
        } while(ret == min || ret == max);
        return ret;
    }
}
