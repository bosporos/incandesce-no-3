package brink.random;

import brink.Brink;

/**
 * Generates normally distributed random numbers.
 */
public class Normal {

    /**
     * Generates random numbers along a Gaussian distribution
     */
    public static double gaussian() {
        return (double)Brink.applet_instance.randomGaussian();
    }

    /**
     * Generates random numbers along a normal distribution centered at `mu` with
     * standard deviation `sigma`.
     */
    public static double normal(double mu, double sigma) {
        return mu + sigma * Normal.gaussian();
    }

    /**
     * Generates random numbers between `lower` and `upper`, normally distributed
     * with mean `mu` and standard deviation `sigma`.
     */
    public static double truncated_normal(double mu,
                                         double sigma,
                                         double lower,
                                         double upper) {
        // I'm sure that there's a better method than a brute-force rejection algorithm like this, but I don't have the math for it.
        double deviate;
        do {
            deviate = Normal.normal(mu, sigma);
        } while(deviate < lower || deviate > upper);
        return deviate;
    }
}
