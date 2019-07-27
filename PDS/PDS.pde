/* NAME: Maximilien Angelo M. Cura
 * DATE: 16 July 2019
 * DESC: I was inspired in part by the writing system of the aliens in the movie Arrival, and I wanted to do something with Perlin
         noise fields, which I've seen here and there but never had the chance to play around with.

         This artwork generates a Poisson disk sampling in a 4-dimensional hyperrectangle, folds its lower dimensions into a torus,
         and then flattens the torus onto the canvas; around each point in the now-torioidal sampling a number of particles, depenedent
         upon the point's position in the fourth dimension of the original hyperrectangular sampling. Each particle is then moved across
         a fixed-magnitude velocity field generated from Perlin noise, using a polar transform to direct the particles movement away
         from the center of the circle, and its path is traced as a pixel exposure whose coloring reflects its proximity to the points
         generated in a second Poisson disk sampling on the 2-dimensional rectangle of the screen.
 */

import brink.Brink;
import brink.random.Uniform;
import brink.random.Pareto;

 /**
  * Processing does not supply the natural number as a constant.
  */
 public static final float E = 2.71828182846;

 /**
  * Canvas width and height, in pixels
  */
 public static final float kCanvasWidth = 1000;
 public static final float kCanvasHeight = 1000;

 /**
  * Fallaway parameter for noiseDetail()
  */
 public static final float kPerlinFallaway = 0.5;

 /**
  */
 public static final float kPerlinOctaves = 4;

 /**
  * Multiplier for the parameters passed to the noise() function, to ensure the appropriate degree of smoothness
  * in the Perlin noise field.
  */
 public static final float kPerlinNoiseSmoothingFactor = 0.00325;

 // /**
 //  * Unfortunate processing.js shim to correct the values coming out of noise(x,y)
 //  * @internal magic constants; values are seemingly arbitrary.
 //  */
 // public final float kPJSShimPerlinInMinimum = 0.26816919850543997;
 // public final float kPJSShimPerlinInMaximum = 0.6814498435318399;
 // public final float kPJSShimPerlinOutMinimum = 0;
 // public final float kPJSShimPerlinOutMaximum = 1;

 public static final int kDefaultLineColor = 0x08FBFAF8;
 public static final int kAccentLineColor = 0x08DB5461;

 public static final int kParticleSystemSteps = 500;

 /**
  * Used in the proximity probability calculation that determines the coloring of each particle's pixel  exposure
  */
 public final float kColoringFalloffDistance = E * 20;

 /**
  * Turn debugging diagnostics on/off.
  */
 public static final boolean kDiagnostics = true;

 /**
  * The 4-dimensional space that the first Poisson disk sampling takes place in.
  */
 Space four_space = new Space(4);

 /**
  * The 2-dimensional space used for all the Vectors that exist in the plane of the canvas
  */
 Space two_space = new Space(2);

 /**
  * The Poisson disk sampling used for generating the particles; conceptually, it's a 4-orthotope (4-dimensional hyperrectangle),
  * whose lower three dimensions are mapped to a 3-dimensional ring torus; each point in the sampling may be described by a vector
  * <r, a, b, c> where r, a describe a point, in polar coordinates, on a cross section of a ring torus with a radius of 250 pixels,
  * b describes the position of the aforementioned cross section along the torus, and c is involved in the calculation of the number
  * of particles that will be generated around that point.
  */
 PoissonDiskSample particle_generation_pds;

 /**
  * The Poisson disk sampling used for determining the color of the particles' pixel exposures; if the point around which the particles
  * were generated was within a circle centred on any of the points in this sampling, with a radius randomly distributed according to
  * the Pareto distribution, then the color of the particle exposures will be red, otherwise white.
  */
 PoissonDiskSample color_determination_pds;

void settings() {
    size((int) kCanvasWidth, (int) kCanvasHeight, P2D);
}

void setup() {
    // Setup for random number classes
    Brink.init(this);

    // Dimensions of the 4-orthotope used for the computation of `particle_generation_pds`
    Vector spatial_dimensions = new Vector(four_space, TWO_PI, TWO_PI, PI * 30, PI);
    // Because all the shapes are center-defined, we have to use this.
    Vector spatial_dimension_center = spatial_dimensions.div(2);

    // Create the Poisson disk samples
    particle_generation_pds = new PoissonDiskSample(
        new Hyperrectangle(four_space, spatial_dimension_center, spatial_dimensions),
        /* In 4d */
        four_space,
        /* Distance between samples; found by trial and error. */
        PI,
        /* Use a high A-value to trim down the number of subdivisions; improves performance */
        32);
    color_determination_pds = new PoissonDiskSample(
        /* The screen */
        new Hyperrectangle(two_space,
            new Vector(two_space, kCanvasWidth / 2, kCanvasHeight / 2),
            new Vector(two_space, kCanvasWidth, kCanvasHeight)),
            two_space,
            /* Larger distance between points so that the whole image isn't just flooded with red*/
            333,
            /* Lower A-value because the PDS for this particular is much simpler */
            1);

    colorMode(RGB, 0xFF, 0xFF, 0xFF, 0xFF);

    if(kDiagnostics) {
        println("== TESTING TRANSFORM SANITY ==");
        Vector scaffold = new Vector(four_space, 2, 2, 3, 2);
        Vector r = new Vector(four_space).from_linear_index(20, scaffold);
        println("Vector (FLI): " +  r);
        println("Vector (LIT): " + r.linear_index_transform(scaffold));
    }

    // Set up the noise properly
    noiseDetail ((int) kPerlinOctaves, kPerlinFallaway);

    noLoop();
}

void draw() {
    // Reset the noise seed each run
    noiseSeed ((int)random (0, MAX_FLOAT));

    if(kDiagnostics) {
        println ("== BEGINNING PDS SAMPLE COMPUTATION ==");
        println ("[INCANDESCE] ==> BEGINNING COMPUTATION OF HYPERparticle_generation_pdsAL PDS SAMPLE");
    }

    // Run the PDS computations

    particle_generation_pds.begin_computation ();
    int iter = 0;
    do {
        // Perform an iteration
        particle_generation_pds.iterate ();
        if(kDiagnostics)
            println (">>>>>>>>>>>>>>>> ITERATION " + iter + " WITH ACTIVE CELLS: " + particle_generation_pds.actives.size ());
        iter++;
        // Cap at 6 iterations to ensure that it actually finishes
    } while (iter < 6 && !particle_generation_pds.computation_finished ());
    if(kDiagnostics)
        println ("            ---> FINISHED");
    ArrayList<Vector> points = particle_generation_pds.get_samples ();
    if(kDiagnostics)
        println ("\t\tSample size: " + points.size ());

    if(kDiagnostics)
        println ("[INCANDESCE] ==> BEGINNING COMPUTATION OF RECTANGULAR PDS SAMPLE");

    // Create the color plane

    color_determination_pds.begin_computation ();
    iter = 0;
    do {
        color_determination_pds.iterate ();
        if(kDiagnostics)
            println (">>>>>>>>>>>>>>>> ITERATION " + iter);
        iter++;
        // Again, cap at 6; usually runs at ~3~4 iterations, though
    } while (iter < 6 && !color_determination_pds.computation_finished ());
    ArrayList<Vector> color_anchors = color_determination_pds.get_samples ();
    if(kDiagnostics) {
        println ("            ---> FINISHED");
        println ("== FINISHED PDS SAMPLE COMPUTATION ==");
    }

    if(kDiagnostics)
        println("== BEGINNING PARTICLE GENERATION ==");
    // Brief note about array sizes: the arrays *will* all be too large; there *will* be significant underused capacity.
    // However, the amount of array capacity that is actually utilised is nondeterministic, so we use the upper bound on
    // the number of particles per point to inform the sizes of the arrays

    // The positions of all the particles; for simplicity's sake, they do not carry velocity or acceleration, only position
    Vector particles[] = new Vector[points.size () * 16];
    // The colors of the particles
    int colors[]       = new int[points.size () * 16];
    // The "base angle" of each particle, i.e. the angle, from the horizontal, of the point that originated the particle
    float angles[]     = new float[points.size () * 16];
    // Because number of particles is nondeterministic, we have to do this tracking externally to the iteration counter(s)
    int total_particles = 0;
    int particle_index  = 0;
    // Diagnostic values
    double ran = PI, rax = PI;
    double rsan = PI, rsax = PI;
    double rsdn = PI, rsdx = PI;
    double rcn = PI, rcx = PI;
    for (int i = 0; i < points.size (); i++) {
        // The position of the cross section in the torus, defined as the angle that it occurs at relative the origin point
        float angle         = (float)points.get (i).tuple[3];
        // Polar coordinates of a position w/in that cross section
        float scatter_angle = (float)points.get (i).tuple[2];
        float scatter_dist  = (float)points.get (i).tuple[1];
        // Number of particles to be generated around a point
        // Biased towards higher numbers
        int num_children    = (int) round(16 - (float) points.get(i).tuple[0] * log((float) pow((float) points.get(i).tuple[0], E / 2)));

        // Diagnostics
        if (angle < ran) ran = angle;
        if (angle > rax) rax = angle;

        if (scatter_angle < rsan) rsan = scatter_angle;
        if (scatter_angle > rsax) rsax = scatter_angle;

        if (scatter_dist < rsdn) rsdn = scatter_dist;
        if (scatter_dist > rsdx) rsdx = scatter_dist;

        if (num_children < rcn) rcn = num_children;
        if (num_children > rcx) rcx = num_children;

        // This is the actual coordinate of the point translated onto the scren
        // Basically, it goes something like:
        // x = Bx + rcos(a) + Sr*cos(a + cos(Sa)^2)
        // y = By + rsin(a) + Sr*sin(a + sin(Sa)^2)
        // Where    Bx = 1/3 * width of canvas, and By = 1/2 * height of canvas
        //          r = width of canvas / 2 = height of canvas / 2
        //          a = `angle`
        //          Sr = `scatter_dist`
        //          Sa = `scatter_angle`
        // This particular results in all of the points being located on or slightly outside the circle defined by the equation
        // [x, y] = [rcos(a), rsin(a)]
        final Vector base = new Vector (two_space,
                                        kCanvasWidth / 3 + kCanvasWidth / 4 * cos (angle) + scatter_dist * cos (angle + cos (scatter_angle) * cos (scatter_angle)),
                                        kCanvasHeight / 2 + kCanvasHeight / 4 * sin (angle) + scatter_dist * sin (angle + sin (scatter_angle) * sin (scatter_angle)));
        // The default color of the particle is white, with 3.125% opacity; using a low opacity allows lines to "stack" and is what gives the sense of gradation
        // I picked up the trick from https://sighack.com/post/getting-creative-with-perlin-noise-fields
        color individual  = kDefaultLineColor;
        for (int k = 0; k < color_anchors.size (); k++) {
            if (color_anchors.get (k).distance (base) < kColoringFalloffDistance * Pareto.pareto (1, 0.1)) {
                // A dark pink-ish red, again with 3.125% opacity
                individual = kAccentLineColor;
            }
        }
        // Generate the particles
        for (int j = 0; j < num_children; j++) {
            particles[particle_index] = base.add (new Vector (two_space, 3 * randomGaussian (), 3 * randomGaussian ()));
            colors[particle_index]    = individual;
            angles[particle_index]    = angle;
            particle_index++;
        }
        total_particles += num_children;
    }

    if(kDiagnostics) {
        println ("[INCANDESCE] ==> GENERATED " + total_particles + " PARTICLES");
        println ("== FINISHED PARTICLE GENERATION ==");
        println ("");
        println ("PARTICLE GENERATION STATS:");
        println ("    RAN/X: " + ran + " " + rax);
        println ("    RSAN/X: " + rsan + " " + rsax);
        println ("    RSDN/X: " + rsdn + " " + rsdx);
        println ("    RCN/X: " + rcn + " " + rcx);

        println ("TTL: " + total_particles);
    }

    if(kDiagnostics)
        println ("== RUNNING PARTICLE SYSTEM ==");

    // Diagnostics
    float rnn = 0.5, rnx = 0.5;

    background(0x0A122A);
    // Slightly stronger stroke weight so
    strokeWeight (3);

    // Run the particle system
    for (int i = 0; i < kParticleSystemSteps; i++) {
        if (kDiagnostics && i % 10 == 9) {
            println ("[INCANDESCE] ==> PSYS ITERATION " + i);
        }
        for (int j = 0; j < total_particles; j++) {
            float noise = noise ((float)particles[j].tuple[0] * kPerlinNoiseSmoothingFactor, (float)particles[j].tuple[1] * kPerlinNoiseSmoothingFactor);

            // Diagnostics
            if (noise < rnn) rnn = noise;
            if (noise > rnx) rnx = noise;

             // The polar transform that gives the image its radiant disposition.
             // I noticed that the Perlin noise field, at least over the dimensions x \in [0, 10] and y \in [0, 10] seems to average at around 2pi * kPerlinFallaway
             // So we subtract that to rotate the average value to the horizontal to the right of the center, and then we add `angles[j]` to rotate the average to point outwards.
            float angle = noise * TWO_PI + angles[j] - (TWO_PI * kPerlinFallaway);
            // Move the particle according to the value of the velocity field at the point
            particles[j].add_to_self (new Vector (two_space, cos (angle), sin (angle)));
            // Set the color
            stroke (colors[j]);
            // Paint
            point ((float)particles[j].tuple[0], (float)particles[j].tuple[1]);
        }
    }

    if(kDiagnostics)

    if(kDiagnostics) {
        println ("RNN/X: " + rnn + " " + rnx);
        println ("== FINISHED ==");
    }
}

/**
 * Regenerate image
 */
void mouseClicked() {
  redraw();
}
