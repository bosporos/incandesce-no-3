import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import brink.*; 
import brink.random.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class PDS extends PApplet {




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

 /**
  * Processing does not supply the natural number as a constant.
  */
 public static final float E = 2.71828182846f;

 /**
  * Canvas width and height, in pixels
  */
 public static final float kCanvasWidth = 1000;
 public static final float kCanvasHeight = 1000;

 /**
  * Fallaway parameter for noiseDetail()
  */
 public static final float kPerlinFallaway = 0.5f;

 /**
  */
 public static final float kPerlinOctaves = 4;

 /**
  * Multiplier for the parameters passed to the noise() function, to ensure the appropriate degree of smoothness
  * in the Perlin noise field.
  */
 public static final float kPerlinNoiseSmoothingFactor = 0.00325f;

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
 public static final boolean kDiagnostics = false;

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

public void settings() {
    size((int) kCanvasWidth, (int) kCanvasHeight, P2D);
}

public void setup() {
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

public void draw() {
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
        // Original ordering: 3, 2, 1, round(16 * 0 / PI)

        // The position of the cross section in the torus, defined as the angle that it occurs at relative the origin point
        float angle         = (float)points.get (i).tuple[0];
        float scatter_angle = (float)points.get (i).tuple[1];
        float scatter_dist  = (float)points.get (i).tuple[2];
        // Calculates the number of particles to be generated around the point
        // Biased towards higher numbers
        int num_children    = (int) round(16 - points.get(i).tuple[3] * log((float) pow((float) points.get(i).tuple[0], E / 2)));

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
                                        333 + 250 * cos (angle) + scatter_dist * cos (angle + cos (scatter_angle) * cos (scatter_angle)),
                                        500 + 250 * sin (angle) + scatter_dist * sin (angle + sin (scatter_angle) * sin (scatter_angle)));
        // The default color of the particle is white, with 3.125% opacity; using a low opacity allows lines to "stack" and is what gives the sense of gradation
        // I picked up the trick from https://sighack.com/post/getting-creative-with-perlin-noise-fields
        int individual  = kDefaultLineColor;
        for (int k = 0; k < color_anchors.size (); k++) {
            if (color_anchors.get (k).distance (base) < kColoringFalloffDistance * Pareto.pareto (1, 0.1f)) {
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
    float rnn = 0.5f, rnx = 0.5f;

    background(0x0A122A);
    // Slightly stronger stroke weight so
    strokeWeight (3);

    // Run the particle system
    for (int i = 0; i < kParticleSystemSteps; i++) {
        if (kDiagnostics && i % 10 == 0) {
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

    if(kDiagnostics) {
        println ("RNN/X: " + rnn + " " + rnx);
        println ("== FINISHED ==");
    }
}

/**
 * Regenerate image
 */
public void mouseClicked() {
  redraw();
}
public class Space {
  public final int dimension;
  
  public Space(int dimension) {
    this.dimension = dimension;
  }
}

/**
 * Vector in N-dimensional Euclidean space
 */
public class Vector {
  Space space;
  public double[] tuple;
  
  public Vector(Space space) {
    this.space = space;
    this.tuple = new double[space.dimension];
    for(int i = 0;i < space.dimension;i++) {
      this.tuple[i] = 0;
    }
  }
  
  public Vector(Space space, double... tuple) {
    this.space = space;
    this.tuple = new double[space.dimension];
    this.set(0, tuple);
  }
  
  public Vector from_linear_index(int index, Vector size) {
    Vector temp = size.clone();
    temp.tuple[this.space.dimension - 1] = 1;
    for(int i = this.space.dimension - 2;i >= 0;i--) {
      temp.tuple[i] = (int)temp.tuple[i] * (int)temp.tuple[i + 1];
    }
    for(int i = 0;i < this.space.dimension - 1 && index > 0;i++) {
      final int diff = (int) (index / temp.tuple[i]);
      this.tuple[i] = diff;
      index %= temp.tuple[i];
    }
    this.tuple[this.space.dimension - 1] = index;
    return this;
  }
  
    
  public int linear_index_transform(Vector size) {
    Vector temp = size.clone();
    temp.tuple[temp.space.dimension - 1] = 1;
    for(int i = temp.space.dimension - 2;i >= 0;i--) {
      temp.tuple[i] = (int)temp.tuple[i] * (int)temp.tuple[i + 1];
    }
    return (int)this.dot(temp);
  }
  
  public Vector cast_self_to_int() {
    for(int i = 0;i < this.space.dimension;i++) {
      this.tuple[i] = (int) this.tuple[i];
    }
    return this;
  }
  
  public Vector cast_to_int() {
    return this.clone().cast_self_to_int();
  }
  
  public Vector round_self() {
    for(int i = 0;i < this.space.dimension;i++) {
      this.tuple[i] = Math.round(this.tuple[i]);
    }
    return this;
  }
  
  public Vector round() {
    return this.clone().round_self();
  }
  
  public Vector fill(double value) {
    for(int i = 0;i < this.space.dimension;i++) {
      this.tuple[i] = value;
    }
    return this;
  }
  
  public Vector set(double... tuple) {
    return this.set(0, tuple);
  }
  
  public Vector set(int offset, double... tuple) {
    for(int i = 0;offset < this.space.dimension;) {
      this.tuple[offset++] = tuple[i++];
    }
    return this;
  }
  
  public double get(int offset) {
    return this.tuple[offset];
  }
  
  public double[] get(int offset, int number) {
    double[] slice = new double[number - (this.space.dimension - offset - number)];
    for(int i = 0;i < number && offset < this.space.dimension;) {
      slice[i++] = this.tuple[offset++]; 
    }
    return slice;
  }
  
  public Vector add_to_self(Vector rhs) {
    for(int i = 0;i < this.space.dimension;i++) {
      this.tuple[i] += rhs.tuple[i];
    }
    return this;
  }
  
  public Vector add(Vector rhs) {
    return this.clone().add_to_self(rhs);
  }
  
  public Vector sub_from_self(Vector rhs) {
    for(int i = 0;i < this.space.dimension;i++) {
      this.tuple[i] -= rhs.tuple[i];
    }
    return this;
  }
  
  public Vector sub(Vector rhs) {
    return this.clone().sub_from_self(rhs);
  }
  
  public Vector mul_self_by(double scalar) {
    for(int i = 0;i < this.space.dimension;i++) {
      this.tuple[i] *= scalar;
    }
    return this;
  }
  
  public Vector mul(double scalar) {
    return this.clone().mul_self_by(scalar);
  }
  
  public Vector div_self_by(double scalar) {
    for(int i = 0;i < this.space.dimension;i++) {
      this.tuple[i] /= scalar;
    }
    return this;
  }
  
  public Vector div(double scalar) {
    return this.clone().div_self_by(scalar);
  }
  
  public double length() {
    double sum = 0;
    for(int i = 0;i < this.space.dimension;i++) {
      sum += this.tuple[i] * this.tuple[i];
    }
    return Math.sqrt(sum);
  }
  
  public double distance(Vector rhs) {
    return this.sub(rhs).length();
  }
  
  public double maximum() {
    double highest = 0;
    for(int i = 0;i < this.space.dimension;i++) {
      if(this.tuple[i] > highest)
        highest = this.tuple[i];
    }
    return highest;
  }
  
  public double minimum() {
    double lowest = this.maximum();
    for(int i = 0;i < this.space.dimension;i++) {
      if(this.tuple[i] < lowest)
        lowest = this.tuple[i];
    }
    return lowest;
  }
  
  public Vector normalize_self() {
    this.div_self_by(this.maximum());
    return this;
  }
  
  public Vector normalize() {
    return this.clone().normalize_self();
  }
  
  public Vector abs_self() {
    for(int i = 0;i < this.space.dimension;i++) {
      if(this.tuple[i] < 0) {
        this.tuple[i] = -this.tuple[i];
      }
    }
    return this;
  }
  
  public Vector abs() {
    return this.clone().abs_self();
  }
  
  public double product() {
    double product = 1;
    for(int i = 0;i < this.space.dimension;i++) {
      product *= this.tuple[i];
    }
    
    return product;
  }
  
  public double sum() {
    double sum = 0;
    for(int i = 0;i < this.space.dimension;i++) {
      sum += this.tuple[i];
    }
    
    return sum;
  }
  
  public double dot(Vector rhs) {
    double sum = 0;
    for(int i = 0;i < this.space.dimension;i++) {
      sum += this.tuple[i] * rhs.tuple[i];
    }
    return sum;
  }
  
  public Vector cross(Vector rhs) {
    if(this.space.dimension != 3 && this.space.dimension != 7) {
      // Pretty sure cross product isn't defined outside of 3 or 7 dimensions
      return new Vector(this.space);
    } else if(this.space.dimension == 3) {
      // https://en.wikipedia.org/wiki/Cross_product
      return new Vector(this.space,
        this.tuple[1]*rhs.tuple[2] - this.tuple[2]*rhs.tuple[1],
        this.tuple[2]*rhs.tuple[0] - this.tuple[0]*rhs.tuple[2],
        this.tuple[0]*rhs.tuple[1] - this.tuple[1]*rhs.tuple[0]
      );
    } else {
      // https://en.wikipedia.org/wiki/Seven-dimensional_cross_product
      return new Vector(this.space,
        this.tuple[1]*rhs.tuple[2] - this.tuple[2]*rhs.tuple[1] + this.tuple[3]*rhs.tuple[4] - this.tuple[4]*rhs.tuple[3] + this.tuple[5]*rhs.tuple[6] - this.tuple[6]*rhs.tuple[5],
        - this.tuple[0]*rhs.tuple[2] + this.tuple[2]*rhs.tuple[0] + this.tuple[3]*rhs.tuple[5] + this.tuple[4]*rhs.tuple[6] - this.tuple[5]*rhs.tuple[3] - this.tuple[6]*rhs.tuple[4],
        this.tuple[0]*rhs.tuple[1] - this.tuple[1]*rhs.tuple[0] + this.tuple[3]*rhs.tuple[6] - this.tuple[4]*rhs.tuple[5] + this.tuple[5]*rhs.tuple[4] - this.tuple[6]*rhs.tuple[3],
        - this.tuple[0]*rhs.tuple[4] - this.tuple[1]*rhs.tuple[5] - this.tuple[2]*rhs.tuple[6] + this.tuple[4]*rhs.tuple[0] + this.tuple[5]*rhs.tuple[1] + this.tuple[6]*rhs.tuple[2],
        this.tuple[0]*rhs.tuple[3] - this.tuple[1]*rhs.tuple[7] + this.tuple[2]*rhs.tuple[5] - this.tuple[3]*rhs.tuple[0] - this.tuple[5]*rhs.tuple[2] + this.tuple[6]*rhs.tuple[1],
        this.tuple[0]*rhs.tuple[6] + this.tuple[1]*rhs.tuple[3] - this.tuple[2]*rhs.tuple[4] - this.tuple[6]*rhs.tuple[0] - this.tuple[3]*rhs.tuple[1] + this.tuple[4]*rhs.tuple[2],
        - this.tuple[0]*rhs.tuple[5] + this.tuple[1]*rhs.tuple[4] + this.tuple[2]*rhs.tuple[3] + this.tuple[5]*rhs.tuple[0] - this.tuple[4]*rhs.tuple[1] - this.tuple[3]*this.tuple[2]
      );
    }
  }

  public boolean equals(Vector rhs) {
    if(this.space.dimension != rhs.space.dimension) {
      return false;
    } else {
      for(int i = 0;i < this.space.dimension;i++) {
        if(this.tuple[i] != rhs.tuple[i]) {
          return false;
        }
      }
      return true;
    }
  }
  
  public boolean less_than(Vector rhs) {
    final int effective_dimension = (int)min(this.space.dimension, rhs.space.dimension);
    for(int i = 0;i < effective_dimension;i++) {
      if(rhs.tuple[i] < this.tuple[i]) {
        return false;
      }
    }
    return true;
  }
  
  public boolean greater_than(Vector rhs) {
    final int effective_dimension = (int)min(this.space.dimension, rhs.space.dimension);
    for(int i = 0;i < effective_dimension;i++) {
      if(rhs.tuple[i] > this.tuple[i]) {
        return false;
      }
    }
    return true;
  }
  
  public Vector clone() {
    Vector copy = new Vector(this.space);
    for(int i = 0;i < this.space.dimension;i++) {
      copy.tuple[i] = this.tuple[i];
    }
    return copy;
  }
  
  /**
   * Debugging purposes.
   */
  public String toString() {
    String representation = "[";
    for(int i = 0;i < this.space.dimension - 1;i++) {
      representation += this.tuple[i] + " ";
    }
    return representation + this.tuple[this.space.dimension - 1] + "]";
  }
}

abstract public class Shape {
  public Space space;
  
  public Shape(Space space) {
    this.space = space;
  }
  
  abstract public boolean contains(Vector point);
}

public class Hyperrectangle extends Shape {
  public Vector center;
  public Vector size;
  
  public Hyperrectangle(Space space, Vector center, Vector size) {
    super(space);
    this.center = center;
    this.size = size;
  }
  
  public boolean contains(Vector point) {
    final Vector offset = this.center.sub(point).abs_self();
    for(int i = 0;i < this.space.dimension;i++) {
      if(offset.tuple[i] > (this.size.tuple[i] / 2)) {
        return false;
      }
    }
    
    return true;
  }
}

public class Hypercube extends Shape {
  public Vector center;
  public double size;
  
  public Hypercube(Space space, Vector center, double size) {
    super(space);
    this.center = center;
    this.size = size;
  }
    
  public double get_longest_diagonal() {
    return Math.sqrt(this.size);
  }
  
  public Vector get_furthest_vertex_from(Vector point) {
    Vector furthest_vertex = new Vector(this.space);
    for(int i = 0;i < this.space.dimension;i++) {
      if(this.center.tuple[i] < point.tuple[i]) {
        furthest_vertex.tuple[i] = this.center.tuple[i] - this.size / 2;
      } else {
        furthest_vertex.tuple[i] = this.center.tuple[i] + this.size / 2;
      }
    }
    return furthest_vertex;
  }
  
  public Hyperrectangle as_hyperrectangle() {
    Vector size = new Vector(this.space);
    for(int i = 0;i < this.space.dimension;i++) {
      size.tuple[i] = this.size;
    }
    return new Hyperrectangle(this.space, center, size);
  }
  
  public boolean contains(Vector point) {
    final Vector offset = this.center.sub(point).abs_self();
    
    for(int i = 0;i < this.space.dimension;i++) {
      if(offset.tuple[i] > (this.size / 2)) {
        return false;
      }
    }
    
    return true;
  }
  
  public String toString() {
    return "[H^CUBE " + this.center + " SIZE " + this.size + "]";
  }
}

public class Hypersphere extends Shape {
  public Vector center;
  public double radius;
  
  public Hypersphere(Space space, Vector center, double radius) {
    super(space);
    this.center = center;
    this.radius = radius;
  }
  
  public boolean contains(Vector point) {
    return this.center.distance(point) < this.radius;
  }
  
  public String toString() {
    return "[H^SPHERE " + this.center + " RAD " + this.radius + "]";
  }
}
public class PoissonDiskSample {

  public Hyperrectangle area;
  Space space;

  public double radius;
  public double base_cell_size;

  public ArrayList<Vector> actives;
  public Hypersphere base_grid[];
  public Vector base_grid_size;

  public int base_grid_offsets[];

  public float param_a;

  public int depth;

  public PoissonDiskSample(Hyperrectangle area, Space space, double radius, float A) {
    this.space = space;
    this.area = area;

    this.radius = radius;
    this.base_cell_size = radius / Math.sqrt((double)space.dimension);

    this.param_a = A;

    this.base_grid_size = area.size.div(this.base_cell_size).cast_self_to_int();
    println("Base grid size (Vec): " + this.base_grid_size);
    println("Base grid size (Prod): " + this.base_grid_size.product());

    this.base_grid = new Hypersphere[(int)this.base_grid_size.product()];

    this.depth = -1;

    // Precompute the base grid offsets of all the base grid cells that could potentially have samples that would render the base grid cell non-disk-free

    Vector search_space = new Vector(space).fill(5);
    int search_space_length = (int)pow(5, space.dimension);
    Vector center = search_space.sub(new Vector(space).fill(3));
    ArrayList<Integer> offsets = new ArrayList<Integer>();
    int insert_orthogonally_adjacent = 0;
    int insert_diagonally_adjacent = 0;
    int insert_nonimmediately_adjacent = 0;
    for (int i = 0; i < search_space_length; i++) {
      Vector position = new Vector(space).from_linear_index(i, search_space);
      position.sub_from_self(center);
      int offset = position.linear_index_transform(this.base_grid_size);
      if (position.abs().sum() == 2) {
        if (position.abs().maximum() == 1) {
          offsets.add(insert_orthogonally_adjacent + insert_diagonally_adjacent++, offset);
        } else {
          offsets.add(insert_orthogonally_adjacent + insert_diagonally_adjacent + insert_nonimmediately_adjacent++, offset);
        }
      } else if (position.abs().sum() == 1) {
        offsets.add(insert_orthogonally_adjacent++, offset);
      }
    }
    Integer offset_Integers[] = (Integer[]) offsets.toArray(new Integer[offsets.size()]);
    this.base_grid_offsets = new int[offset_Integers.length];
    //println("Base grid offsets:");
    for (int i = 0; i < this.base_grid_offsets.length; i++) {
      this.base_grid_offsets[i] = offset_Integers[i].intValue();
      //println("Offset " + i + " -> " + this.base_grid_offsets[i]);
    }
  }

  public void begin_computation() {
    this.depth = 0;

    for (int i = 0; i < this.base_grid.length; i++) {
      this.base_grid[i] = null;
    }
    this.actives = new ArrayList<Vector>((int)this.base_grid_size.product());
    for (int i = 0; i < this.base_grid.length; i++) {
      this.actives.add(i, new Vector(this.space).from_linear_index(i, this.base_grid_size));
      //println("Active " + i + " -> (ACTIVE " + new Vector(this.space).from_linear_index(i, this.base_grid_size) + ")");
    }
  }

  public void iterate() {
    final double effective_cell_size = this.base_cell_size / pow(2, this.depth);

    //println("Throwing darts... ");
    final int dart_throws = (int) (this.param_a * this.actives.size());
    for (int i = 0; i < dart_throws; i++) {
      int active = (int)Uniform.uniform(0, this.actives.size() - 1);
      Vector dart = new Vector(this.space);
      for (int j = 0; j < this.space.dimension; j++) {
        dart.tuple[j] = Uniform.uniform_exclude_both(0, effective_cell_size);
      }
      //print("Dart: (BOX " + dart + " ON " + this.actives.get(active).mul(effective_cell_size) + ")");
      dart.add_to_self(this.actives.get(active).mul(effective_cell_size));
      //print(" -> (DART " + dart + ")");
      final int base_index = this.base_grid_offset_for_position(dart, this.depth);
      //println(" -> " + base_index);
      //print("    [DART] -> ");
      if (this.base_grid[base_index] != null) {
        //println("(DEAD SAME)");
        this.actives.remove(active);
      } else {
        if (!this.is_covered(dart, this.depth)) {
          //println("(PLACE)");
          this.base_grid[base_index] = new Hypersphere(this.space, dart, this.radius);
          this.actives.remove(active);
        } else {
          //println("(DEAD COVER)");
        }
      }
    }

    //println("Subdividing... ");
    this.subdivide_actives();

    this.depth++;
  }

  public boolean computation_finished() {
    return this.actives.size() == 0;
  }

  public void subdivide_actives() {
    for (int i = 0; i < this.actives.size(); i++) {
      if (this.is_covered(i, this.depth)) {
        //println(" REM (ACTIVE " + this.actives.get(i) + ")");
        this.actives.remove(i);
      }
    }
    
    //println("REMOVE END");

    final double child_size = this.base_cell_size / pow(2, this.depth + 1);
    final Vector center_corrector = new Vector(this.space).fill(child_size/2);

    int prior_size = this.actives.size();
    for (int i = 0; i < prior_size; i++) {
      final int base_index = this.base_grid_offset_for_active_cell(this.actives.get(0), depth);
      Hypercube[] children = this.subdivide_active(0, this.depth);
      this.actives.remove(0);
      for (Hypercube child : children) {
        if (!this.is_covered(child, base_index)) {
          // Use sub_from_self because the `child` will be deleted after this, and sub() would just be a waste of an allocation
          Vector v = child.center.sub_from_self(center_corrector).div_self_by(child_size).round_self();
          //print(v);
          this.actives.add(v.cast_self_to_int());
          //println(" CHILD -> " + child);
        }
      }
    }
  }

  public Hypercube[] subdivide_active(int active, int depth) {
    Hypercube[] cells = new Hypercube[(int) pow(2, this.space.dimension)];
    final Hypercube parent = this.get_active_cell(this.actives.get(active), depth);
    final double new_size = parent.size / 2;

    final Vector transform_space = new Vector(this.space).fill(2);
    final Vector transform_space_offset = new Vector(this.space).fill(-0.5f);

    for (int i = 0; i < cells.length; i++) {
      Vector v = new Vector(this.space).from_linear_index(i, transform_space).add_to_self(transform_space_offset).mul_self_by(new_size).add_to_self(parent.center);
      //println("SUBDIV " + parent + " (" + i + ") -> " + v);
      cells[i] = new Hypercube(this.space, v, new_size);
    }

    return cells;
  }

  public ArrayList<Vector> get_samples() {
    ArrayList<Vector> points = new ArrayList<Vector>();
    for (int i = 0; i < this.base_grid.length; i++) {
      if (this.base_grid[i] != null) {
        points.add(this.base_grid[i].center);
      }
    }
    return points;
  }

  public boolean is_covered(Hypersphere sphere, Vector cell, int depth) {
    return sphere.contains(this.get_active_cell(cell, depth).get_furthest_vertex_from(sphere.center));
  }

  public boolean is_covered(Vector coords, int depth) {
    final int base_index = this.base_grid_offset_for_position(coords.cast_to_int(), depth);
    //print("(COVER " + coords + " -> BASE " + base_index + ") ");
    for (int i = 0; i < this.base_grid_offsets.length; i++) {
      final int base_offset = base_index + this.base_grid_offsets[i];
      if (base_offset >= 0 && base_offset < this.base_grid.length) {
        if (this.base_grid[base_offset] != null) {
          if (this.base_grid[base_offset].contains(coords)) {
            //print(this.base_grid[base_offset]);
            return true;
          }
        }
      }
    }
    return false;
  }

  public boolean is_covered(int active, int depth) {
    return this.is_covered(this.get_active_cell(this.actives.get(active), depth), this.base_grid_offset_for_active_cell(this.actives.get(active), depth));
  }

  public boolean is_covered(Hypercube cell, int base_index) {
    for (int i = 0; i < this.base_grid_offsets.length; i++) {
      final int base_offset = base_index + this.base_grid_offsets[i];
      if (base_offset >= 0 && base_offset < this.base_grid.length) {
        if (this.base_grid[base_offset] != null) {
          if (this.base_grid[base_offset].contains(cell.get_furthest_vertex_from(this.base_grid[base_offset].center))) {
            //println("COVERED BASE " + base_index + " -> " + this.base_grid[base_offset] + " INTERSECT " + cell + " @ " + cell.get_furthest_vertex_from(this.base_grid[base_offset].center));
            return true;
          } else {
            //println("UNCOVERED BASE " + base_index + " -> " + this.base_grid[base_offset] + " INTERSECT " + cell + " @ " + cell.get_furthest_vertex_from(this.base_grid[base_offset].center));
          }
        }
      }
    }
    return false;
  }

  public Hypercube get_active_cell(Vector coords, int depth) {
    final double unit = this.base_cell_size / pow(2, depth);

    Vector size = new Vector(this.space).fill(unit);
    Vector center = coords.mul(unit).add_to_self(size.div(2));

    return new Hypercube(this.space, center, unit);
  }

  public int base_grid_offset_for_position(Vector coords, int depth) {
    final Vector tmp = coords.div(this.base_cell_size).cast_self_to_int();
    //print(" (BASE POS " + tmp + ") ");
    return tmp.linear_index_transform(this.base_grid_size);
  }

  public int base_grid_offset_for_active_cell(Vector coords, int depth) {
    Vector v = coords.div(pow(2, depth)).cast_self_to_int();
    //print(" {BGOFAC: " + v + "} ");
    return v.linear_index_transform(this.base_grid_size);
  }
}
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "PDS" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
