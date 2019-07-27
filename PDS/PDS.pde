import brink.*;
import brink.algorithm.optimization.*;
import brink.draw.*;
import brink.graph.*;
import brink.random.*;

Space three_space = new Space(4);
Space two_space = new Space(2);
PoissonDiskSample toroid;
PoissonDiskSample color_plane;
ArrayList<Vector> points;

void setup() {
  Brink.init(this);

  size(1000, 1000, P2D);

  Vector spatial_dimensions = new Vector(three_space, TWO_PI, TWO_PI, PI * 20, PI);
  Vector center = spatial_dimensions.div(2);

  toroid = new PoissonDiskSample(new Hyperrectangle(three_space, center, spatial_dimensions), three_space, PI, 32);
  color_plane = new PoissonDiskSample(new Hyperrectangle(two_space, new Vector(two_space, 500, 500), new Vector(two_space, 1000, 1000)), two_space, 333, 1);
  
  println("== TESTING TRANSFORM SANITY ==");
  Vector scaffold = new Vector(three_space, 2, 2, 3, 2);
  Vector r = new Vector(three_space).from_linear_index(20, scaffold);
  println("Vector (FLI): " +  r);
  println("Vector (LIT): " + r.linear_index_transform(scaffold));

  noLoop();
  //frameRate(60);
}

float angle = 0;

void draw() {
  noiseSeed((int)random(0, MAX_FLOAT));
  
  background(0x0A122A);
  
  println("== BEGINNING PDS SAMPLE COMPUTATION ==");
  println("[INCANDESCE] ==> BEGINNING COMPUTATION OF HYPERTOROIDAL PDS SAMPLE");
  
  toroid.begin_computation();
  int iter = 0;
  do {
    toroid.iterate();
    println(">>>>>>>>>>>>>>>> ITERATION " + iter + " WITH ACTIVE CELLS: " + toroid.actives.size());
    iter++;
  } while (iter < 6 && !toroid.computation_finished());
  println("            ---> FINISHED");
  points = toroid.get_samples();
  println("\t\tSample size: " + points.size());
  
  println("[INCANDESCE] ==> BEGINNING COMPUTATION OF RECTANGULAR PDS SAMPLE");
  color_plane.begin_computation();
  iter = 0;
  do {
    color_plane.iterate();
    println(">>>>>>>>>>>>>>>> ITERATION " + iter);
    iter++;
  } while (iter < 6 && !color_plane.computation_finished());
  ArrayList<Vector> color_anchors = color_plane.get_samples();
  println("            ---> FINISHED");
  println("== FINISHED PDS SAMPLE COMPUTATION ==");

  Vector particles[] = new Vector[points.size() * 16];
  int colors[] = new int[points.size() * 16];
  float angles[] = new float[points.size() * 16];
  int ttl = 0;
  int idx = 0;
  double ran = PI, rax = PI;
  double rsan = PI, rsax = PI;
  double rsdn = PI, rsdx = PI;
  double rcn = PI, rcx = PI;
  for (int i = 0; i < points.size(); i++) {
    float angle = (float) points.get(i).tuple[3];
    float scatter_angle = (float) points.get(i).tuple[2];
    float scatter_dist = (float) points.get(i).tuple[1];
    int num_children = (int) round((float) (16 * points.get(i).tuple[0] / PI));
    //int num_children = 16;
    if(points.get(i).tuple[3] < ran) ran = points.get(i).tuple[3];
    if(points.get(i).tuple[3] > rax) rax = points.get(i).tuple[3];
    
    if(points.get(i).tuple[2] < rsan) rsan = points.get(i).tuple[2];
    if(points.get(i).tuple[2] > rsax) rsax = points.get(i).tuple[2];
    
    if(points.get(i).tuple[1] < rsdn) rsdn = points.get(i).tuple[1];
    if(points.get(i).tuple[1] > rsdx) rsdx = points.get(i).tuple[1];
    
    if(points.get(i).tuple[0] < rcn) rcn = points.get(i).tuple[0];
    if(points.get(i).tuple[0] > rcx) rcx = points.get(i).tuple[0];
    
    final Vector base = new Vector(two_space,
                                   333 + 250 * cos(angle) + scatter_dist * cos(angle + cos(scatter_angle) * cos(scatter_angle)),
                                   500 + 250 * sin(angle) + scatter_dist * sin(angle + sin(scatter_angle) * sin(scatter_angle)));
    color individual = 0x08FBFAF8;
    for(int k = 0;k < color_anchors.size();k++) {
      if(color_anchors.get(k).distance(base) < 54.36563657 * Pareto.pareto(1, 0.1)) {
        individual = 0x08DB5461;
      }
    }
    for(int j = 0;j < num_children;j++) {
      angles[idx] = angle;
      particles[idx] = base.add(new Vector(two_space, 3 * randomGaussian(), 3 * randomGaussian()));
      colors[idx] = individual;
      idx++;
    }
    ttl += num_children;
  }
  
  println("RAN/X: " + ran + " " + rax);
  println("RSAN/X: " + rsan + " " + rsax);
  println("RSDN/X: " + rsdn + " " + rsdx);
  println("RCN/X: " + rcn + " " + rcx);
  
  println("TTL: " + ttl);

  strokeWeight(3);
  
  println("== RUNNING PARTICLE SYSTEM ==");
  
  float rnn = 0.5, rnx = 0.5;

  for (int i = 0; i < 500; i++) {
    if(i % 10 == 0) {
      println("[INCANDESCE] ==> PSYS ITERATION " + i);
    }
    for(int j = 0;j < ttl;j++) {
      float noise = noise((float)particles[j].tuple[0] * 0.01, (float)particles[j].tuple[1] * 0.01);
      if(noise < rnn) rnn = noise;
      if(noise > rnx) rnx = noise;
      float angle = noise * TWO_PI - PI + angles[j];
      particles[j].add_to_self(new Vector(two_space, cos(angle), sin(angle)));
      stroke(colors[j]);
      point((float)particles[j].tuple[0], (float)particles[j].tuple[1]);
    }
  }
  
  println("RNN/X: " + rnn + " " + rnx);
  println("== FINISHED ==");
}

void mouseClicked() {
  redraw();
}
