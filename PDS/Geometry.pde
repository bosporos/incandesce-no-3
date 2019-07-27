// Everything should be relatively self-explanatory with the exception of from_linear_index and linear_index_transform,
// so I tried to explain & comment those two methods as best I could.


/**
 * Represents an n-dimensional Euclidean space.
 */

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

    /**
     * Inverse operation of linear_index_transform(scaffold).
     *
     * Given P being an N-cube of dimensions D = <d[1], d[2], ..., d[N]>, composed of q = product(D)
     * N-cubes of dimensions <1, 1, ..., 1> C[0], C[1], ... C[q-1] indexed respectively as 0..q-1 let
     * from_linear_index(m, D) give the coordinates V = <v[1], v[2], ... v[N]> of C[m] within P.
     */
    public Vector from_linear_index(int index, Vector scaffold) {
        // Because it's passed by reference, we do this so as not to damage it
        Vector temp = scaffold.clone();
        // Essentially
        // let t[1..n]
        //  t[i - 1] = d[i] * d[i + 1] * ... * d[N] where i in 1..n
        //  t[N] = 1
        temp.tuple[this.space.dimension - 1] = 1;
        for(int i = this.space.dimension - 2;i >= 0;i--) {
            temp.tuple[i] = (int)temp.tuple[i] * (int)temp.tuple[i + 1];
        }
        // Then let
        //  v[i] = floor(m[i] / t[i]) for i in 1..n-1
        // where
        //  m[i + 1] = m[i] mod t[i] for i in 1..n-1
        //  m[1] = m
        for(int i = 0;i < this.space.dimension - 1 && index > 0;i++) {
            final int diff = (int) (index / temp.tuple[i]);
            this.tuple[i] = diff;
            index %= temp.tuple[i];
        }
        // And finally let v[N] = m[N]
        this.tuple[this.space.dimension - 1] = index;
        return this;
    }

    /**
     * Inverse operation of from_linear_index(index, scaffold).
     *
     * Given P being an N-cube of dimensions D = <d[1], d[2], ... d[N]>, composed of q = product(D)
     * N-cubes of dimensions <1, 1, ..., 1> C[0], C[1], ... C[q-1] indexed respectively as 0..q-1 let
     * V.linear_index_transform(D) give the index m of C[m] with coordinates V = <v[1], v[2], ... v[N]>
     * within P.
     */
    public int linear_index_transform(Vector scaffold) {
        // Because it's passed by reference we do this so as not to damage it
        Vector temp = scaffold.clone();
        // Essentially:
        // let t[1..n]
        //  t[i - 1] = d[i] * d[i+1] * ... * d[N]
        //  t[N] = 1
        temp.tuple[temp.space.dimension - 1] = 1;
        for(int i = temp.space.dimension - 2;i >= 0;i--) {
            temp.tuple[i] = (int)temp.tuple[i] * (int)temp.tuple[i + 1];
        }
        // And then do the dot product
        // t[1] * v[1] + t[2] * v[2] + ... + t[N] * v[N]
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
