/**
 * Implementation of Ebeida et al.'s algorithm for n-dimensional Poisson disk sampling
 * in a hyperrectangular space.
 */
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
        for (int i = 0; i < this.base_grid_offsets.length; i++) {
            this.base_grid_offsets[i] = offset_Integers[i].intValue();
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
        }
    }

    public void iterate() {
        final double effective_cell_size = this.base_cell_size / pow(2, this.depth);

        final int dart_throws = (int) (this.param_a * this.actives.size());
        for (int i = 0; i < dart_throws; i++) {
            int active = (int)Uniform.uniform(0, this.actives.size() - 1);
            Vector dart = new Vector(this.space);
            for (int j = 0; j < this.space.dimension; j++) {
                dart.tuple[j] = Uniform.uniform_exclude_both(0, effective_cell_size);
            }
            dart.add_to_self(this.actives.get(active).mul(effective_cell_size));
            final int base_index = this.base_grid_offset_for_position(dart, this.depth);
            if (this.base_grid[base_index] != null) {
                this.actives.remove(active);
            } else {
                if (!this.is_covered(dart, this.depth)) {
                    this.base_grid[base_index] = new Hypersphere(this.space, dart, this.radius);
                    this.actives.remove(active);
                }
            }
        }

        this.subdivide_actives();

        this.depth++;
    }

    public boolean computation_finished() {
        return this.actives.size() == 0;
    }

    public void subdivide_actives() {
        for (int i = 0; i < this.actives.size(); i++) {
            if (this.is_covered(i, this.depth)) {
                this.actives.remove(i);
            }
        }

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
                    this.actives.add(v.cast_self_to_int());
                }
            }
        }
    }

    public Hypercube[] subdivide_active(int active, int depth) {
        Hypercube[] cells = new Hypercube[(int) pow(2, this.space.dimension)];
        final Hypercube parent = this.get_active_cell(this.actives.get(active), depth);
        final double new_size = parent.size / 2;

        final Vector transform_space = new Vector(this.space).fill(2);
        final Vector transform_space_offset = new Vector(this.space).fill(-0.5);

        for (int i = 0; i < cells.length; i++) {
            Vector v = new Vector(this.space).from_linear_index(i, transform_space).add_to_self(transform_space_offset).mul_self_by(new_size).add_to_self(parent.center);
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
        for (int i = 0; i < this.base_grid_offsets.length; i++) {
            final int base_offset = base_index + this.base_grid_offsets[i];
            if (base_offset >= 0 && base_offset < this.base_grid.length) {
                if (this.base_grid[base_offset] != null) {
                    if (this.base_grid[base_offset].contains(coords)) {
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
                        return true;
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
        return tmp.linear_index_transform(this.base_grid_size);
    }

    public int base_grid_offset_for_active_cell(Vector coords, int depth) {
        Vector v = coords.div(pow(2, depth)).cast_self_to_int();
        return v.linear_index_transform(this.base_grid_size);
    }
}
