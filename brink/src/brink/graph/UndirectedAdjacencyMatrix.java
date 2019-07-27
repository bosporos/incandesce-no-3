package brink.graph;

import java.util.ArrayList;

/**
 * An undirected adjacency matrix with a fixed number of vertices
 *
 * Technical note: since undirected graphs have symmetric adjacency matrices, we can optimize!
 *
 * For n vertices, total number of necessary adjacencies will be the nth triangular number
 * Thus, we can use the formula a(n) = floor((n+1)^3/(n+2))/2 for the nth triangular number,
 * courtesy of Gary Detlefs (representation shown here adjusted for starting offset of 1), [https://oeis.org/A000217]
 * for the total number of vertices
 * Then, individual adjacency values can be accessed with an extension of this formula;
 * given two vertices b, c where b <= c
 *          offset in array = a(b - 1) + c
 */
public class UndirectedAdjacencyMatrix<T> implements UndirectedGraph<T> {

    /**
     * Vertex values
     */
    public ArrayList<T> vertices;

    /**
     * Actual adjacency matrix.
     * See the technical note in {@see UndirectedAdjacencyMatrix<T>} for more information
     * about the indexing scheme.
     */
    public boolean[] adjacencies;

    /**
     * The number of vertices in the graph.
     */
    public final int dimension;

    public UndirectedAdjacencyMatrix(int dimension) {
        this.vertices = new ArrayList<>(dimension);
        for(int i = 0;i < dimension;i++) {
            this.vertices.add(null);
        }
        this.dimension = dimension;
        dimension++;
        this.adjacencies = new boolean[
            // floor(n^3)/(n+1)/2
            ((dimension*dimension*dimension)/(dimension+1))/2
        ];
    }

    @Override
    public void connect(int lhs, int rhs) {
        final int p = (int)Math.min(lhs, rhs);
        final int s = (int)Math.max(lhs, rhs) + 1;
        this.adjacencies[s*s*s/(s+1)/2+p] = true;
    }

    @Override
    public void disconnect(int lhs, int rhs) {
        final int p = (int)Math.min(lhs, rhs);
        final int s = (int)Math.max(lhs, rhs) + 1;
        this.adjacencies[s*s*s/(s+1)/2+p] = false;
    }

    @Override
    public boolean connected(int lhs, int rhs) {
        final int p = (int)Math.min(lhs, rhs);
        final int s = (int)Math.max(lhs, rhs) + 1;
        return this.adjacencies[s*s*s/(s+1)/2+p];
    }

    @Override
    public int[] get_neighbours(int vertex) {
        int num_neighbours = 0;
        for(int i = 0;i < this.dimension;i++) {
            if(this.connected(vertex, i)) {
                num_neighbours++;
            }
        }
        int[] neighbours = new int[num_neighbours];
        int i = 0;
        for(int neighbour_index = 0;neighbour_index < num_neighbours;neighbour_index++) {
            if(this.connected(vertex, i++)) {
                neighbours[neighbour_index] = i;
            }
        }
        return neighbours;
    }

    @Override
    public int size() {
        return this.dimension;
    }

    @Override
    public T get_vertex(int vertex) {
        return this.vertices.get(vertex);
    }

    @Override
    public void set_vertex(int vertex, T value) {
        this.vertices.set(vertex, value);
    }
}
