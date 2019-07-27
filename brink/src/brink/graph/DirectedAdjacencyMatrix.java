package brink.graph;

import java.util.ArrayList;

public class DirectedAdjacencyMatrix<T> implements DirectedGraph<T> {
    public ArrayList<T> vertices;
    public boolean[][] adjacencies;
    final public int dimension;

    public DirectedAdjacencyMatrix(int dimension) {
        this.vertices = new ArrayList<>(dimension);
        for(int i = 0;i < dimension;i++) {
            this.vertices.add(null);
        }
        this.dimension = dimension;
        this.adjacencies = new boolean[this.dimension][this.dimension];
    }

    @Override
    public void connect(int lhs, int rhs) {
        this.adjacencies[lhs][rhs] = true;
    }

    @Override
    public void disconnect(int lhs, int rhs) {
        this.adjacencies[lhs][rhs] = false;
    }

    @Override
    public boolean connected(int lhs, int rhs) {
        return this.adjacencies[lhs][rhs];
    }

    @Override
    public int[] get_neighbours(int vertex) {
        int num_neighbours = 0;
        for(int i = 0;i < this.dimension;i++) {
            if(this.adjacencies[vertex][i]) {
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
