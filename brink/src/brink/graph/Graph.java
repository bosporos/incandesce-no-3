package brink.graph;

/**
 * Abstract Graph type containing vertices of type T.
 *
 * For convenience of access, all vertices are refered to by an index;
 * the actual value of the vertex (i.e. a value of type T) is only accessed
 * through {@see get_vertex()} and {@see set_vertex()}.
 */
public interface Graph<T> {

    /**
     * Create an edge between the vertices `lhs` and `rhs`.
     */
    public void connect(int lhs, int rhs);

    /**
     * Destroy the edge between vertices `lhs` and `rhs`, if such exists.
     */
    public void disconnect(int lhs, int rhs);

    /**
     * Determines if the vertices `lhs` and `rhs` are connected by an edge.
     */
    public boolean connected(int lhs, int rhs);

    /**
     * Get an array of the indices of the vertices neighbouring `vertex`.
     */
    public int[] get_neighbours(int vertex);

    /**
     * Get the number of nodes in the graph.
     */
    public int size();

    /**
     * Get the value of the vertex with index `vertex`.
     *
     * Note that because Java treats all objects as references, the "value" returned
     * will in fact be a reference to the one contained within this Graph;
     *
     * <code>
     * class Vertex { public int inner = 0; }
     *
     * // ... snip ...
     *
     * Vertex one = graph.get_vertex(index);
     * Vertex two = graph.get_vertex(index);
     * two.inner = 1;
     * println("One: " + one.inner); // will print "One: 1"
     * two.inner = 3;
     * Vertex three = graph.get_vertex(index);
     * println("Three: " + three.inner) // will print "Three: 3"
     * </code>
     */
    public T get_vertex(int vertex);

    /**
     * Set the value of the vertex with index `vertex` to `value`.
     *
     * Note that Java treats all objects as references, so `value` will be stored
     * as a *reference* within the Graph.
     */
    public void set_vertex(int vertex, T value);
}
