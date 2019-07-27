package brink.graph;

import java.util.ArrayList;

/**
 * Represents a subgraph of an undirected graph.
 */
public class UndirectedSubgraph<T> {

    /**
     * The graph this subgraph is a part of.
     */
    public UndirectedGraph<T> graph;

    /**
     * The indices of the vertices of the parent graph that are included in this
     * subgraph. Because Java doesn't allow primitive-specialized generics, we
     * must use Integer instead of int.
     */
    public ArrayList<Integer> vertices;

    /**
     * Create a subgraph from a graph and a list of vertex indices.
     */
    public UndirectedSubgraph(UndirectedGraph<T> graph, ArrayList<Integer> vertices) {
        this.graph = graph;
        // We want a shallow copy of the list of vertices; a deep copy would be preferable, but I can't find an elegant way to do it
        // In any case, better than using a reference
        this.vertices = (ArrayList<Integer>)vertices.clone();
    }

    /**
     * Check if the subgraph is independent.
     * @return true if independent, false otherwise
     */
    public boolean is_independent() {
        for(int i = 0;i < this.vertices.size();i++) {
            for(int j = i;j < this.vertices.size();j++) {
                if(this.graph.connected(this.vertices.get(i),this.vertices.get(j))) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Check if the subgraph would be independent if the vertex with the index `rhs`
     * were to be added to it.
     *
     * @return true if independent, false otherwise
     */
    public boolean is_independent_of(int rhs) {
        for(int i = 0;i < this.vertices.size();i++) {
            if(this.graph.connected(this.vertices.get(i), rhs)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Return a copy of the subgraph that includes another vertex `rhs`.
     *
     * If the specified vertex is already present in the subgraph, then a copy of
     * the subgraph will be returned.
     */
    public UndirectedSubgraph<T> extend(Integer rhs) {
        UndirectedSubgraph<T> subgraph = new UndirectedSubgraph<T>(this.graph, this.vertices);
        subgraph.extend_self(rhs);
        return subgraph;
    }

    /**
     * Returns a new subgraph that does not include the vertex `rhs`.
     *
     * If the specified vertex is not present in the subgraph, then a copy of the
     * subgraph will be returned.
     */
    public UndirectedSubgraph<T> contract(Integer rhs) {
        UndirectedSubgraph<T> subgraph = new UndirectedSubgraph<T>(this.graph, this.vertices);
        subgraph.contract_self(rhs);
        return subgraph;
    }

    /**
     * Add the vertex with the index `rhs` to the subgraph, if not already present.
     */
    public void extend_self(Integer rhs) {
        if(!this.vertices.contains(rhs)) {
            this.vertices.add(rhs);
        }
    }

    /**
     * Remove the vertex with the index `rhs` from the subgraph, if present.
     */
    public void contract_self(Integer rhs) {
        if(this.vertices.contains(rhs)) {
            this.vertices.remove(rhs);
        }
    }

    /**
     * Works like {@see extend()}, however, will return null if the resulting
     * subgraph is not independent.
     */
    public UndirectedSubgraph<T> extend_if_independent(Integer rhs) {
        // Use is_independent_of instead of is_independent 'cause it's way cheaper.
        if(!this.is_independent_of(rhs)) {
            return null;
        }
        return this.extend(rhs);
    }

    /**
     * Works like {@see extend_self()}, however, will return false if the resulting
     * subgraph is not independent. Will return true otherwise.
     */
    public boolean extend_self_if_independent(Integer rhs) {
        // Use is_independent_of instead of is_independent 'cause it's way cheaper.
        if(!this.is_independent_of(rhs)) {
            return false;
        }
        this.extend_self(rhs);
        return true;
    }

    /**
     * Get the number of vertices in the subgraph.
     */
    public int size() {
        return this.vertices.size();
    }

    /**
     * @internal Get the index of the vertex of th graph that is represented by the `index`th
     * item of the vertices list of this subgraph.
     */
    public int get_vertex(int index) {
        return (int)this.vertices.get(index);
    }
}
