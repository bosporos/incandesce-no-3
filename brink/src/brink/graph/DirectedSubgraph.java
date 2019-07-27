package brink.graph;

import java.util.ArrayList;

public class DirectedSubgraph<T> {
    public DirectedGraph<T> graph;
    public ArrayList<Integer> vertices;

    public DirectedSubgraph(DirectedGraph<T> graph, ArrayList<Integer> vertices) {
        this.graph = graph;
        this.vertices = (ArrayList<Integer>)vertices.clone();
    }

    public boolean is_independent() {
        for(int i = 0;i < this.vertices.size();i++) {
            for(int j = 0;j < this.vertices.size();j++) {
                if(this.graph.connected(this.vertices.get(i),this.vertices.get(j))) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean is_independent_of(int rhs) {
        for(int i = 0;i < this.vertices.size();i++) {
            if(this.graph.connected(this.vertices.get(i), rhs) || this.graph.connected(rhs, this.vertices.get(i))) {
                return false;
            }
        }
        return true;
    }

    public DirectedSubgraph<T> extend(Integer rhs) {
        DirectedSubgraph<T> subgraph = new DirectedSubgraph<T>(this.graph, this.vertices);
        subgraph.extend_self(rhs);
        return subgraph;
    }

    public DirectedSubgraph<T> contract(Integer rhs) {
        DirectedSubgraph<T> subgraph = new DirectedSubgraph<T>(this.graph, this.vertices);
        subgraph.contract_self(rhs);
        return subgraph;
    }

    public void extend_self(Integer rhs) {
        this.vertices.add(rhs);
    }

    public void contract_self(Integer rhs) {
        this.vertices.remove(rhs);
    }

    public DirectedSubgraph<T> extend_if_independent(Integer rhs) {
        if(this.vertices.contains(rhs) || !this.is_independent_of(rhs)) {
            return null;
        }
        return this.extend(rhs);
    }

    public boolean extend_self_if_independent(Integer rhs) {
        if(this.vertices.contains(rhs) || !this.is_independent_of(rhs)) {
            return false;
        }
        this.extend_self(rhs);
        return true;
    }

    public int size() {
        return this.vertices.size();
    }

    public int get_vertex(int index) {
        return (int)this.vertices.get(index);
    }
}
