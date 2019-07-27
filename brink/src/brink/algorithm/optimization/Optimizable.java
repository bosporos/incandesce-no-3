package brink.algorithm.optimization;

public interface Optimizable {
    public double cost();
    public Optimizable neighbouring_solution() throws NoAdjacentSolutionException;
}
