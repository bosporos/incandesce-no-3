package brink.algorithm.optimization;

import brink.random.Uniform;

public class SimulatedAnnealer {
    public double temperature;
    public double temperature_gradient;
    public double temperature_min;

    public int neighbours_on_temperature_plateau;

    public SimulatedAnnealer(double temp, double temp_gradient, double temp_min, int plateau_neighbours) {
        this.temperature = temp;
        this.temperature_gradient = temp_gradient;
        this.temperature_min = temp_min;
        this.neighbours_on_temperature_plateau = plateau_neighbours;
    }

    public Optimizable anneal(Optimizable base_solution, OptimizationTarget target) {
        Optimizable best_solution = base_solution;
        Optimizable old_solution = base_solution;

        while(this.temperature > this.temperature_min) {
            for(int i = 0;i < this.neighbours_on_temperature_plateau;i++) {
                try {
                    Optimizable new_solution = old_solution.neighbouring_solution();
                    if(target == OptimizationTarget.MAXIMIZE) {
                        // Maximize the value of cost()
                        if(Math.exp((new_solution.cost() - old_solution.cost()) / this.temperature) > Uniform.uniform(0.0, 1.0)) {
                            old_solution = new_solution;
                            if(old_solution.cost() > best_solution.cost()) {
                                best_solution = old_solution;
                            }
                        }
                    } else if(target == OptimizationTarget.MINIMIZE) {
                        // Minimize the value of cost()
                        if(Math.exp((old_solution.cost() - new_solution.cost()) / this.temperature) > Uniform.uniform(0.0, 1.0)) {
                            old_solution = new_solution;
                            if(old_solution.cost() < best_solution.cost()) {
                                best_solution = old_solution;
                            }
                        }
                    }
                } catch(NoAdjacentSolutionException nase) {
                    i = this.neighbours_on_temperature_plateau;
                    this.temperature = this.temperature_min;
                    // Break the for, which'll take temp = temp_min * temp_grad so that temp < temp_min, forcing an exit out w/ best_solution
                    break;
                }
            }

            this.temperature *= this.temperature_gradient;
        }

        return best_solution;
    }
}
