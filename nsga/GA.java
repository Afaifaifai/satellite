import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class GA {
    int[][] access_states;
    double[][] coords, distances;
    public GA(int[][] access_states, double[][] coords, double[][] distances, Settings s) {
        this.access_states = access_states;
        this.coords = coords;
        this.distances = distances;

        Chromosome[] population = new Chromosome[s.num_population];
        for (int i = 0; i < s.num_population; i++)
            population[i] = random_generate();

        Arrays.sort(population, (o1, o2) -> Double.compare(o1.fitness, o2.fitness));
        for (int i = 0; i < s.ga_iteration; i++) {
            
        }
    }

    public Chromosome random_generate() {
        int begin = 0, end = 1;
        List<Integer> path = new ArrayList<>(), waiting = new ArrayList<>();
        path.add(begin);
        for (int i = 2; i < access_states.length; i++) {
            waiting.add(i);
        }

        Collections.shuffle(waiting);
        for (int i : waiting) {
            if (i != end)
                path.add(i);
            else
                break;
        }
        path.add(end);

        return new Chromosome(path);
    }

    class Chromosome {
        List<Integer> path;
        double fitness;
        public Chromosome(List<Integer> path) {
            this.path = path;
            fitness = calculate_fitness();
        }

        private double calculate_fitness() {
            return 0.0;
        }

        private int choice(int[] array) {
            int idx = ThreadLocalRandom.current().nextInt(array.length);
            return array[idx];
        }
    }
}


