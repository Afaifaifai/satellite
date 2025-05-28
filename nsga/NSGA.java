import java.util.*;

public class NSGA {
    /* =========  Genetic Algorithm core  ========= */
    static class GeneticAlgorithm {
        /* GA parameters */
        int popSize = 100;
        double crossoverRate = 0.80;
        double mutationRate  = 0.15;
        int maxGen = 500;

        final Graph g;
        final int start, end;
        final Random rnd = new Random();
        List<Chromosome> population = new ArrayList<>();

        GeneticAlgorithm(Graph g, int s, int e) { this.g = g; this.start = s; this.end = e; }

        /* ----- initialisation ----- */
        void initPopulation() {
            while (population.size() < popSize) {
                Chromosome c = randomChromosome();
                c.evaluate(g);
                population.add(c);
            }
        }

        /** Random feasible path (simple DFS walk). */
        Chromosome randomChromosome() {
            boolean[] visited = new boolean[g.n];
            List<Integer> p = new ArrayList<>();
            p.add(start);
            visited[start] = true;
            int cur = start;

            while (cur != end) {
                List<Integer> candidates = new ArrayList<>();
                for (int nxt : g.adj[cur]) {
                    if (!visited[nxt] || nxt == end) candidates.add(nxt);
                }
                if (candidates.isEmpty()) {
                    // dead‑end → restart
                    return randomChromosome();
                }
                cur = candidates.get(rnd.nextInt(candidates.size()));
                visited[cur] = true;
                p.add(cur);
            }
            return new Chromosome(p);
        }

        /* ----- selection ----- */
        Chromosome tournamentSelect() {
            Chromosome a = population.get(rnd.nextInt(popSize));
            Chromosome b = population.get(rnd.nextInt(popSize));
            return a.fitness < b.fitness ? a : b;
        }

        /* ----- evolution loop ----- */
        Chromosome run() {
            initPopulation();
            for (int gen = 0; gen < maxGen; gen++) {
                List<Chromosome> next = new ArrayList<>();
                while (next.size() < popSize) {
                    Chromosome p1 = tournamentSelect();
                    Chromosome p2 = tournamentSelect();
                    Chromosome child = (rnd.nextDouble() < crossoverRate)
                            ? crossover(p1, p2)
                            : new Chromosome(p1.path);
                    if (rnd.nextDouble() < mutationRate) mutate(child);
                    child.evaluate(g);
                    next.add(child);
                }
                population = next;
            }
            return population.stream().min(Comparator.comparingDouble(c -> c.fitness)).orElse(null);
        }

        /* ----- genetic operators (TODO: customise) ----- */

        /**
         * VARIABLE‑LENGTH crossover.
         * TODO implement Common‑Subpath or other operator that guarantees:
         *  1) no repeated cities  2) path remains connected from start to end.
         * Currently returns the better parent as placeholder.
         */
        Chromosome crossover(Chromosome p1, Chromosome p2) {
            Chromosome better = (p1.fitness < p2.fitness) ? p1 : p2;
            return new Chromosome(better.path);
        }

        /**
         * Mutation operator.
         * TODO implement insertion / deletion / sub‑path replacement, ensuring feasibility.
         */
        void mutate(Chromosome c) {
            // no‑op placeholder
        }
    }

    /* =========  demo main() ========= */
    public static void main(String[] args) {
        Graph g = new Graph(6);
        g.addEdge(0, 1, 2); g.addEdge(1, 2, 2); g.addEdge(2, 5, 2);
        g.addEdge(0, 3, 3); g.addEdge(3, 4, 3); g.addEdge(4, 5, 3);

        GeneticAlgorithm ga = new GeneticAlgorithm(g, 0, 5);
        Chromosome best = ga.run();
        System.out.println("Best found path → " + best);
    }
}



/* =========  Chromosome (variable‑length path)  ========= */
class Chromosome {
    final List<Integer> path;   // includes start & end, cities unique
    double fitness = Double.POSITIVE_INFINITY;

    Chromosome(List<Integer> p) { this.path = new ArrayList<>(p); }

    /** Sum of edge weights along the path. */
    void evaluate(Graph g) {
        double cost = 0.0;
        for (int i = 0; i < path.size() - 1; i++) {
            cost += g.weight(path.get(i), path.get(i + 1));
        }
        this.fitness = cost;
    }

    @Override public String toString() { return path + " : " + fitness; }
}



/* =========  Graph structure  ========= */
class Graph {
    final int n;                      // number of vertices
    final double[][] dist;            // distance matrix (0 = no edge)
    final List<Integer>[] adj;        // adjacency list

    @SuppressWarnings("unchecked")
    Graph(int n) {
        this.n = n;
        dist = new double[n][n];
        adj  = new List[n];
        for (int i = 0; i < n; i++) adj[i] = new ArrayList<>();
    }

    void addEdge(int u, int v, double w) {
        dist[u][v] = dist[v][u] = w;
        adj[u].add(v);
        adj[v].add(u);
    }

    boolean connected(int u, int v) { return dist[u][v] > 0; }
    double  weight   (int u, int v) { return dist[u][v]; }
}