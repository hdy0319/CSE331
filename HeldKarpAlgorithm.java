import java.util.*;
import java.io.*;
import java.nio.file.*;

public class HeldKarpAlgorithm {
    private final int n;
    private final double[][] dist;
    private final Map<Long, double[]> dp = new HashMap<>();
    private final Map<Long, int[]> parent = new HashMap<>();

    public HeldKarpAlgorithm(double[][] dist) {
        this.n = dist.length;
        this.dist = dist;

        if (n >= 64) throw new IllegalArgumentException("n must be <= 63 (long bitmask limit)");
    }

    private void computeDP() {
        dp.put(1L, new double[n]);
        Arrays.fill(dp.get(1L), Double.POSITIVE_INFINITY);
        dp.get(1L)[0] = 0;

        for (long mask = 1; mask < (1L << n); mask++) {
            if ((mask & 1L) == 0) continue;

            for (int last = 0; last < n; last++) {
                if ((mask & (1L << last)) == 0) continue;
                if (mask == (1L << last)) continue;

                long prevMask = mask ^ (1L << last);
                if (!dp.containsKey(prevMask)) continue;

                double[] prevDp = dp.get(prevMask);
                double[] currDp = dp.computeIfAbsent(mask, k -> {
                    double[] row = new double[n];
                    Arrays.fill(row, Double.POSITIVE_INFINITY);
                    return row;
                });

                int[] currParent = parent.computeIfAbsent(mask, k -> new int[n]);

                for (int prev = 0; prev < n; prev++) {
                    if ((prevMask & (1L << prev)) == 0) continue;
                    double cost = prevDp[prev] + dist[prev][last];
                    if (cost < currDp[last]) {
                        currDp[last] = cost;
                        currParent[last] = prev;
                    }
                }
            }
        }
    }

    private List<Integer> reconstructTour() {
        long fullMask = (1L << n) - 1;
        double[] lastCosts = dp.get(fullMask);
        double minCost = Double.POSITIVE_INFINITY;
        int lastCity = -1;

        for (int last = 1; last < n; last++) {
            double cost = lastCosts[last] + dist[last][0];
            if (cost < minCost) {
                minCost = cost;
                lastCity = last;
            }
        }

        List<Integer> tour = new ArrayList<>();
        long mask = fullMask;
        int curr = lastCity;

        while (true) {
            tour.add(curr);
            if (mask == (1L << curr)) break;
            int prev = parent.get(mask)[curr];
            mask ^= (1L << curr);
            curr = prev;
        }

        Collections.reverse(tour);
        tour.add(0);
        return tour;
    }

    public List<Integer> solve() {
        if (n > 25) System.err.println("Warning: n=" + n + " may be too large for exact Held-Karp");
        computeDP();
        return reconstructTour();
    }

    public double tourCost(List<Integer> tour) {
        double c = 0;
        for (int i = 0; i < tour.size() - 1; i++)
            c += dist[tour.get(i)][tour.get(i + 1)];
        return c;
    }

    public static double[][] readTSPLIB(String file) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line; int dim = 0; List<double[]> pts = new ArrayList<>(); boolean read = false;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("DIMENSION")) dim = Integer.parseInt(line.split(":")[1].trim());
                else if (line.equals("NODE_COORD_SECTION")) read = true;
                else if (line.equals("EOF")) break;
                else if (read && !line.isEmpty()) {
                    String[] p = line.split("\\s+");
                    int id = Integer.parseInt(p[0]) - 1;
                    while (pts.size() <= id) pts.add(null);
                    pts.set(id, new double[]{Double.parseDouble(p[1]), Double.parseDouble(p[2])});
                }
            }
            if (dim == 0 || pts.size() != dim)
                throw new IllegalArgumentException("Invalid TSPLIB file");
            double[][] d = new double[dim][dim];
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    d[i][j] = (i == j ? 0 : Math.hypot(
                        pts.get(i)[0] - pts.get(j)[0],
                        pts.get(i)[1] - pts.get(j)[1]));
            return d;
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.err.println("사용법: java -cp . HeldKarpAlgorithm <TSPLIB 파일>");
            return;
        }

        Path inPath = Paths.get(args[0]);
        String fileName = inPath.getFileName().toString();
        int dot = fileName.lastIndexOf('.');
        String stem = (dot == -1) ? fileName : fileName.substring(0, dot);

        double[][] dist = HeldKarpAlgorithm.readTSPLIB(inPath.toString());
        HeldKarpAlgorithm solver = new HeldKarpAlgorithm(dist);

        long t0 = System.nanoTime();
        List<Integer> tour = solver.solve();
        long t1 = System.nanoTime();

        double cost = solver.tourCost(tour);
        double ms = (t1 - t0) / 1e6;

        System.out.println("Tour  : " + tour);
        System.out.printf("Cost  : %.3f%n", cost);
        System.out.printf("Time  : %.3f ms%n", ms);

        Path outDir = Paths.get("result");
        Files.createDirectories(outDir);
        Path outFile = outDir.resolve(stem + "_time.txt");

        try (var w = Files.newBufferedWriter(outFile, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
            w.write(String.format("%.3f%n", ms));
        }
    }
}
