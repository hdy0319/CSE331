import java.util.*;
import java.io.*;
import java.nio.file.*;

/**
 * Held-Karp Dynamic Programming TSP Algorithm
 *   - Exact solution via DP with bitmasks
 *   - Time complexity: O(n²2ⁿ)
 *   - Space complexity: O(n2ⁿ)
 */
public class HeldKarpAlgorithm {

    /*────────────────────────────┐
     │  Instance fields           │
     └────────────────────────────*/
    private final int n;
    private final double[][] dist;
    private final double[][] dp;      // dp[mask][last] = min cost to visit cities in mask ending at last
    private final int[][] parent;     // parent[mask][last] = previous city in optimal path

    public HeldKarpAlgorithm(double[][] dist){
        this.n=dist.length; this.dist=dist;
        int states=1<<n;
        this.dp=new double[states][n];
        this.parent=new int[states][n];
        for(double[] row:dp) Arrays.fill(row,Double.POSITIVE_INFINITY);
    }

    /*────────────────────────────┐
     │  DP solver                 │
     └────────────────────────────*/
    private void computeDP(){
        // Base case: start from city 0
        dp[1][0]=0;
        
        // Fill DP table
        for(int mask=1;mask<(1<<n);mask++){
            if((mask&1)==0) continue; // must include city 0
            for(int last=0;last<n;last++){
                if((mask&(1<<last))==0) continue;
                if(mask==(1<<last)) continue; // only city 0
                
                // Try all previous cities
                int prevMask=mask^(1<<last);
                for(int prev=0;prev<n;prev++){
                    if((prevMask&(1<<prev))==0) continue;
                    double cost=dp[prevMask][prev]+dist[prev][last];
                    if(cost<dp[mask][last]){
                        dp[mask][last]=cost;
                        parent[mask][last]=prev;
                    }
                }
            }
        }
    }

    /*────────────────────────────┐
     │  Reconstruct tour          │
     └────────────────────────────*/
    private List<Integer> reconstructTour(){
        int fullMask=(1<<n)-1;
        
        // Find best last city
        double minCost=Double.POSITIVE_INFINITY;
        int lastCity=-1;
        for(int last=1;last<n;last++){
            double cost=dp[fullMask][last]+dist[last][0];
            if(cost<minCost){
                minCost=cost;
                lastCity=last;
            }
        }
        
        // Reconstruct path
        List<Integer> tour=new ArrayList<>();
        int mask=fullMask;
        int curr=lastCity;
        
        while(mask!=0){
            tour.add(curr);
            if(mask==(1<<curr)) break;
            int prev=parent[mask][curr];
            mask^=(1<<curr);
            curr=prev;
        }
        
        Collections.reverse(tour);
        tour.add(0); // close the tour
        return tour;
    }

    /*───────────────────────────┐
     │  Public API               │
     └───────────────────────────*/
    public List<Integer> solve(){
        if(n>20) System.err.println("Warning: n="+n+" is large for Held-Karp (exponential complexity)");
        computeDP();
        return reconstructTour();
    }
    
    public double tourCost(List<Integer> tour){ double c=0; for(int i=0;i<tour.size()-1;i++) c+=dist[tour.get(i)][tour.get(i+1)]; return c; }

    /*────────────────────────────┐
     │  TSPLIB (EUCLIDEAN 2‑D)    │
     └────────────────────────────*/
    public static double[][] readTSPLIB(String file) throws IOException{
        try(BufferedReader br=new BufferedReader(new FileReader(file))){
            String line; int dim=0; List<double[]> pts=new ArrayList<>(); boolean read=false;
            while((line=br.readLine())!=null){
                line=line.trim();
                if(line.startsWith("DIMENSION")) dim=Integer.parseInt(line.split(":")[1].trim());
                else if(line.equals("NODE_COORD_SECTION")) read=true;
                else if(line.equals("EOF")) break;
                else if(read && !line.isEmpty()){
                    String[] p=line.split("\\s+"); int id=Integer.parseInt(p[0])-1;
                    while(pts.size()<=id) pts.add(null);
                    pts.set(id,new double[]{Double.parseDouble(p[1]),Double.parseDouble(p[2])});
                }
            }
            if(dim==0||pts.size()!=dim) throw new IllegalArgumentException("Invalid TSPLIB file");
            double[][] d=new double[dim][dim];
            for(int i=0;i<dim;i++) for(int j=0;j<dim;j++) d[i][j]=(i==j?0:Math.hypot(pts.get(i)[0]-pts.get(j)[0],pts.get(i)[1]-pts.get(j)[1]));
            return d; }
    }

    /*───────────────────────────┐
     │  Demo                     │
     └───────────────────────────*/
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
        List<Integer> tour = solver.solve();        // 투어 구하고
        long t1 = System.nanoTime();

        double cost = solver.tourCost(tour);        // 비용 계산
        double ms   = (t1 - t0) / 1e6;

        // 콘솔 출력 — 눈으로 바로 확인
        System.out.println("Tour  : " + tour);
        System.out.printf ("Cost  : %.3f%n", cost);
        System.out.printf ("Time  : %.3f ms%n", ms);

        // 결과 폴더와 파일(누적)
        Path outDir  = Paths.get("result");
        Files.createDirectories(outDir);
        Path outFile = outDir.resolve(stem + "_time.txt");

        try (var w = Files.newBufferedWriter(
                outFile,
                StandardOpenOption.CREATE,
                StandardOpenOption.APPEND)) {
            w.write(String.format("%.3f%n", ms));   // 시간만 한 줄 추가
        }
    }
}