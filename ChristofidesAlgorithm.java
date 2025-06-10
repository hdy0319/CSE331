import java.util.*;
import java.io.*;
import java.nio.file.*;

/**
 * Major steps
 *   1. Kruskal MST
 *   2. Collect odd‑degree vertices
 *   3. Minimum‑weight perfect matching on the induced subgraph (Weighted Blossom)
 *   4. Eulerian multigraph → Euler tour → shortcut
 *   5. Return Hamiltonian tour and its cost
 */
public class ChristofidesAlgorithm {

    /*────────────────────────────┐
     │  Basic structures          │
     └────────────────────────────*/
    static class Edge implements Comparable<Edge> {
        int u,v; double w;
        Edge(int u,int v,double w){this.u=u;this.v=v;this.w=w;}
        public int compareTo(Edge o){return Double.compare(w,o.w);} }
    static class UnionFind{int[] p,r;UnionFind(int n){p=new int[n];r=new int[n];for(int i=0;i<n;i++)p[i]=i;}int f(int x){return p[x]==x?x:(p[x]=f(p[x]));}boolean u(int a,int b){a=f(a);b=f(b);if(a==b)return false;if(r[a]<r[b])p[a]=b;else if(r[a]>r[b])p[b]=a;else{p[b]=a;r[a]++;}return true;}}

    /*────────────────────────────┐
     │  Instance fields           │
     └────────────────────────────*/
    private final int n;               // total vertices
    private final double[][] coord;     // metric distances (symmetric)
    private final List<List<Integer>> adj; // multigraph adjacency (duplicates allowed)

    public ChristofidesAlgorithm(double[][] coord){
        this.n = coord.length;
        this.coord = coord;
        this.adj = new ArrayList<>(n);
        for (int i = 0; i < n; i++) adj.add(new ArrayList<>());
    }

    private double dist(int i, int j){
        if (i == j) return 0;
        double dx = coord[i][0] - coord[j][0];
        double dy = coord[i][1] - coord[j][1];
        return Math.hypot(dx, dy);
    }

    /*────────────────────────────┐
     │ 1. Minimum-spanning tree   │
     └────────────────────────────*/
    private void buildMST(){
        // Prim (dense, O(n²) 메모리 O(n))
        boolean[] in = new boolean[n];        // MST 포함 여부
        double[]  key = new double[n];        // 현재 가장 가까운 거리
        int[]     parent = new int[n];        // 연결될 정점

        Arrays.fill(key, Double.POSITIVE_INFINITY);
        key[0] = 0; parent[0] = -1;            // 시작 정점 0

        for (int cnt = 0; cnt < n; cnt++){
            // ① key 최소인 아직-미포함 정점 u 찾기 (O(n))
            int u = -1;
            double best = Double.POSITIVE_INFINITY;
            for (int v = 0; v < n; v++)
                if (!in[v] && key[v] < best){ best = key[v]; u = v; }

            in[u] = true;
            if (parent[u] != -1) addEdge(u, parent[u]);   // MST 간선 추가

            // ② u를 통해 이웃 정점 거리 갱신
            for (int v = 0; v < n; v++) if (!in[v]){
                double d = dist(u, v);        // ← 온-더-플라이 거리 계산
                if (d < key[v]){
                    key[v] = d;
                    parent[v] = u;
                }
            }
        }
    }
    /*────────────────────────────┐
     │ 2. Odd‑degree vertices     │
     └────────────────────────────*/
    private List<Integer> oddVertices(){
        List<Integer> odd=new ArrayList<>();
        for(int v=0;v<n;v++) if(adj.get(v).size()%2==1) odd.add(v);
        return odd;
    }

    /*────────────────────────────┐
     │ 3. MWPM via Blossom        │
     └────────────────────────────*/
    private void perfectMatching(List<Integer> odd){
        if (odd.isEmpty()) return;
        int k = odd.size();
        if ((k & 1) == 1) throw new IllegalStateException("odd set size must be even");
    
        // ── Blossom에 on-the-fly 거리 함수 넘기기 ──
        WeightedBlossom wb = new WeightedBlossom(
            k,
            (i, j) -> dist(odd.get(i), odd.get(j))     // lambda = WeightFn
        );
    
        int[] mate = wb.solve();
        boolean[] done = new boolean[k];
        for (int i = 0; i < k; i++) if (!done[i]){
            int j = mate[i];
            addEdge(odd.get(i), odd.get(j));
            done[i] = done[j] = true;
        }
    }

    /*────────────────────────────┐
     │ 4. Euler & shortcut        │
     └────────────────────────────*/
    private List<Integer> euler(){
        List<Deque<Integer>> g=new ArrayList<>(n);
        for(int i=0;i<n;i++)g.add(new ArrayDeque<>(adj.get(i)));
        List<Integer> circuit=new ArrayList<>();
        Deque<Integer> st=new ArrayDeque<>();st.push(0);
        while(!st.isEmpty()){
            int v=st.peek();
            if(g.get(v).isEmpty()) circuit.add(st.pop());
            else{
                int u=g.get(v).pop(); g.get(u).remove(v); st.push(u);
            }
        }
        Collections.reverse(circuit); return circuit;
    }
    private List<Integer> shortcut(List<Integer> euler){
        boolean[] seen=new boolean[n]; List<Integer> tour=new ArrayList<>();
        for(int v:euler) if(!seen[v]){tour.add(v);seen[v]=true;}
        tour.add(tour.get(0)); return tour; }
    private void addEdge(int u,int v){adj.get(u).add(v);adj.get(v).add(u);}    

    /*───────────────────────────┐
     │  Public API               │
     └───────────────────────────*/
    public List<Integer> solve(){
        buildMST();
        List<Integer> odd=oddVertices();
        perfectMatching(odd);
        return shortcut(euler());
    }
    
    public double tourCost(List<Integer> tour){
        double c = 0;
        for (int i = 0; i < tour.size() - 1; i++)
            c += dist(tour.get(i), tour.get(i + 1));   // ⟵ 여기만 교체
        return c;
    }

    /*────────────────────────────┐
     │  TSPLIB (EUCLIDEAN 2‑D)    │
     └────────────────────────────*/
    public static double[][] readTSPLIB(String file) throws IOException{
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line; int dim = 0; List<double[]> pts = new ArrayList<>(); boolean read = false;
            while ((line = br.readLine()) != null){
                line = line.trim();
                if (line.startsWith("DIMENSION")) dim = Integer.parseInt(line.split(":")[1].trim());
                else if (line.equals("NODE_COORD_SECTION")) read = true;
                else if (line.equals("EOF")) break;
                else if (read && !line.isEmpty()){
                    String[] p = line.split("\\s+");
                    int id = Integer.parseInt(p[0]) - 1;
                    while (pts.size() <= id) pts.add(null);
                    pts.set(id, new double[]{Double.parseDouble(p[1]), Double.parseDouble(p[2])});
                }
            }
            if (dim == 0 || pts.size() != dim)
                throw new IllegalArgumentException("Invalid TSPLIB file");
            return pts.toArray(new double[dim][2]);   // ⟵ 좌표만 리턴
        }
    }

    /*───────────────────────────┐
     │  Demo                     │
     └───────────────────────────*/
    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            System.err.println("사용법: java -cp . ChristofidesRunner <TSPLIB 파일>");
            return;
        }

        Path inPath = Paths.get(args[0]);
        String fileName = inPath.getFileName().toString();
        int dot = fileName.lastIndexOf('.');
        String stem = (dot == -1) ? fileName : fileName.substring(0, dot);

        double[][] dist = ChristofidesAlgorithm.readTSPLIB(inPath.toString());
        ChristofidesAlgorithm solver = new ChristofidesAlgorithm(dist);

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

    /*──────────────────────────────────────────────┐
     │  Weighted Blossom (O(N³)) – full version     │
     └──────────────────────────────────────────────*/
    /**
     * Complete weighted Edmonds–Gabow blossom algorithm (dense graph, symmetric weights).
     *
     * This implementation follows the structure of Kolmogorov’s Blossom V (2015) but is
     * fully self‑contained and trimmed for clarity.  It guarantees a perfect matching for
     * even |V| when the metric is finite and symmetric; otherwise it throws.
     *
     * Complexity:   O(N³)  time,  O(N²)  memory (dense).
     */
    @FunctionalInterface
    interface WeightFn { double w(int i, int j); }

    static class WeightedBlossom {
        /*── Core fields ──*/
        private final int N;
        private final WeightFn W;          // ← 함수로 대체
        private final int[] mate, label, parent, base, Q;
        private final boolean[] inq, blossom, seen;
        private final double[] dual;
        private static final double INF = 1e100, EPS = 1e-7;
        private int qh, qt;

        /*── Init ──*/
        WeightedBlossom(int n, WeightFn w) {
            N = n;  W = w;
            if ((N & 1) == 1) throw new IllegalArgumentException("N must be even");
            mate = new int[N]; Arrays.fill(mate, -1);
            label = new int[N]; parent = new int[N]; base = new int[N];
            inq = new boolean[N]; blossom = new boolean[N]; seen = new boolean[N];
            Q = new int[N];   dual = new double[N];
            for (int i = 0; i < N; i++) dual[i] = minIncident(i);
        }

        /*── Helpers ──*/
        private double minIncident(int v){
            double m = INF;
            for (int u = 0; u < N; u++) if (u != v) m = Math.min(m, W.w(v,u));
            return m;
        }
        private double redCost(int u,int v){ return W.w(u,v) - dual[u] - dual[v]; }
        private boolean tight(int u,int v){ return Math.abs(redCost(u,v)) < EPS; }
        private void enqueue(int v){ Q[qt++] = v; inq[v] = true; label[v] = 0; }
    
        /*──────── Find Lowest Common Ancestor (of blossoms) ────────*/
        private int lca(int x, int y) {
            Arrays.fill(seen, false);
            while (true) {
                x = base[x]; seen[x] = true;
                if (mate[x] == -1) break;
                x = parent[mate[x]];
            }
            while (true) {
                y = base[y];
                if (seen[y]) return y;
                y = parent[mate[y]];
            }
        }
    
        /*──────── Mark blossom from v to base b, w is the child on the path ────────*/
        private void markBlossom(int v, int b, int w) {
            while (base[v] != b) {
                blossom[base[v]] = blossom[base[mate[v]]] = true;
                parent[v] = w; w = mate[v]; v = parent[w];
            }
        }
    
        /*──────── Contract blossom with base b ────────*/
        private void contract(int v, int u) {
            int b = lca(v, u);
            Arrays.fill(blossom, false);
            markBlossom(v, b, u);
            markBlossom(u, b, v);
            for (int i = 0; i < N; i++) if (blossom[base[i]]) {
                base[i] = b;
                parent[i] = -1; // reset parent
                if (!inq[i]) enqueue(i);
            }
        }
    
        /*──────── Augment along the path from root to vertex v ────────*/
        private void augmentPath(int v) {
            while (v != -1) {
                int pv = parent[v];
                int nv = (pv == -1 ? -1 : mate[pv]);
                mate[v] = pv;
                if (pv != -1) mate[pv] = v;
                v = nv;
            }
        }
    
        /*── BFS ──*/
        private boolean bfs(int root){
            Arrays.fill(label,-1); Arrays.fill(parent,-1);
            for(int i=0;i<N;i++) base[i]=i;
            qh=qt=0; Arrays.fill(inq,false);
            enqueue(root);
            while(qh<qt){
                int v=Q[qh++];
                for(int u=0;u<N;u++) if(u!=v&&base[v]!=base[u]&&mate[v]!=u&&tight(v,u)){
                    if(label[u]==-1){
                        label[u]=1; parent[u]=v;
                        if(mate[u]==-1){ augmentPath(u); return true; }
                        label[mate[u]]=0; enqueue(mate[u]);
                    }else if(label[u]==0){ contract(v,u); }
                }
            }
            return false;
        }
    
        /*──────── Dual variable adjustment (Hungarian style) ────────*/
        private void dualAdjust() {
            double delta = INF;
            // δ1 : outer–outer (다른 트리)  rc/2
            for (int v = 0; v < N; v++) if (label[v] == 0)
                for (int u = 0; u < N; u++) if (label[u] == 0 &&
                        base[v] != base[u] && mate[v] != u) {
                    double rc = redCost(v, u);
                    if (rc > EPS)              // ← 필터
                        delta = Math.min(delta, rc * 0.5);
                }
            
            // δ2 : outer–unreached
            for (int v = 0; v < N; v++) if (label[v] == 0)
                for (int u = 0; u < N; u++) if (label[u] == -1 &&
                        base[v] != base[u]) {
                    double rc = redCost(v, u);
                    if (rc > EPS)              // ← 필터
                        delta = Math.min(delta, rc);
                }
            
            // δ3 : inner vertex ↔ mate(inner)
            for (int v = 0; v < N; v++) if (label[v] == 1) {
                int u = mate[v];
                double rc = redCost(v, u);
                if (rc > EPS)                  // ← 필터
                    delta = Math.min(delta, rc);
            }

            if (delta >= INF/2) throw new IllegalStateException("No perfect matching");
        
            // 듀얼 업데이트
            for (int v = 0; v < N; v++) {
                switch (label[v]) {
                    case 0  -> dual[v] += delta;      // outer
                    case 1  -> dual[v] -= delta;      // inner
                    default -> {}                     // unreached
                }
            }
            if (delta < EPS) {
                // 수치 오류나 δ 종류 누락이 의심될 때 강제 중단
                throw new IllegalStateException("δ collapsed to zero; check dualAdjust()");
            }

        }
    
        /*──────── Main solver ────────*/
        int[] solve() {
            int matched = 0;
            int iteration = 0;
            while (matched < N) {
                System.out.println("Blossom iteration " + (++iteration) + ", matched: " + matched + "/" + N);
                boolean progress = false;
                for (int v = 0; v < N; v++) if (mate[v] == -1) {
                    if (bfs(v)) { 
                        matched += 2; 
                        progress = true; 
                    }
                }
                if (!progress) {
                    dualAdjust();
                }
            }
            return mate;
        }
    }
}