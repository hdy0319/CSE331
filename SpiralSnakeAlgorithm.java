import java.util.*;
import java.io.*;
import java.nio.file.*;

/**
 * Spiral Snake Algorithm (SSA) – lightweight TSP heuristic
 *   1. 중심 좌표 계산 + 대략적 스케일 추정
 *   2. 중심에 가장 가까운 k개 중 랜덤 시작 → 뱀이 나선형으로 회전하며 이동
 *   3. 각도에 가장 근접하고 반경이 적절한 미방문 노드를 "먹음"
 *   4. 마지막 2‑opt 꼬리 정리
 *
 *   메모리 O(n), 시간 O(n²) (KD‑tree 붙이면 O(n log n)까지 단축 가능)
 */
public class SpiralSnakeAlgorithm {

    /*────────────────────────────┐
     │  뱀 모델                   │
     └────────────────────────────*/
    static class Snake {
        final List<Integer> body = new ArrayList<>();
        double angle = 0;               // 현재 각도 (rad)
    }

    /*────────────────────────────┐
     │  필드                      │
     └────────────────────────────*/
    private final int n;
    private final double[][] coord;     // [i][0]=x, [i][1]=y
    private final Snake snake = new Snake();
    private final boolean[] eaten;      // 방문 여부

    private double cx, cy;              // 중심
    private double bbox;                // bounding‑box 한 변 길이
    private final Random rng;

    /*────────────────────────────┐
     │  생성자                    │
     └────────────────────────────*/
    public SpiralSnakeAlgorithm(double[][] coord){
        this(coord, System.nanoTime());
    }
    public SpiralSnakeAlgorithm(double[][] coord, long seed){
        this.n      = coord.length;
        this.coord  = coord;
        this.eaten  = new boolean[n];
        this.rng    = new Random(seed);
        initCenter();
    }

    /*────────────────────────────┐
     │  거리 헬퍼                 │
     └────────────────────────────*/
    private double dist(int i,int j){
        double dx = coord[i][0]-coord[j][0];
        double dy = coord[i][1]-coord[j][1];
        return Math.hypot(dx,dy);
    }

    /*────────────────────────────┐
     │  1. 중심 & 스케일          │
     └────────────────────────────*/
    private void initCenter(){
        double minX=Double.POSITIVE_INFINITY, minY=minX;
        double maxX=Double.NEGATIVE_INFINITY, maxY=maxX;
        for(var p:coord){
            cx+=p[0]; cy+=p[1];
            if(p[0]<minX)minX=p[0]; if(p[0]>maxX)maxX=p[0];
            if(p[1]<minY)minY=p[1]; if(p[1]>maxY)maxY=p[1];
        }
        cx/=n; cy/=n;
        bbox=Math.max(maxX-minX,maxY-minY);
    }

    /*────────────────────────────┐
     │  2. 랜덤 시작점            │
     └────────────────────────────*/
    private int pickStart(){
        // 중심에 가장 가까운 m개(<=5) 후보 중 랜덤 뽑기
        record Node(int idx,double d){}
        PriorityQueue<Node> pq=new PriorityQueue<>(Comparator.comparingDouble(a->-a.d));
        for(int i=0;i<n;i++){
            double d=Math.hypot(coord[i][0]-cx,coord[i][1]-cy);
            pq.offer(new Node(i,d));
            if(pq.size()>5) pq.poll();
        }
        var list=new ArrayList<Node>();
        while(!pq.isEmpty()) list.add(pq.poll());
        return list.get(rng.nextInt(list.size())).idx();
    }

    /*────────────────────────────┐
     │  3. 각도+반경 기반 선택   │
     └────────────────────────────*/
    private int selectNext(double targetAngle,double radius){
        int best=-1; double bestScore=Double.POSITIVE_INFINITY;
        for(int i=0;i<n;i++) if(!eaten[i]){
            double dx=coord[i][0]-cx, dy=coord[i][1]-cy;
            double theta=Math.atan2(dy,dx);
            double angDiff=Math.abs(theta-targetAngle);
            if(angDiff>Math.PI) angDiff=2*Math.PI-angDiff;
            double rad=Math.hypot(dx,dy);
            double score=angDiff*1000 + Math.abs(rad-radius)*0.1; // 가중합
            if(score<bestScore){bestScore=score;best=i;}
        }
        return best;
    }

    /*────────────────────────────┐
     │  4. 뱀 이동                │
     └────────────────────────────*/
    private void generateTour(){
        double radius=bbox*0.01;            // 초깃 반경 = 1 %
        double spiralRate=1.2;              // 한 바퀴마다 반경 확대배수

        int start=pickStart();
        snake.body.add(start); eaten[start]=true;

        while(snake.body.size()<n){
            // Δθ: 방문이 늘어날수록 증가폭 작게
            double deltaTheta=Math.PI/(8 + snake.body.size()/5000.0);
            snake.angle+=deltaTheta;
            if(snake.angle>2*Math.PI){
                snake.angle-=2*Math.PI;
                radius*=spiralRate;
            }
            int next=selectNext(snake.angle,radius);
            if(next==-1) next=findNearest();
            snake.body.add(next); eaten[next]=true;
        }
    }

    /*────────────────────────────┐
     │  5. 보조: 최근접 미방문    │
     └────────────────────────────*/
    private int findNearest(){
        int head=snake.body.get(snake.body.size()-1);
        int best=-1; double bestD=Double.POSITIVE_INFINITY;
        for(int i=0;i<n;i++) if(!eaten[i]){
            double d=dist(head,i);
            if(d<bestD){bestD=d;best=i;}
        }
        return best;
    }

    /*────────────────────────────┐
     │  6. 꼬리 소형 2‑opt       │
     └────────────────────────────*/
    private void tailOpt(List<Integer> tour){
        int tail=15;
        int m=tour.size();
        int start=Math.max(1,m-tail-1);
        boolean improved=true;
        while(improved){
            improved=false;
            for(int i=start;i<m-2;i++)
                for(int j=i+2;j<m-1;j++){
                    double before=dist(tour.get(i-1),tour.get(i))
                                 +dist(tour.get(j),tour.get(j+1));
                    double after =dist(tour.get(i-1),tour.get(j))
                                 +dist(tour.get(i),tour.get(j+1));
                    if(after<before-1e-9){
                        Collections.reverse(tour.subList(i,j+1));
                        improved=true;
                    }
                }
        }
    }

    /*───────────────────────────┐
     │  Public API               │
     └───────────────────────────*/
    public List<Integer> solve(){
        generateTour();
        List<Integer> tour=new ArrayList<>(snake.body);
        tour.add(tour.get(0));            // 닫기
        tailOpt(tour);
        return tour;
    }

    public double tourCost(List<Integer> tour){
        double c=0; for(int i=0;i<tour.size()-1;i++) c+=dist(tour.get(i),tour.get(i+1));
        return c;
    }

    /*────────────────────────────┐
     │  TSPLIB 파서               │
     └────────────────────────────*/
    public static double[][] readTSPLIB(String file) throws IOException{
        try(BufferedReader br=new BufferedReader(new FileReader(file))){
            String line; int dim=0; List<double[]> pts=new ArrayList<>(); boolean read=false;
            while((line=br.readLine())!=null){
                line=line.trim();
                if(line.startsWith("DIMENSION")) dim=Integer.parseInt(line.split(":")[1].trim());
                else if(line.equals("NODE_COORD_SECTION")) read=true;
                else if(line.equals("EOF")) break;
                else if(read&&!line.isEmpty()){
                    String[] p=line.split("\\s+"); int id=Integer.parseInt(p[0])-1;
                    while(pts.size()<=id) pts.add(null);
                    pts.set(id,new double[]{Double.parseDouble(p[1]),Double.parseDouble(p[2])});
                }
            }
            if(dim==0||pts.size()!=dim) throw new IllegalArgumentException("Invalid TSPLIB file");
            return pts.toArray(new double[dim][2]);
        }
    }

    /*───────────────────────────┐
     │  Demo                     │
     └───────────────────────────*/
    public static void main(String[] args) throws Exception{
        if(args.length==0){
            System.err.println("사용법: java -cp . SpiralSnakeAlgorithm <TSPLIB 파일>");
            return; }

        Path inPath=Paths.get(args[0]); String fileName=inPath.getFileName().toString();
        int dot=fileName.lastIndexOf('.'); String stem=(dot==-1)?fileName:fileName.substring(0,dot);

        double[][] coord=readTSPLIB(inPath.toString());
        SpiralSnakeAlgorithm solver=new SpiralSnakeAlgorithm(coord);

        long t0=System.nanoTime();
        List<Integer> tour=solver.solve();
        long t1=System.nanoTime();

        double cost=solver.tourCost(tour);
        double ms=(t1-t0)/1e6;

        // 콘솔 출력
        System.out.println("Tour  : "+tour);        // ← 누락됐던 투어 출력 추가
        System.out.printf("Cost  : %.3f%n",cost);
        System.out.printf("Time  : %.3f ms%n",ms);

        // 결과 저장
        Path outDir=Paths.get("result"); Files.createDirectories(outDir);
        Path outFile=outDir.resolve(stem+"_time.txt");
        try(var w=Files.newBufferedWriter(outFile,StandardOpenOption.CREATE,StandardOpenOption.APPEND)){
            w.write(String.format("%.3f%n",ms)); }
    }
}