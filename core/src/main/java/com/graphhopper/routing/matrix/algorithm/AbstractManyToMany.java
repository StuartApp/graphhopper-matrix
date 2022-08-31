package com.graphhopper.routing.matrix.algorithm;

import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.cursors.IntCursor;
import com.carrotsearch.hppc.cursors.IntIntCursor;
import com.carrotsearch.hppc.cursors.IntObjectCursor;
import com.carrotsearch.hppc.cursors.ObjectCursor;
import com.graphhopper.coll.GHIntObjectHashMap;
import com.graphhopper.routing.matrix.BucketEntry;
import com.graphhopper.routing.matrix.DistanceMatrix;
import com.graphhopper.routing.matrix.MatrixEntry;
import com.graphhopper.routing.matrix.RankedNode;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.*;
import com.graphhopper.storage.index.Snap;
import com.graphhopper.util.DistanceCalcEarth;
import com.graphhopper.util.StopWatch;
import com.graphhopper.util.shapes.GHPoint;

import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

public abstract class AbstractManyToMany implements MatrixAlgorithm {


    protected RoutingCHGraph graph;

    protected Weighting weighting;

    protected RoutingCHEdgeExplorer inEdgeExplorer;
    protected RoutingCHEdgeExplorer outEdgeExplorer;
    protected CHEdgeFilter levelEdgeFilter;
    protected CHEdgeFilter sbiLevelEdgeFilter;

    protected IntObjectMap<IntObjectMap<BucketEntry>> buckets;

    protected boolean alreadyRun = false;

    protected int size;
    protected int maxVisitedNodes = Integer.MAX_VALUE;
    protected int visitedNodes = 0;

    protected int maxNodes;

    protected IntObjectMap<MatrixEntry> map;
    protected PriorityQueue<MatrixEntry> heap;

    protected IntDoubleMap tentativeWeights;

    protected IntSet visited = new IntHashSet();

    double[] targetsMaxDistance;
    double[] sourcesMaxDistance;

    double DISTANCE_MULT = 1.5;
    double MIN_DISTANCE =  15000; //15km
    double MAX_DISTANCE = 300000; //300km

    public AbstractManyToMany(QueryRoutingCHGraph graph){

        this.graph = graph;
        this.weighting = graph.getWrappedWeighting();
        this.inEdgeExplorer = graph.createInEdgeExplorer();
        this.outEdgeExplorer = graph.createOutEdgeExplorer();
        this.maxNodes = graph.getBaseGraph().getBaseGraph().getNodes();
        this.levelEdgeFilter = new CHEdgeFilter() {

            @Override
            public boolean accept(RoutingCHEdgeIteratorState edgeState) {

                int base = edgeState.getBaseNode();
                int adj = edgeState.getAdjNode();

                // always accept virtual edges, see #288
                if (base >= maxNodes || adj >= maxNodes) return true;

                // minor performance improvement: shortcuts in wrong direction are disconnected, so no need to exclude them
                if (edgeState.isShortcut()) return true;

                return graph.getLevel(base) <= graph.getLevel(adj);

            }
        };

        this.sbiLevelEdgeFilter = new CHEdgeFilter() {

            @Override
            public boolean accept(RoutingCHEdgeIteratorState edgeState) {

                int base = edgeState.getBaseNode();
                int adj = edgeState.getAdjNode();

                if(base == adj) return false;

                // always accept virtual edges, see #288
                if (base >= maxNodes || adj >= maxNodes) return true;

                // minor performance improvement: shortcuts in wrong direction are disconnected, so no need to exclude them
                if (edgeState.isShortcut()) return true;

                return graph.getLevel(base) < graph.getLevel(adj);
            }
        };

        this.size = Math.min(Math.max(200, graph.getNodes() / 10), 150_000);
        this.buckets = new GHIntObjectHashMap<>(size);

        this.map = new GHIntObjectHashMap<>(size);
        this.heap = new PriorityQueue<>(size);
        this.tentativeWeights = new IntDoubleHashMap(100);
    }

    private void calculateMaxDistance(List<Snap>  sources, List<Snap> targets,DistanceMatrix matrix){

        this.targetsMaxDistance = new double[targets.size()];
        this.sourcesMaxDistance = new double[sources.size()];

        DistanceCalcEarth distanceCalc = new DistanceCalcEarth();

        int idxTarget =0;
        while(idxTarget < targets.size()){

            double targetLat = targets.get(idxTarget).getQueryPoint().lat;
            double targetLon = targets.get(idxTarget).getQueryPoint().lon;

            double max = 0;

            int idxSource = 0;
            while(idxSource < sources.size()){

                double sourceLat = sources.get(idxSource).getQueryPoint().lat;
                double sourceLon = sources.get(idxSource).getQueryPoint().lon;

                double meters = distanceCalc.calcDist(targetLat, targetLon, sourceLat, sourceLon);
                if(meters > max){
                    max = meters;
                }

                double maxSource = sourcesMaxDistance[idxSource];
                if(meters > maxSource){
                    sourcesMaxDistance[idxSource] = meters;
                }

                idxSource++;
            }

            targetsMaxDistance[idxTarget]  = max;

            idxTarget++;
        }

    }

    @Override
    public DistanceMatrix calcMatrix(List<Snap>  sources, List<Snap> targets){

        checkAlreadyRun();

        DistanceMatrix matrix = new DistanceMatrix(sources.size(),targets.size());
        IntObjectMap<IntArrayList> targetIdxsNodes = new GHIntObjectHashMap<>(targets.size());

        StopWatch watch = new StopWatch();
        watch.start();
        simultaneousBucketInitialization(targets);
        watch.stop();
        System.out.println("SBI: " + watch.getTimeString());


        calculateMaxDistance(sources,targets,matrix);

        //Backward
        int idxTarget =0;
        while(idxTarget < targets.size()){
            int targetClosestNode = targets.get(idxTarget).getClosestNode();

            //Avoid iterate over the same node two times
            if(!targetIdxsNodes.containsKey(targetClosestNode)){
                IntArrayList a = new IntArrayList();
                a.add(idxTarget);
                targetIdxsNodes.put(targetClosestNode,a);
                //backwardSearch(targets.get(idxTarget), idxTarget);
            }else{
                targetIdxsNodes.get(targetClosestNode).add(idxTarget);
            }

            idxTarget++;
        }

        //Forward
        int idxSource =0;
        while(idxSource < sources.size()){
            forwardSearch(sources.get(idxSource),idxSource,matrix,targetIdxsNodes);
            idxSource++;
        }


        return matrix;
    }

    protected void checkAlreadyRun() {
        if (alreadyRun) throw new IllegalStateException("Create a new instance per call");
        alreadyRun = true;
    }

    protected void backwardSearch( Snap targetSnap, int targetIdx){

        visited.clear();

        int target = targetSnap.getClosestNode();

        this.map.clear();
        this.heap.clear();

        MatrixEntry currEdge = new MatrixEntry(target,0,0,0);
        heap.add(currEdge);

        // For the first step though we need all edges, so we need to ignore this filter.
        CHEdgeFilter tmpEdgeFilter = this.levelEdgeFilter;
        this.levelEdgeFilter = CHEdgeFilter.ALL_EDGES;

        boolean run;

        double maxDistanceInMeters = targetsMaxDistance[targetIdx] * DISTANCE_MULT;
        if(maxDistanceInMeters < MIN_DISTANCE){
            maxDistanceInMeters = MIN_DISTANCE;
        }else if(maxDistanceInMeters > MAX_DISTANCE){
            maxDistanceInMeters = MAX_DISTANCE;
        }

        do {
            visitedNodes++;
            currEdge = heap.poll();

            if(currEdge.distance > maxDistanceInMeters)
                break;

            int currNode = currEdge.adjNode;

            if(visited.contains(currNode)){
                run = !heap.isEmpty();
                continue;
            }

            RoutingCHEdgeIterator iterator = inEdgeExplorer.setBaseNode(currNode);

            while(iterator.next()){

                if(!accept(iterator,currEdge))
                    continue;

                final double weight = calcWeight(iterator,currEdge,true);

                if(!Double.isFinite(weight))
                    continue;

                final int origEdgeId = getOrigEdgeId(iterator, true);
                final int traversalId = getTraversalId(iterator, origEdgeId, true);
                MatrixEntry entry = map.get(traversalId);

                if (entry == null) {

                    final double distance = calcDistance(iterator, currEdge);
                    final long time = calcTime(iterator, currEdge, true);

                    entry = new MatrixEntry(iterator.getEdge(), origEdgeId, iterator.getAdjNode(), weight, time, distance);
                    map.put(traversalId, entry);
                    heap.add(entry);

                    saveToBucket(entry, iterator, target);

                } else if (entry.getWeightOfVisitedPath() > weight) {


                    final double distance = calcDistance(iterator, currEdge);
                    final long time = calcTime(iterator, currEdge, true);

                    heap.remove(entry);
                    entry.edge = iterator.getEdge();
                    entry.weight = weight;
                    entry.distance = distance;
                    entry.time = time;
                    heap.add(entry);

                    saveToBucket(entry, iterator, target);

                }
            }

            this.levelEdgeFilter = tmpEdgeFilter;

            visited.add(currNode);
            run = (visitedNodes <= maxVisitedNodes) && !heap.isEmpty();
        }while(run);
    }

    protected abstract int getTraversalId(RoutingCHEdgeIteratorState edge, int origEdgeId, Boolean reverse);

    protected int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse) {
        return reverse ? edge.getOrigEdgeFirst() : edge.getOrigEdgeLast();
    }

    protected boolean accept(RoutingCHEdgeIteratorState edge, MatrixEntry currEdge) {

        if(edge.getEdge() == getIncomingEdge(currEdge))
            return false;
        else
            return levelEdgeFilter == null || levelEdgeFilter.accept(edge);
    }

    protected abstract double calcWeight(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge, boolean reverse);

    protected abstract long calcTime(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge, boolean reverse);

    protected int getIncomingEdge(MatrixEntry entry) {
        return entry.edge;
    }

    protected abstract double calcDistance(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge);

    private double calcWeight(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId){

        double edgeWeight = edgeState.getWeight(reverse);
        final int origEdgeId = reverse ? edgeState.getOrigEdgeLast() : edgeState.getOrigEdgeFirst();
        double turnCosts = reverse
                ? graph.getTurnWeight(origEdgeId, edgeState.getBaseNode(), prevOrNextEdgeId)
                : graph.getTurnWeight(prevOrNextEdgeId, edgeState.getBaseNode(), origEdgeId);
        return edgeWeight + turnCosts;

    }

    private long calcTime(RoutingCHEdgeIteratorState edgeState, Boolean reverse, int prevOrNextEdgeId){

        long time = edgeState.getTime(reverse);
        int origEdgeId;
        if(reverse){
            origEdgeId = edgeState.getOrigEdgeLast();
        }else{
            origEdgeId = edgeState.getOrigEdgeFirst();
        }
        long turnCost;
        if(reverse){
            turnCost = weighting.calcTurnMillis(origEdgeId,edgeState.getBaseNode(),prevOrNextEdgeId);
        }else{
            turnCost = weighting.calcTurnMillis(prevOrNextEdgeId,edgeState.getBaseNode(),origEdgeId);
        }

        return time + turnCost;
    }


    protected double calcDistance(RoutingCHEdgeIteratorState iter){
        return iter.getDistance();
    }

    protected void saveToBucket(MatrixEntry entry, RoutingCHEdgeIteratorState iter, int target){
        int node = iter.getAdjNode();

        if(node != target){

            IntObjectMap<BucketEntry> bucketDistances = buckets.get(node);

            if(bucketDistances == null){

                bucketDistances = new GHIntObjectHashMap<>();
                bucketDistances.put(target,new BucketEntry(entry.weight,entry.time ,entry.distance,0));
                buckets.put(node,bucketDistances);
            }else {

                BucketEntry targetEntry = bucketDistances.get(target);
                if (targetEntry == null || targetEntry.weight > entry.weight) {
                    bucketDistances.put(target, new BucketEntry(entry.weight, entry.time, entry.distance,0));
                }
            }
        }
    }

    protected void forwardSearch(Snap sourceSnap, int idxSource, DistanceMatrix dm, IntObjectMap<IntArrayList> targets){

        int source = sourceSnap.getClosestNode();

        this.visited.clear();
        this.map.clear();
        this.heap.clear();
        this.tentativeWeights.clear();

        //For the case when the shortest path is the direct path between source and target
        saveBestPath(source,idxSource,targets,source,0,0,0,dm);

        MatrixEntry currEdge = new MatrixEntry(source,0,0,0);
        map.put(source,currEdge);
        heap.add(currEdge);

        // For the first step though we need all edges, so we need to ignore this filter.
        CHEdgeFilter tmpEdgeFilter = this.levelEdgeFilter;
        this.levelEdgeFilter = CHEdgeFilter.ALL_EDGES;

        boolean run;

        double maxDistanceInMeters = sourcesMaxDistance[idxSource] * DISTANCE_MULT;
        if(maxDistanceInMeters < MIN_DISTANCE){
            maxDistanceInMeters = MIN_DISTANCE;
        }else if(maxDistanceInMeters > MAX_DISTANCE){
            maxDistanceInMeters = MAX_DISTANCE;
        }

        do {
            visitedNodes++;

            currEdge = heap.poll();

            if(currEdge.distance > maxDistanceInMeters)
                break;

            int currNode = currEdge.adjNode;

            if(visited.contains(currNode)){
                run = !heap.isEmpty();
                continue;
            }

            final RoutingCHEdgeIterator iterator = outEdgeExplorer.setBaseNode(currNode);
            while (iterator.next()) {

                if (!accept(iterator, currEdge))
                    continue;

                final double weight = calcWeight(iterator, currEdge, false);
                if (!Double.isFinite(weight))
                    continue;

                final int origEdgeId = getOrigEdgeId(iterator, false);
                final int traversalId = getTraversalId(iterator, origEdgeId, false);
                MatrixEntry entry = map.get(traversalId);

                if (entry == null) {

                    final double distance = calcDistance(iterator, currEdge);
                    final long time = calcTime(iterator, currEdge, false);

                    entry = new MatrixEntry(iterator.getEdge(), origEdgeId, iterator.getAdjNode(), weight, time, distance);
                    map.put(traversalId, entry);
                    heap.add(entry);

                    saveBestPath(source,idxSource,targets,iterator.getAdjNode(),weight,time,distance,dm);

                } else if (entry.getWeightOfVisitedPath() > weight) {

                    final double distance = calcDistance(iterator, currEdge);
                    final long time = calcTime(iterator, currEdge, false);

                    heap.remove(entry);
                    entry.edge = iterator.getEdge();
                    entry.weight = weight;
                    entry.distance = distance;
                    entry.time = time;
                    heap.add(entry);

                    saveBestPath(source,idxSource,targets,iterator.getAdjNode(),weight,time,distance,dm);
                }
            }

            this.levelEdgeFilter = tmpEdgeFilter;

            visited.add(currNode);

            run = (visitedNodes <= maxVisitedNodes) && !heap.isEmpty();
        }while(run);

    }

    private void saveBestPath( int sourceNode,int idxSource, IntObjectMap<IntArrayList> targets,int currNode, double currentEdgeWeight,
                               long currentEdgeTime, double currentEdgeDistance, DistanceMatrix dm){

        final IntObjectMap<BucketEntry> bucketEntries = buckets.get(currNode);
        if(bucketEntries != null){

            for( IntObjectCursor<BucketEntry> next : bucketEntries){

                int target = next.key;

                if(sourceNode == target) continue;

                final double savedWeight = tentativeWeights.get(target);
                final double currentWeight = currentEdgeWeight + next.value.weight;

                if(savedWeight == 0.0){

                    final long time = currentEdgeTime + next.value.time;
                    final double distance = currentEdgeDistance + next.value.distance;
                    tentativeWeights.put(target,currentWeight);


                    for(IntCursor idxNext : targets.get(target)){
                        dm.setCell(idxSource,idxNext.value,distance,time);
                    }

                } else if(currentWeight < savedWeight){

                    final long time = currentEdgeTime + next.value.time;
                    final double distance = currentEdgeDistance + next.value.distance;
                    tentativeWeights.put(target,currentWeight);

                    for(IntCursor idxNext : targets.get(target)){
                        dm.setCell(idxSource,idxNext.value,distance,time);
                    }
                }
            }
        }
    }

    @Override
    public int getVisitedNodes() {
        return visitedNodes;
    }

    @Override
    public void setMaxVisitedNodes(int numberOfNodes){
        this.maxVisitedNodes = numberOfNodes;
    }


    /* SBI Methods */

    public class Vertex {
        public int vertex;
        public double weight;
        public long time;
        public double distance;


        public Vertex(int vertex, double weight, long time, double distance) {
            this.vertex = vertex;
            this.weight = weight;
            this.time = time;
            this.distance = distance;
        }
    }

    private void simultaneousBucketInitialization(List<Snap> targets){

        int nodes = graph.getNodes();
        //System.out.println("Graph Nodes:" + nodes);

        IntObjectMap<ObjectArrayList<Vertex>> upVertices = new IntObjectHashMap(nodes);
        IntObjectMap<ObjectArrayList<Vertex>> downVertices = new IntObjectHashMap<>(nodes);

        IntSet targetsNodes = new IntHashSet();
        IntSet nodesAdded = new IntHashSet();
        IntIntHashMap terminals = new IntIntHashMap();
        System.out.println("Targets:" + targets.size());

        double [] distances = new double[targets.size()];
        long [] times = new long[targets.size()];
        double [] weights = new double[targets.size()];

        PriorityQueue<RankedNode> queue = new PriorityQueue<>(size);

        //Targets Initialization
        //we initialize a priority queue Q ordered by minimum rank rank(v) with all targets
        int idxTarget =0;
        while(idxTarget < targets.size()){
            int node = targets.get(idxTarget).getClosestNode();
            int rank = graph.getLevel(node);

            distances[idxTarget] = Double.POSITIVE_INFINITY;
            times[idxTarget] = Long.MAX_VALUE;
            weights[idxTarget] = Double.POSITIVE_INFINITY;

            nodesAdded.add(node);

            //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
            //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.
            IntObjectMap<BucketEntry> bucketTargets = new GHIntObjectHashMap<>();
            bucketTargets.put(node,new BucketEntry(0,0,0, idxTarget));
            buckets.put(node,bucketTargets);

            if(rank == Integer.MAX_VALUE){
                rank = 0;
            }
            queue.add(new RankedNode(node,rank));

            System.out.println("Node Target: " + node + " ->" + rank);

            idxTarget++;
        }

        int nodesCount = 0;

        //Main loop
        while(!queue.isEmpty()){

            //In each iteration, a vertex v with minimum rank is removed from Q and settled
            RankedNode rn = queue.poll();
            int v = rn.node;
            System.out.println("**** Node: " + v +  "(" + graph.getLevel(v) +") ****");
            nodesCount++;

            //Initialize DownVertices
            RoutingCHEdgeIterator downIterator = inEdgeExplorer.setBaseNode(v);
            while(downIterator.next()){

                int u = downIterator.getAdjNode();


                boolean accepted =  this.sbiLevelEdgeFilter.accept(downIterator);

                if(!nodesAdded.contains(u) && accepted){
                    nodesAdded.add(u);
                    int rank = graph.getLevel(u);
                    if(rank == Integer.MAX_VALUE)
                        rank = 0;
                    queue.add(new RankedNode(u,rank));
                }

                if(accepted){

                    double weight = calcWeight(downIterator,true,v);
                    double distance = calcDistance(downIterator);
                    long time = calcTime(downIterator,true,v);

                    ObjectArrayList<Vertex> dVertices = downVertices.get(u);
                    if(dVertices == null){
                        dVertices = new ObjectArrayList<>();
                        downVertices.put(u,dVertices);
                    }

                    System.out.println("In: " + u  + " : " + weight);
                    dVertices.add(new Vertex(v,weight,time,distance));
                }
            }

            //Initialize UpVertices
            /*
            RoutingCHEdgeIterator upIterator = outEdgeExplorer.setBaseNode(v);
            while(upIterator.next()){
                if(this.sbiLevelEdgeFilter.accept(upIterator)){
                    double weight = calcWeight(upIterator,false,v);

                    ObjectArrayList<Vertex> uVertices = upVertices.get(upIterator.getAdjNode());
                    if(uVertices == null){
                        uVertices = new ObjectArrayList<>();
                        upVertices.put(upIterator.getAdjNode(),uVertices);
                    }
                    uVertices.add(new Vertex(v,weight,0,0));
                }
            }

             */

            //Discover bucket entries to copy
            ObjectArrayList<Vertex> downList = downVertices.get(v);
            if(downList != null){
                for(ObjectCursor<Vertex> vertex : downList){

                    int w = vertex.value.vertex;
                    double weight = vertex.value.weight;

                    IntObjectMap<BucketEntry> bucketTargets = buckets.get(w);
                    if(bucketTargets != null){

                        for( IntObjectCursor<BucketEntry> next : bucketTargets){

                            int t = next.key;
                            int targetIdx = next.value.targetIdx;
                            double pathWeight = next.value.weight;

                            terminals.putIfAbsent(t,targetIdx);

                            double currentWeight = pathWeight + weight;
                            if(currentWeight < weights[targetIdx]){
                                weights[targetIdx] = currentWeight;
                                distances[targetIdx] = next.value.distance + vertex.value.distance;
                                times[targetIdx] = next.value.time + vertex.value.time;
                            }
                        }
                    }
                }
            }


            //Apply the restrospective prunning algorithm
            /*
            ObjectArrayList<Vertex> upList = upVertices.get(v);
            if(upList != null){
                for(ObjectCursor<Vertex> vertex : upList){
                    IntObjectMap<BucketEntry> bucketTargets = buckets.get(vertex.value.vertex);
                    if(bucketTargets != null){
                        for( IntObjectCursor<BucketEntry> next : bucketTargets){
                            int target = next.key;
                            if(next.value.weight > (vertex.value.weight) + weights[next.value.targetIdx]){
                                //TODO - Enable prunning
                                //System.out.println("Prunning");
                                //bucketTargets.remove(target);
                            }
                        }
                    }
                }
            }

             */

            //Copy bucket entries to the current vertex
            for(IntIntCursor tCursor : terminals){

                int targetIdx = tCursor.value;
                int target = tCursor.key;
                IntObjectMap<BucketEntry> b = buckets.get(v);

                if(b == null){
                    b = new IntObjectHashMap<>();
                    buckets.put(v,b);
                }
                System.out.println("Bucket: " + target + " ->" + weights[targetIdx]);
                b.put(target,new BucketEntry(weights[targetIdx],times[targetIdx],distances[targetIdx],targetIdx));
                weights[targetIdx] = Double.POSITIVE_INFINITY;
                distances[targetIdx] = Double.POSITIVE_INFINITY;
                times[targetIdx] = Long.MAX_VALUE;
            }

            terminals.clear();

        }
        System.out.println("Nodes Processed:" + nodesCount);
        System.out.println("Bucket Entries: " + buckets.size());

    }
}