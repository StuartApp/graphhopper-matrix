package com.graphhopper.routing.matrix.algorithm;

import com.graphhopper.routing.matrix.BucketEntry;
import com.graphhopper.routing.matrix.DistanceMatrix;
import com.graphhopper.routing.matrix.MatrixEntry;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.*;
import com.graphhopper.storage.index.Snap;
import com.graphhopper.util.DistanceCalcEarth;
import com.graphhopper.util.StopWatch;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.ObjectHeapPriorityQueue;

import java.util.List;
import java.util.PriorityQueue;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public abstract class AbstractManyToMany implements MatrixAlgorithm {


    protected RoutingCHGraph graph;

    protected Weighting weighting;

    protected RoutingCHEdgeExplorer inEdgeExplorer;
    protected RoutingCHEdgeExplorer outEdgeExplorer;
    protected CHEdgeFilter levelEdgeFilter;

    protected Int2ObjectMap<Int2ObjectMap<BucketEntry>> bucket;

    protected boolean alreadyRun = false;

    protected int size;
    protected int maxVisitedNodes = Integer.MAX_VALUE;
    protected int visitedNodes = 0;

    protected int maxNodes;

    protected PriorityQueue<MatrixEntry> heap;
    protected Int2ObjectMap<MatrixEntry> map;
    protected Int2DoubleMap tentativeWeights;

    double[] targetsMaxDistance;
    double[] sourcesMaxDistance;

    double DISTANCE_MULT = 2;

    public AbstractManyToMany(QueryRoutingCHGraph graph){

        this.graph = graph;
        this.weighting = graph.getWrappedWeighting();
        this.inEdgeExplorer = graph.createInEdgeExplorer();
        this.outEdgeExplorer = graph.createOutEdgeExplorer();
        this.maxNodes = graph.getBaseGraph().getBaseGraph().getNodes();
        this.levelEdgeFilter = new CHEdgeFilter() {

            @Override
            public boolean accept(RoutingCHEdgeIteratorState edgeState) {

                // always accept virtual edges, see #288
                if (edgeState.getBaseNode() >= maxNodes || edgeState.getAdjNode() >= maxNodes) return true;

                // minor performance improvement: shortcuts in wrong direction are disconnected, so no need to exclude them
                if (edgeState.isShortcut()) return true;

                return graph.getLevel(edgeState.getBaseNode()) <= graph.getLevel(edgeState.getAdjNode());

            }
        };

        this.size = Math.min(Math.max(200, graph.getNodes() / 10), 150_000);
        this.bucket = new Int2ObjectOpenHashMap<>();
        this.heap = new PriorityQueue<>();
        this.map = new Int2ObjectOpenHashMap<>();
        this.tentativeWeights = new Int2DoubleOpenHashMap();
    }

    private void calculateMaxDistance(List<Snap>  sources, List<Snap> targets){

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

        Int2ObjectMap<IntArrayList> targetIdxsNodes = new Int2ObjectOpenHashMap<>(targets.size());

        StopWatch maxWatch = new StopWatch();
        maxWatch.start();
        calculateMaxDistance(sources,targets);
        maxWatch.stop();
        System.out.println("Time Elapsed Max: " + maxWatch.getMillis());

        StopWatch watch = new StopWatch();
        watch.start();
        //Backward
        int idxTarget =0;
        while(idxTarget < targets.size()){
            int targetClosestNode = targets.get(idxTarget).getClosestNode();

            //Avoid iterate over the same node two times
            if(!targetIdxsNodes.containsKey(targetClosestNode)){
                targetIdxsNodes.put(targetClosestNode,IntArrayList.of(idxTarget));
                backwardSearch(targets.get(idxTarget),idxTarget);
            }else{
                targetIdxsNodes.get(targetClosestNode).add(idxTarget);
            }

            idxTarget++;
        }
        watch.stop();
        System.out.println("Time Elapsed Backward: " + watch.getMillis());

        StopWatch watchf = new StopWatch();
        watchf.start();
        //Forward
        int idxSource =0;
        while(idxSource < sources.size()){
            forwardSearch(sources.get(idxSource),idxSource,matrix,targetIdxsNodes);
            idxSource++;
        }
        watchf.stop();
        System.out.println("Time Elapsed Forward: " + watchf.getMillis());

        return matrix;
    }

    protected void checkAlreadyRun() {
        if (alreadyRun) throw new IllegalStateException("Create a new instance per call");
        alreadyRun = true;
    }

    protected void backwardSearch( Snap targetSnap, int targetIdx){

        IntSet visited = new IntOpenHashSet();

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

                    saveToBucket(entry, iterator,target);
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

    protected void saveToBucket(MatrixEntry entry, RoutingCHEdgeIteratorState iter, int target){
        int node = iter.getAdjNode();

            //if(graph.getLevel(iter.getBaseNode()) < graph.getLevel(node)){

                Int2ObjectMap<BucketEntry> bucketDistances = bucket.get(node);

                if(bucketDistances == null){
                    bucketDistances = new Int2ObjectOpenHashMap<>();
                    bucketDistances.put(target,new BucketEntry(entry.weight,entry.time ,entry.distance));
                    bucket.put(node,bucketDistances);
                }else {

                    BucketEntry targetEntry = bucketDistances.get(target);
                    if (targetEntry == null || targetEntry.weight > entry.weight) {
                        bucketDistances.put(target, new BucketEntry(entry.weight, entry.time, entry.distance));
                    }
                }
           // }
      }

    protected void forwardSearch(Snap sourceSnap, int idxSource, DistanceMatrix dm, Int2ObjectMap<IntArrayList> targets){

        int source = sourceSnap.getClosestNode();

        IntSet visited = new IntOpenHashSet();
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

        double maxDistanceInMeters = sourcesMaxDistance[idxSource] * DISTANCE_MULT;

        boolean run;

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

    private void saveBestPath( int sourceNode,int idxSource,  Int2ObjectMap<IntArrayList> targets,int currNode,
                               double currentEdgeWeight,
                               long currentEdgeTime, double currentEdgeDistance, DistanceMatrix dm){

        Int2ObjectMap<BucketEntry> bucketEntries = bucket.get(currNode);
        if(bucketEntries != null){

            for( Int2ObjectMap.Entry<BucketEntry> next : bucketEntries.int2ObjectEntrySet()){

                int target = next.getIntKey();

                if(sourceNode == target) continue;

                final double savedWeight = tentativeWeights.get(target);
                final double currentWeight = currentEdgeWeight + next.getValue().weight;

                if(savedWeight == 0.0){

                    final long time = currentEdgeTime + next.getValue().time;
                    final double distance = currentEdgeDistance + next.getValue().distance;
                    tentativeWeights.put(target,currentWeight);

                    for(Integer idxNext : targets.get(target)){
                        dm.setCell(idxSource,idxNext,distance,time);
                    }

                } else if(currentWeight < savedWeight){

                    final long time = currentEdgeTime + next.getValue().time;
                    final double distance = currentEdgeDistance + next.getValue().distance;
                    tentativeWeights.put(target,currentWeight);

                    for(Integer idxNext : targets.get(target)){
                        dm.setCell(idxSource,idxNext,distance,time);
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
}