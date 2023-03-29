package com.graphhopper.routing.matrix.algorithm;


import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.procedures.IntObjectProcedure;
import com.carrotsearch.hppc.procedures.ObjectProcedure;
import com.graphhopper.routing.matrix.*;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.storage.*;
import com.graphhopper.storage.index.Snap;
import com.graphhopper.util.PairingUtils;


import java.util.*;

/*
 * Implements Simultaneous Bucket Initialization Many-to-Many algorithm (Theo Wieland)
 * https://i11www.iti.kit.edu/_media/teaching/theses/ba_wieland22.pdf
 */

public class ManyToManySBI implements MatrixAlgorithm {

    protected RoutingCHGraph graph;
    protected RoutingCHGraph graphNoVirtualNodes;
    protected RoutingCHEdgeExplorer inEdgeExplorer;
    protected RoutingCHEdgeExplorer inEdgeExplorerNoVirtual;
    protected RoutingCHEdgeExplorer outEdgeExplorer;
    protected RoutingCHEdgeExplorer outEdgeExplorerNoVirtual;

    protected CHEdgeFilter sbiLevelEdgeFilter;

    protected boolean alreadyRun = false;

    protected int size;
    protected int maxVisitedNodes = Integer.MAX_VALUE;
    protected int visitedNodes = 0;
    protected int maxNodes;
    protected PriorityQueue<RankedNode> heap;
    protected IntObjectMap<ObjectArrayList<DownEdge>> downEdges;
    protected IntObjectMap<ObjectArrayList<UpEdge>> upEdges;
    protected IntObjectMap<IntObjectMap<Bucket>> backwardBuckets;
    protected IntObjectMap<IntObjectMap<Bucket>> forwardBuckets;
    protected LongDoubleHashMap shortRoutes;
    IntSet traversedNodes;


    public ManyToManySBI(QueryRoutingCHGraph graph, RoutingCHGraph graphNoVirtualNodes) {

        this.graph = graph;
        this.inEdgeExplorer = graph.createInEdgeExplorer();
        this.outEdgeExplorer = graph.createOutEdgeExplorer();
        this.maxNodes = graph.getBaseGraph().getBaseGraph().getNodes();

        this.graphNoVirtualNodes = graphNoVirtualNodes;
        this.inEdgeExplorerNoVirtual = graphNoVirtualNodes.createInEdgeExplorer();
        this.outEdgeExplorerNoVirtual = graphNoVirtualNodes.createOutEdgeExplorer();

        this.traversedNodes = new IntHashSet();


        this.sbiLevelEdgeFilter = new CHEdgeFilter() {

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

        this.size = Math.min(Math.max(200, graph.getNodes() / 10), 150_000);
        this.heap = new PriorityQueue<>(size);

        downEdges = new IntObjectScatterMap<>(size);
        upEdges= new IntObjectScatterMap<>(size);
        backwardBuckets = new IntObjectScatterMap<>(size);
        forwardBuckets = new IntObjectScatterMap<>(size);
        shortRoutes = new LongDoubleScatterMap();
    }

    @Override
    public DistanceMatrix calcMatrix(MatrixSnapResult sources, MatrixSnapResult targets) {

        checkAlreadyRun();

        int numTargets = targets.size();
        int numSources = sources.size();
        DistanceMatrix matrix = new DistanceMatrix(numSources, numTargets,sources.getPointsNotFound(),targets.getPointsNotFound());

        //Backward

        int idxTarget = 0;
        while (idxTarget < numTargets) {
            if(targets.isFound(idxTarget)){
                int closestNode = targets.get(idxTarget).getClosestNode();
                findInitialNodesBackward(closestNode, idxTarget,null);
            }
            idxTarget++;
        }

        backward();

        //Reset collections for forward
        this.heap.clear();
        this.traversedNodes.clear();
        this.downEdges.clear();
        this.upEdges.clear();

        //Forward
        int idxSource = 0;
        while (idxSource < numSources) {
            if(sources.isFound(idxSource)) {
                int closestNode = sources.get(idxSource).getClosestNode();
                findInitialNodesForward(closestNode, idxSource, matrix);
            }
            idxSource++;
        }

        forward(matrix);

        return matrix;

    }

    protected void checkAlreadyRun() {
        if (alreadyRun) throw new IllegalStateException("Create a new instance per call");
        alreadyRun = true;
    }

    @Override
    public int getVisitedNodes() {
        return visitedNodes;
    }

    @Override
    public void setMaxVisitedNodes(int numberOfNodes) {
        this.maxVisitedNodes = numberOfNodes;
    }

    @Override
    public String getName() {
        return "many_to_many_sbi";
    }

    private boolean isVirtual(int node) {
        return node >= maxNodes;
    }

    private int getLevel(int node) {
        if (isVirtual(node)) {
            return 0;
        } else {
            return graph.getLevel(node);
        }
    }

    private void addInitialBucket(int node,int idx, IntObjectMap<IntObjectMap<Bucket>> buckets){
        Bucket bucket = new Bucket(idx,0,0,0);
        IntObjectMap<Bucket> bucketsDistance = obtainNodeBuckets(node,buckets);
        bucketsDistance.put(idx,bucket);
        copyBucketsToCurrentVertex(node,bucketsDistance, buckets);
    }

    private  IntObjectMap<Bucket> obtainNodeBuckets(int node, IntObjectMap<IntObjectMap<Bucket>> buckets){
        IntObjectMap<Bucket>  nodeBuckets = buckets.get(node);
        if(nodeBuckets == null){
            nodeBuckets = new IntObjectScatterMap<>();
        }

        return nodeBuckets;
    }

    private void addInitialNode(int idx, int nonVirtualNode, RoutingCHEdgeExplorer explorer,
                                IntObjectMap<IntObjectMap<Bucket>> buckets, boolean reverse,DistanceMatrix dm){

        boolean noAccessibleNodes = true;
        RoutingCHEdgeIterator iter = explorer.setBaseNode(nonVirtualNode);

        while (iter.next()) {

            if(!sbiLevelEdgeFilter.accept(iter)){
                continue;
            }

            int adjNode = iter.getAdjNode();
            double weight = iter.getWeight(reverse);
            boolean isVirtualAdj = isVirtual(adjNode);

            if (weight == Double.POSITIVE_INFINITY || isVirtualAdj) {
                continue;
            }

            noAccessibleNodes = false;
        }

        heap.add(new RankedNode(nonVirtualNode,getLevel(nonVirtualNode),noAccessibleNodes));
        traversedNodes.add(nonVirtualNode);
        addInitialBucket(nonVirtualNode,idx,buckets);
        if(reverse == false){
            discoverInitialShortPaths(idx,nonVirtualNode,dm);
        }
    }

    private void addInitialVirtualNode(int idx,int virtualNode,RoutingCHEdgeExplorer explorer,
                                       IntObjectMap<IntObjectMap<Bucket>> buckets, boolean reverse,
                                       DistanceMatrix dm){

        IntSet processed = new IntHashSet();

        int closestLevel = getLevel(virtualNode);
        RankedNode initial = new RankedNode(virtualNode,closestLevel,false);
        Deque<RankedNode> queue = new ArrayDeque<>();
        queue.add(initial);
        addInitialBucket(virtualNode,idx,buckets);
        if(reverse == false){
            discoverInitialShortPaths(idx,virtualNode,dm);
        }

        this.traversedNodes.add(virtualNode);

        while (!queue.isEmpty()) {

            RankedNode current = queue.poll();

            RoutingCHEdgeIterator iter = explorer.setBaseNode(current.node);
            while (iter.next()) {

                int adjNode = iter.getAdjNode();
                int adjLevel = getLevel(adjNode);

                double weight = iter.getWeight(reverse);
                boolean isVirtualAdj = isVirtual(adjNode);

                if (weight == Double.POSITIVE_INFINITY) {
                    continue;
                }

                double distance = iter.getDistance();
                long time = iter.getTime(reverse);

                addDownVertex(current.node,adjNode,weight,time,distance);

                if (isVirtualAdj && adjNode != virtualNode && !processed.contains(adjNode)) {
                    processed.add(adjNode);
                    queue.add(new RankedNode(adjNode,adjLevel,false));
                } else if (adjNode != virtualNode && !processed.contains(adjNode)) {
                    if(!this.traversedNodes.contains(adjNode)){
                        heap.add(new RankedNode(adjNode,adjLevel,false));
                        this.traversedNodes.add(adjNode);
                    }
                }
            }
            if(current.node != virtualNode){
                IntObjectMap<Bucket> nodeBucketsDistance = obtainNodeBuckets(current.node,buckets);
                discoverBucketsEntriesToCopy(current.node, buckets,nodeBucketsDistance);
                copyBucketsToCurrentVertex(current.node,nodeBucketsDistance, buckets);
            }

        }
    }

    private void findInitialNodesBackward(int closestNode, int idx, DistanceMatrix dm) {
        findInitialNodes(closestNode,idx,inEdgeExplorer,backwardBuckets,true,dm);
    }

    private void findInitialNodesForward(int closestNode, int idx, DistanceMatrix dm) {
        discoverInitialShortPaths(idx,closestNode,dm);
        findInitialNodes(closestNode,idx,outEdgeExplorer,forwardBuckets,false,dm);
    }

    private void findInitialNodes(int closestNode, int idx,RoutingCHEdgeExplorer explorer,
                                  IntObjectMap<IntObjectMap<Bucket>> buckets, boolean reverse, DistanceMatrix dm){

        boolean closestIsVirtual = isVirtual(closestNode);

        if(!closestIsVirtual && !this.traversedNodes.contains(closestNode)){
            addInitialNode(idx,closestNode,explorer,buckets,reverse,dm);
        }else{
            addInitialVirtualNode(idx,closestNode,explorer,buckets,reverse,dm);
        }
    }

    private void backward(){

        while (!heap.isEmpty()) {

            RankedNode current = heap.poll();
            IntObjectMap<Bucket> bucketsDistances = obtainNodeBuckets(current.node,backwardBuckets);
            processDownEdgesForCurrentNode(current,inEdgeExplorerNoVirtual,true,false);
            processUpEdgesForCurrentNode(current,outEdgeExplorerNoVirtual,true);
            discoverBucketsEntriesToCopy(current.node, backwardBuckets,bucketsDistances);
            retrospectivePruning(current.node, backwardBuckets,bucketsDistances);
            copyBucketsToCurrentVertex(current.node,bucketsDistances, backwardBuckets);
        }
    }

    private void forward(DistanceMatrix matrix) {

        while (!heap.isEmpty()) {
            RankedNode current = heap.poll();
            IntObjectMap<Bucket> bucketsDistances = obtainNodeBuckets(current.node,forwardBuckets);
            processDownEdgesForCurrentNode(current,outEdgeExplorerNoVirtual,false,false);
            processUpEdgesForCurrentNode(current,inEdgeExplorerNoVirtual,false);
            discoverBucketsEntriesToCopy(current.node, forwardBuckets,bucketsDistances);
            retrospectivePruning(current.node, forwardBuckets,bucketsDistances);
            copyBucketsToCurrentVertex(current.node,bucketsDistances, forwardBuckets);
            discoverShortPaths(current.node,bucketsDistances,matrix);
        }
    }

    private void processDownEdgesForCurrentNode(RankedNode current,
                                                RoutingCHEdgeExplorer explorer,
                                                boolean reverse, boolean noAccessibleNodes){

        RoutingCHEdgeIterator iter = explorer.setBaseNode(current.node);
        while (iter.next()) {

            if (!this.sbiLevelEdgeFilter.accept(iter) && !noAccessibleNodes) {
                continue;
            }

            double weight = iter.getWeight(reverse);

            if (Double.isInfinite(weight)) {
                continue;
            }

            long time = iter.getTime(reverse);
            double distance = iter.getDistance();

            if(!traversedNodes.contains(iter.getAdjNode()) ){
                this.visitedNodes++;
                RankedNode adj = new RankedNode(iter.getAdjNode(),getLevel(iter.getAdjNode()),false);
                traversedNodes.add(iter.getAdjNode());
                heap.add(adj);
            }

            addDownVertex(current.node,iter.getAdjNode(),weight,time,distance);
        }
    }

    private void processUpEdgesForCurrentNode(RankedNode current,
                                              RoutingCHEdgeExplorer explorer,
                                              boolean reverse){

        RoutingCHEdgeIterator iter = explorer.setBaseNode(current.node);
        while (iter.next()) {

            if (!this.sbiLevelEdgeFilter.accept(iter)) {
                continue;
            }

            double weight = iter.getWeight(reverse);

            if (Double.isInfinite(weight)) {
                continue;
            }

            addUpVertex(current.node,iter.getAdjNode(),weight);
        }
    }

    private void addUpVertex(int node, int adjNode, double weight) {
        ObjectArrayList<UpEdge> ups = upEdges.get(adjNode);
        if(ups == null){
            ups = new ObjectArrayList<>();
            upEdges.put(adjNode,ups);
        }
        ups.add(new UpEdge(node,weight));
    }

    private void addDownVertex(int node, int adjNode, double weight, long time, double distance) {

        ObjectArrayList<DownEdge> downs = downEdges.get(adjNode);
        if(downs == null){
            downs = new ObjectArrayList<>();
            downEdges.put(adjNode,downs);
        }
        downs.add(new DownEdge(node,weight,time,distance));
    }


    protected void discoverBucketsEntriesToCopy(int current, IntObjectMap<IntObjectMap<Bucket>> buckets,
                                                IntObjectMap<Bucket> bucketsDistance) {

        ObjectArrayList<DownEdge> downs = downEdges.get(current);
        if (downs != null) {

            downs.forEach(new ObjectProcedure<DownEdge>() {
                              @Override
                              public void apply(DownEdge edge) {
                                  IntObjectMap<Bucket> edgeBuckets = buckets.get(edge.baseNode);

                                  if (edgeBuckets != null) {
                                      edgeBuckets.forEach(new IntObjectProcedure<Bucket>() {
                                          @Override
                                          public void apply(int terminalId, Bucket bucket) {
                                              saveBucket(terminalId, bucket, edge, bucketsDistance);
                                          }
                                      });
                                  }
                              }
                          }
            );
        }

    }

    private void saveBucket(int terminalId, Bucket bucket, DownEdge edge,
                            IntObjectMap<Bucket> bucketsDistance) {
        Bucket current = bucketsDistance.get(terminalId);
        double newWeight = edge.weight + bucket.weight;

        if(current == null || current.weight > newWeight){

            bucketsDistance.put(terminalId,new Bucket(newWeight,
                    edge.time + bucket.time,
                    edge.distance + bucket.distance,
                    bucket.idx));
        }
    }

    private void retrospectivePruning(int currentNode, IntObjectMap<IntObjectMap<Bucket>> buckets,
                                      IntObjectMap<Bucket> bucketsDistance){

        ObjectArrayList<UpEdge> ups = upEdges.get(currentNode);
        if(ups != null){
            ups.forEach(new ObjectProcedure<UpEdge>() {
                @Override
                public void apply(UpEdge upEdge) {
                    IntObjectMap<Bucket> nodeBuckets = buckets.get(upEdge.baseNode);
                    if(nodeBuckets != null){
                        nodeBuckets.forEach(new IntObjectProcedure<Bucket>() {
                            @Override
                            public void apply(int i, Bucket bucket) {
                                Bucket currentNodeDistance = bucketsDistance.get(i);
                                if(currentNodeDistance != null){
                                    if(bucket.weight > (upEdge.weight + currentNodeDistance.weight)){
                                        nodeBuckets.remove(i);
                                    }
                                }
                            }
                        });
                    }
                }
            });
        }
    }

    private void copyBucketsToCurrentVertex(int node, IntObjectMap<Bucket> terminalsDistance,
                                            IntObjectMap<IntObjectMap<Bucket>> edgesTerminals) {
        if(terminalsDistance.size() > 0){
            edgesTerminals.put(node,terminalsDistance);
        }
    }

    private void discoverShortPaths(int node, IntObjectMap<Bucket> sourcesDistances,DistanceMatrix matrix) {

        IntObjectMap<Bucket> backTerminals = backwardBuckets.get(node);
        if(backTerminals != null){
            backTerminals.forEach(new IntObjectProcedure<Bucket>() {
                @Override
                public void apply(int targetIdx, Bucket target) {
                    sourcesDistances.forEach(new IntObjectProcedure<Bucket>() {
                        @Override
                        public void apply(int sourceIdx, Bucket source) {
                            long key = PairingUtils.pair(sourceIdx,targetIdx);
                            double currentWeight = shortRoutes.getOrDefault(key,Double.MAX_VALUE);
                            double newWeight = source.weight + target.weight;
                            if(currentWeight > newWeight){
                                shortRoutes.put(key,newWeight);
                                double distance = source.distance + target.distance;
                                long time = source.time + target.time;
                                matrix.setCell(sourceIdx,targetIdx,distance,time);
                            }
                        }
                    });

                }
            });
        }
    }

    private void discoverInitialShortPaths(int sourceIdx,int node,DistanceMatrix matrix) {

        IntObjectMap<Bucket> backTerminals = backwardBuckets.get(node);
        if(backTerminals != null){
            backTerminals.forEach(new IntObjectProcedure<Bucket>() {
                @Override
                public void apply(int targetIdx, Bucket target) {
                    long key = PairingUtils.pair(sourceIdx,targetIdx);
                    shortRoutes.put(key,target.weight);
                    matrix.setCell(sourceIdx,targetIdx,target.distance,target.time);
                }
            });
        }
    }
}