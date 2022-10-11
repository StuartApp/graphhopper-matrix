package com.graphhopper.routing.matrix.algorithm;

import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.procedures.IntIntProcedure;
import com.carrotsearch.hppc.procedures.IntObjectProcedure;
import com.carrotsearch.hppc.procedures.IntProcedure;
import com.carrotsearch.hppc.procedures.ObjectProcedure;
import com.graphhopper.coll.GHIntObjectHashMap;
import com.graphhopper.routing.matrix.*;
import com.graphhopper.routing.querygraph.QueryRoutingCHGraph;
import com.graphhopper.routing.weighting.Weighting;
import com.graphhopper.storage.*;
import com.graphhopper.storage.index.Snap;
import com.graphhopper.util.PairingUtils;

import java.util.*;

/*
 * Implements Simultaneous Bucket Initialization Many-to-Many algorithm (Theo Wieland)
 * https://i11www.iti.kit.edu/_media/teaching/theses/ba_wieland22.pdf
 */

public abstract class AbstractManyToMany implements MatrixAlgorithm {


    protected RoutingCHGraph graph;
    protected RoutingCHGraph graphNoVirtualNodes;

    protected Weighting weighting;
    protected Weighting noVirtualWeighting;

    protected RoutingCHEdgeExplorer inEdgeExplorer;
    protected RoutingCHEdgeExplorer inEdgeExplorerNoVirtual;
    protected RoutingCHEdgeExplorer outEdgeExplorer;
    protected RoutingCHEdgeExplorer outEdgeExplorerNoVirtual;

    protected CHEdgeFilter sbiLevelEdgeFilter;

    protected IntObjectMap<IntObjectMap<BucketEntry>> backwardBuckets;
    protected IntObjectMap<IntObjectMap<BucketEntry>> forwardBuckets;

    IntObjectMap<MatrixEntry> bestWeightMap;

    protected boolean alreadyRun = false;

    protected int size;
    protected int maxVisitedNodes = Integer.MAX_VALUE;
    protected int visitedNodes = 0;

    protected int maxNodes;

    protected PriorityQueue<MatrixEntry> heap;

    protected LongDoubleMap tentativeWeights;

    protected IntSet visited = new IntHashSet();

    IntObjectMap<IntObjectHashMap<Vertex>> inVertices;
    IntObjectMap<ObjectArrayList<PruningVertex>> prunningVertices;
    IntObjectMap<IntObjectHashMap<Vertex>> outVertices;

    IntObjectHashMap<Vertex> shortcutSelfVertices;

    IntSet nodesAdded;
    IntSet traversed;
    IntIntHashMap terminals;

    public AbstractManyToMany(QueryRoutingCHGraph graph, RoutingCHGraph graphNoVirtualNodes) {

        this.graph = graph;
        this.weighting = graph.getWrappedWeighting();
        this.inEdgeExplorer = graph.createInEdgeExplorer();
        this.outEdgeExplorer = graph.createOutEdgeExplorer();
        this.maxNodes = graph.getBaseGraph().getBaseGraph().getNodes();

        this.graphNoVirtualNodes = graphNoVirtualNodes;
        this.noVirtualWeighting = graphNoVirtualNodes.getWeighting();
        this.inEdgeExplorerNoVirtual = graphNoVirtualNodes.createInEdgeExplorer();
        this.outEdgeExplorerNoVirtual = graphNoVirtualNodes.createOutEdgeExplorer();

        this.nodesAdded = new IntHashSet();
        this.traversed = new IntHashSet();
        this.terminals = new IntIntHashMap();

        this.bestWeightMap = new IntObjectHashMap<>();

        this.sbiLevelEdgeFilter = new CHEdgeFilter() {

            @Override
            public boolean accept(RoutingCHEdgeIteratorState edgeState) {

                int base = edgeState.getBaseNode();
                int adj = edgeState.getAdjNode();

                // always accept virtual edges, see #288
                if (base >= maxNodes || adj >= maxNodes) return true;

                // minor performance improvement: shortcuts in wrong direction are disconnected, so no need to exclude them
                if (edgeState.isShortcut()) return true;

                if (base == adj) return false;

                return graph.getLevel(base) <= graph.getLevel(adj);
            }
        };

        this.size = Math.min(Math.max(200, graph.getNodes() / 10), 150_000);
        this.backwardBuckets = new GHIntObjectHashMap<>(size);
        this.forwardBuckets = new GHIntObjectHashMap<>(size);

        this.heap = new PriorityQueue<>(size);
        this.tentativeWeights = new LongDoubleHashMap(100);

        this.prunningVertices = new IntObjectHashMap<>(this.graph.getNodes());
        this.inVertices = new IntObjectHashMap<>(this.graph.getNodes());
        this.outVertices = new IntObjectHashMap<>(this.graph.getNodes());
        this.shortcutSelfVertices = new IntObjectHashMap<>();

    }

    @Override
    public DistanceMatrix calcMatrix(List<Snap> sources, List<Snap> targets) {

        checkAlreadyRun();

        DistanceMatrix matrix = new DistanceMatrix(sources.size(), targets.size());
        IntObjectMap<IntArrayList> targetIdxsNodes = new GHIntObjectHashMap<>(targets.size());

        //Backward
        int idxTarget = 0;
        while (idxTarget < targets.size()) {
            int targetClosestNode = targets.get(idxTarget).getClosestNode();

            //Find Buckets
            findInitialNodesBackward(targets.get(idxTarget), idxTarget);

            //Avoid iterate over the same node two times
            if (!targetIdxsNodes.containsKey(targetClosestNode)) {
                IntArrayList a = new IntArrayList();
                a.add(idxTarget);
                targetIdxsNodes.put(targetClosestNode, a);
            } else {
                targetIdxsNodes.get(targetClosestNode).add(idxTarget);
            }

            idxTarget++;
        }

        System.out.println("############################## BACKWARD");

        simultaneousBucketInitializationBackward(targets);

        //Reset collections for forward
        this.heap.clear();
        this.nodesAdded.clear();
        this.terminals.clear();
        this.prunningVertices.clear();

        //Forward
        int idxSource = 0;
        IntIntMap processedSources = new IntIntHashMap();
        List<IntIntPair> cloneResults = new ArrayList<>();
        while (idxSource < sources.size()) {
            //If we have n sources pointing to the same closestNode, we don't want to calculate both, we will copy the result of
            //one to the other
            int closestNode = sources.get(idxSource).getClosestNode();
            if (processedSources.containsKey(closestNode)) {
                cloneResults.add(new IntIntPair(processedSources.get(closestNode), idxSource));
            } else {
                findInitialNodesForward(sources.get(idxSource), idxSource, targetIdxsNodes, matrix);
                processedSources.put(closestNode, idxSource);
            }

            idxSource++;
        }

        System.out.println("############################## FORWARD");

        simultaneousBucketInitializationForward(sources, matrix, targetIdxsNodes);

        cloneResults.stream().forEach(pair -> matrix.copyResult(pair.source, pair.target));

        return matrix;
    }

    protected void checkAlreadyRun() {
        if (alreadyRun) throw new IllegalStateException("Create a new instance per call");
        alreadyRun = true;
    }

    protected int getOtherNode(int edge, int node) {
        return graph.getBaseGraph().getOtherNode(edge, node);
    }

    protected abstract int getTraversalId(RoutingCHEdgeIteratorState edge, Boolean reverse);

    protected abstract int getOrigEdgeId(RoutingCHEdgeIteratorState edge, boolean reverse);

    protected abstract double calcWeight(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                                         boolean reverse, boolean accumulate);

    protected abstract long calcTime(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                                     boolean reverse, boolean accumulate);

    protected abstract double calcDistance(RoutingCHEdgeIteratorState iter, MatrixEntry currEdge,
                                           boolean accumulate);

    @Override
    public int getVisitedNodes() {
        return visitedNodes;
    }

    @Override
    public void setMaxVisitedNodes(int numberOfNodes) {
        this.maxVisitedNodes = numberOfNodes;
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

    private void findInitialNodesBackward(Snap snap, int idx) {

        int closestNode = snap.getClosestNode();
        IntSet processed = new IntHashSet();

        //System.out.println("Backward Init ********************** : " + closestNode);

        //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
        //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.
        IntObjectMap<BucketEntry> bucketTargets = new GHIntObjectHashMap<>();
        bucketTargets.put(closestNode, new BucketEntry(0, 0, 0, idx));
        backwardBuckets.put(closestNode, bucketTargets);

        MatrixEntry initial = new MatrixEntry(closestNode, 0, 0, 0, getLevel(closestNode));
        Deque<MatrixEntry> queue = new ArrayDeque<>();
        queue.add(initial);

        while (!queue.isEmpty()) {

            MatrixEntry current = queue.poll();
            int baseNode = current.adjNode;

            RoutingCHEdgeIterator downIterator = inEdgeExplorer.setBaseNode(baseNode);
            while (downIterator.next()) {

                int adjNode = downIterator.getAdjNode();
                int origEdgeId = getOrigEdgeId(downIterator, true);
                int traversalId = getTraversalId(downIterator, true);

                double weight = calcWeight(downIterator, current, true, true);

                if (weight == Double.POSITIVE_INFINITY) {
                    continue;
                }

                boolean isVirtualAdj = isVirtual(adjNode);
                double distance = calcDistance(downIterator, current, true);
                long time = calcTime(downIterator, current, true, true);

                MatrixEntry entry = new MatrixEntry(downIterator.getEdge(), origEdgeId, adjNode, current.adjNode, getLevel(adjNode), weight, time, distance);

                if (isVirtualAdj && adjNode != closestNode && !processed.contains(adjNode)) {
                    processed.add(adjNode);
                    queue.add(entry);
                    addVertices(inVertices, closestNode, adjNode, weight, time, distance,downIterator.isShortcut());
                } else if (adjNode != closestNode && !processed.contains(adjNode)) {
                    heap.add(entry);
                    nodesAdded.add(traversalId);
                    addVertices(inVertices, closestNode, adjNode, weight, time, distance,downIterator.isShortcut());
                }
            }
        }

    }

    private void findInitialNodesForward(Snap snap, int idx,
                                         IntObjectMap<IntArrayList> targets, DistanceMatrix dm) {

        int closestNode = snap.getClosestNode();
        IntSet processed = new IntHashSet();

        //During the initialization, we add pair (t, 0) to B(t) for each t ∈ T. By adding
        //such a pair, we indicate that t can be reached from t, with a shortest path of length zero.
        IntObjectMap<BucketEntry> bucketTargets = new GHIntObjectHashMap<>();
        bucketTargets.put(closestNode, new BucketEntry(0, 0, 0, idx));
        forwardBuckets.put(closestNode, bucketTargets);

        MatrixEntry initial = new MatrixEntry(closestNode, 0, 0, 0, getLevel(closestNode));

        Deque<MatrixEntry> queue = new ArrayDeque<>();
        queue.add(initial);

        while (!queue.isEmpty()) {

            MatrixEntry current = queue.poll();
            int baseNode = current.adjNode;

            RoutingCHEdgeIterator outIterator = outEdgeExplorer.setBaseNode(baseNode);
            while (outIterator.next()) {

                int adjNode = outIterator.getAdjNode();
                int traversalId = getTraversalId(outIterator, false);
                int origEdgeId = getOrigEdgeId(outIterator, false);

                boolean isVirtualAdj = isVirtual(adjNode);

                double weight = calcWeight(outIterator, current, false, true);

                if (weight == Double.POSITIVE_INFINITY) {
                    continue;
                }

                if (targets.containsKey(adjNode)) {
                    long uniqueId = PairingUtils.pair(closestNode, adjNode);

                    final double savedWeight = tentativeWeights.get(uniqueId);

                    if ((savedWeight == 0.0 || (weight < savedWeight)) && closestNode != adjNode) {

                        double distance = calcDistance(outIterator, current, true);
                        long time = calcTime(outIterator, current, false, true);

                        tentativeWeights.put(uniqueId, weight);

                        targets.get(adjNode).forEach(new IntProcedure() {
                            @Override
                            public void apply(int i) {
                                dm.setCell(idx, i, distance, time);
                            }
                        });
                    }
                }

                double distance = calcDistance(outIterator, current, true);
                long time = calcTime(outIterator, current, false, true);

                MatrixEntry entry = new MatrixEntry(outIterator.getEdge(), origEdgeId, adjNode, current.adjNode, getLevel(adjNode), weight, time, distance);

                if (isVirtualAdj && adjNode != closestNode && !processed.contains(adjNode)) {
                    processed.add(adjNode);
                    queue.add(entry);
                } else if (adjNode != closestNode && !processed.contains(adjNode)) {
                    this.heap.add(entry);
                    //System.out.println("Add Heap:" + entry);
                    this.nodesAdded.add(traversalId);
                    addVertices(outVertices, closestNode, adjNode, weight, time, distance, outIterator.isShortcut());
                }
            }
        }
    }

    private void initializeVertices(MatrixEntry currentEdge, RoutingCHEdgeExplorer explorer,
                                    IntObjectMap<IntObjectHashMap<Vertex>> vertices, boolean reverse) {

        RoutingCHEdgeIterator iter = explorer.setBaseNode(currentEdge.adjNode);
        System.out.println("################ Current:" + currentEdge.adjNode + " - " + currentEdge.edge);

        while (iter.next()) {

            if (!this.sbiLevelEdgeFilter.accept(iter)) {
                continue;
            }

            final double weight = calcWeight(iter, currentEdge, reverse, false);
            final int origEdgeId = getOrigEdgeId(iter, reverse);
            final int traversalId = getTraversalId(iter, reverse);

            if (Double.isInfinite(weight)) {
                continue;
            }

            if(traversed.contains(traversalId)){
                continue;
            }


            int adjNode = iter.getAdjNode();
            int adjRank = getLevel(adjNode);

            double distance = calcDistance(iter, currentEdge, false);
            long time = calcTime(iter, currentEdge, reverse, false);

            System.out.println("Edge:" + iter.getEdge());
            MatrixEntry entry = new MatrixEntry(iter.getEdge(), origEdgeId, adjNode, currentEdge.adjNode, adjRank, weight, time, distance);
            heap.add(entry);
            addVertices(vertices, currentEdge.adjNode, adjNode, weight, time, distance, iter.isShortcut());
            traversed.add(traversalId);
        }
    }

    private void addVertices(IntObjectMap<IntObjectHashMap<Vertex>> vertices,
                             int baseNode, int adjNode, double weight, long time, double distance,
                             boolean shortcut) {

        if(baseNode == adjNode && shortcut){
            Vertex v = shortcutSelfVertices.get(baseNode);
            if(v == null || v.weight > weight){
                System.out.println("Add Self Shortcut Vertex " + baseNode + " " + weight);
                shortcutSelfVertices.put(adjNode,new Vertex(baseNode, weight, time, distance));
            }

        }else{
            IntObjectHashMap<Vertex> dVertices = vertices.get(adjNode);
            if (dVertices == null) {
                dVertices = new IntObjectHashMap<>();
                vertices.put(adjNode, dVertices);
            }

            Vertex vertex = dVertices.get(baseNode);
            if(vertex == null){
                vertex = new Vertex(baseNode, weight, time, distance);
                dVertices.put(baseNode,vertex);
                System.out.println("Add vertex " + baseNode + "-->" + adjNode + " " + weight);
            }else if(vertex.weight > weight){
                System.out.println("Add vertex " + baseNode + "-->" + adjNode + " " + weight);
                vertex = new Vertex(baseNode, weight, time, distance);
                dVertices.put(baseNode,vertex);
            }
        }
    }

    private void initializePruningVertices(int baseNode, RoutingCHEdgeExplorer explorer, boolean reverse) {

        /*
        RoutingCHEdgeIterator upIterator = explorer.setBaseNode(baseNode);
        while (upIterator.next()) {
            int adjNode = upIterator.getAdjNode();
            boolean accept = this.sbiLevelEdgeFilter.accept(upIterator);

            if (accept) {
                ObjectArrayList<PruningVertex> uVertices = prunningVertices.get(adjNode);
                if (uVertices == null) {
                    uVertices = new ObjectArrayList<>();
                    prunningVertices.put(adjNode, uVertices);
                }

                double weight = calcWeight(upIterator, reverse, baseNode);
                if (weight < Double.POSITIVE_INFINITY) {
                    uVertices.add(new PruningVertex(baseNode, weight));
                }
            }
        }

         */
    }

    private void applyRetrospectivePruningAlgorithm(int baseNode, double[] weights,
                                                    IntObjectMap<IntObjectMap<BucketEntry>> buckets) {

        ObjectArrayList<PruningVertex> upList = prunningVertices.get(baseNode);

        if (upList != null) {
            upList.forEach(new ObjectProcedure<PruningVertex>() {
                @Override
                public void apply(PruningVertex vertex) {
                    int w = vertex.vertex;
                    double weight = vertex.weight;

                    IntObjectMap<BucketEntry> wBuckets = buckets.get(w);
                    if (wBuckets != null) {
                        wBuckets.forEach(new IntObjectProcedure<BucketEntry>() {
                            @Override
                            public void apply(int i, BucketEntry wBucket) {
                                if (wBucket.weight > weight + weights[wBucket.idx]) {
                                    wBuckets.remove(i);
                                }
                            }
                        });
                    }
                }
            });
        }
    }

    private void discoverBucketEntriesToCopy(int baseNode, IntObjectMap<IntObjectHashMap<Vertex>> vertices,
                                             double[] distances, long[] times, double[] weights,
                                             IntObjectMap<IntObjectMap<BucketEntry>> buckets,
                                             IntIntHashMap terminals) {


        Vertex shortcut = shortcutSelfVertices.get(baseNode);

        //final double shortcutWeight = (shortcut != null) ? shortcut.weight : 0;
        //final double shortcutDistance= (shortcut != null) ? shortcut.distance : 0;
        //final long shortcutTime= (shortcut != null) ? shortcut.time : 0;

        IntObjectHashMap<Vertex> downList = vertices.get(baseNode);
        if (downList != null) {

            downList.forEach(new IntObjectProcedure<Vertex>() {

                @Override
                public void apply(int i, Vertex vertex) {
                    int w = vertex.base;
                    double weight = vertex.weight;

                    IntObjectMap<BucketEntry> bucketTargets = buckets.get(w);
                    if (bucketTargets != null) {

                        bucketTargets.forEach(new IntObjectProcedure<BucketEntry>() {
                            @Override
                            public void apply(int t, BucketEntry entry) {

                                int targetIdx = entry.idx;
                                double pathWeight = entry.weight;

                                terminals.putIfAbsent(t, targetIdx);

                                double currentWeight = pathWeight + weight;
                                System.out.println("Bucket Current:" + weights[targetIdx]);
                                System.out.println("Bucket : " + baseNode + " -> " + w + " -> " + entry.weight + " - " + vertex.weight);
                                if (currentWeight < weights[targetIdx] || weights[targetIdx] == 0) {

                                    System.out.println("Bucket Saved: " + baseNode + " -> " + w + " -> " + entry.weight + " - " + vertex.weight);

                                    weights[targetIdx] = currentWeight;
                                    distances[targetIdx] = entry.distance + vertex.distance;
                                    times[targetIdx] = entry.time + vertex.time;
                                }
                            }
                        });
                    }
                }
            });
        }


    }

    private void copyBucketsEntriesToCurrentVertex(int baseNode,
                                                   IntObjectMap<IntObjectMap<BucketEntry>> buckets,
                                                   double[] distances, long[] times, double[] weights, boolean saveBestPath,
                                                   DistanceMatrix dm, IntObjectMap<IntArrayList> targets,
                                                   IntIntHashMap terminals) {
        terminals.forEach(new IntIntProcedure() {
            @Override
            public void apply(int source, int sourceIdx) {
                IntObjectMap<BucketEntry> b = buckets.get(baseNode);

                if (b == null) {
                    b = new IntObjectHashMap<>();
                    buckets.put(baseNode, b);
                }

                BucketEntry entry = new BucketEntry(weights[sourceIdx], times[sourceIdx], distances[sourceIdx], sourceIdx);
                b.put(source, entry);

                    System.out.println("Saving Bucket:" + baseNode + " - " +  source + " -> " + entry);

                if (saveBestPath) {
                    saveBestPath(source, sourceIdx, baseNode, entry, targets, dm);
                }

                weights[sourceIdx] = Double.POSITIVE_INFINITY;
                distances[sourceIdx] = Double.POSITIVE_INFINITY;
                times[sourceIdx] = Long.MAX_VALUE;
            }
        });
    }

    private void simultaneousBucketInitializationBackward(List<Snap> targets) {

        IntIntHashMap terminals = new IntIntHashMap();

        double[] distances = new double[targets.size()];
        long[] times = new long[targets.size()];
        double[] weights = new double[targets.size()];

        System.out.println("##################################");
        System.out.println("BACKWARD");

        int previousNode = -1;

        while (!heap.isEmpty()) {

            //System.out.println("Heap:" + heap);

            MatrixEntry rn = heap.poll();
            int baseNode = rn.adjNode;

            if(baseNode != previousNode && previousNode != -1){

                //initializePruningVertices(baseNode, outEdgeExplorerNoVirtual, false);

                discoverBucketEntriesToCopy(previousNode, inVertices, distances, times, weights, backwardBuckets, terminals);

                //applyRetrospectivePruningAlgorithm(previousNode, weights, backwardBuckets);

                copyBucketsEntriesToCurrentVertex(previousNode, backwardBuckets, distances, times, weights, false, null, null, terminals);

                terminals.clear();
            }

            //if(rn.adjNode == 2865619){
              //  System.out.println("## Backward " + rn);
            //}
            //System.out.println(" ** Processing Backward " + baseNode + " - " + rn.level);
            //System.out.println(" Processing Backward " + rn);

            initializeVertices(rn, inEdgeExplorerNoVirtual, inVertices, true);
            previousNode = baseNode;

        }

    }

    

    private void simultaneousBucketInitializationForward(List<Snap> sources,
                                                         DistanceMatrix dm, IntObjectMap<IntArrayList> targets) {

        System.out.println("##################################");
        System.out.println("FORWARD");
        double[] distances = new double[sources.size()];
        long[] times = new long[sources.size()];
        double[] weights = new double[sources.size()];

        int previousNode = -1;

        while (!heap.isEmpty()) {

            MatrixEntry rn = heap.poll();
            int baseNode = rn.adjNode;
            //System.out.println("** Processing Forward: " + baseNode);
            //System.out.println(" Processing Forward " + rn);

            if(baseNode != previousNode && previousNode != -1){
                //initializePruningVertices(baseNode, inEdgeExplorerNoVirtual, true);

                discoverBucketEntriesToCopy(previousNode, outVertices, distances, times, weights, forwardBuckets, terminals);

                //applyRetrospectivePruningAlgorithm(previousNode, weights, forwardBuckets);

                copyBucketsEntriesToCurrentVertex(previousNode, forwardBuckets, distances, times, weights,
                        true, dm, targets, terminals);

                terminals.clear();
            }

            initializeVertices(rn, outEdgeExplorerNoVirtual, outVertices, false);
            previousNode = baseNode;

        }
    }

    private void saveBestPath(int sourceNode, int idxSource, int currNode, BucketEntry forwardEntry,
                              IntObjectMap<IntArrayList> targets, DistanceMatrix dm) {

        final IntObjectMap<BucketEntry> backwardEntries = backwardBuckets.get(currNode);


        if (backwardEntries != null) {

            //System.out.println("$$$$$ Saving: " + currNode + " - " + backwardEntries);

            backwardEntries.forEach(new IntObjectProcedure<BucketEntry>() {
                @Override
                public void apply(int target, BucketEntry entry) {

                    if (sourceNode != target) {

                        long uniqueId = PairingUtils.pair(sourceNode, target);

                        final double savedWeight = tentativeWeights.get(uniqueId);
                        final double currentWeight = forwardEntry.weight + entry.weight;

                        if (savedWeight == 0.0 || (currentWeight < savedWeight)) {


                                System.out.println(currNode + " $$$$$ Current:" + currentWeight + " Saved: " + savedWeight);


                            final long time = forwardEntry.time + entry.time;
                            final double distance = forwardEntry.distance + entry.distance;
                            tentativeWeights.put(uniqueId, currentWeight);

                            //System.out.println(" $$$$$ Saved Distance:" + distance);

                            IntArrayList targetsIdxs = targets.get(target);
                            int[] buffer = targetsIdxs.buffer;
                            int size = targetsIdxs.size();

                            for (int i = 0; i < size; i++) {
                                dm.setCell(idxSource, buffer[i], distance, time);
                            }
                        }
                    }
                }
            });
        }
    }
}