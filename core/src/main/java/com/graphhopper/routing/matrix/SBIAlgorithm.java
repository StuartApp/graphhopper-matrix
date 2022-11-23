package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.*;
import com.carrotsearch.hppc.procedures.IntObjectProcedure;
import com.carrotsearch.hppc.procedures.IntProcedure;
import com.carrotsearch.hppc.procedures.ObjectProcedure;
import com.graphhopper.storage.RoutingCHGraph;
import com.graphhopper.util.PairingUtils;

public class SBIAlgorithm {

    RoutingCHGraph graph;

    IntObjectMap<IntArrayList> targetIndexesNodes;
    IntObjectMap<IntArrayList> sourcesIndexesNodes;

    IntObjectMap<NodeTerminals> nodeTerminals;

    DistanceMatrix dm;

    LongDoubleMap shortestRoutes;


    //IntDoubleMap forwardComparator = new IntDoubleHashMap();
    //IntIntMap forwardComparatorOrder = new IntIntHashMap();
    //IntDoubleMap backwardComparator = new IntDoubleHashMap();
    //IntIntMap backwardComparatorOrder = new IntIntHashMap();


    public SBIAlgorithm(RoutingCHGraph graph, DistanceMatrix dm ) {

        this.graph = graph;
        this.targetIndexesNodes = new IntObjectHashMap<>(dm.getNumberOfDestinations());
        this.sourcesIndexesNodes = new IntObjectHashMap<>(dm.getNumberOfOrigins());

        int size = Math.min(Math.max(200, graph.getNodes() / 10), 150_000);
        this.nodeTerminals = new IntObjectHashMap<>(size);
        this.dm = dm;
        this.shortestRoutes = new LongDoubleHashMap(dm.getNumberOfDestinations() * dm.getNumberOfOrigins());

        /*
        this.forwardComparator.put(2966296,8010.0);
        this.forwardComparator.put(2966296,69626.0);
        this.forwardComparator.put(150209,18021.0);
        this.forwardComparator.put(490,10150.0);

        this.forwardComparatorOrder.put(2966296,2966296);
        this.forwardComparatorOrder.put(150209,2966296);
        this.forwardComparatorOrder.put(490,150209);


        this.backwardComparator.put(2966296,229917.0);
        this.backwardComparator.put(1612521,162879.0);
        this.backwardComparator.put(65376,63344.0);
        this.backwardComparator.put(284578,28668.0);
        this.backwardComparator.put(284580,2042.0);

        this.backwardComparatorOrder.put(1612521,2966296);
        this.backwardComparatorOrder.put(65376,1612521);
        this.backwardComparatorOrder.put(284578,65376);
        this.backwardComparatorOrder.put(284580,284578);

         */

    }

    public void addSource(int sourceNode, int sourceIdx){

        IntArrayList sourceIdxs = sourcesIndexesNodes.get(sourceNode);
        if(sourceIdxs == null){
            sourceIdxs = new IntArrayList();
            sourcesIndexesNodes.put(sourceNode,sourceIdxs);
            initializeTerminals(sourceNode,sourceIdx,false);
            findRoutes(sourceNode,null);
        }

        sourceIdxs.add(sourceIdx);
    }

    public void addTarget(int targetNode, int targetIdx){

        IntArrayList targetIdxs = targetIndexesNodes.get(targetNode);
        if(targetIdxs == null){
            targetIdxs = new IntArrayList();
            targetIndexesNodes.put(targetNode,targetIdxs);
            initializeTerminals(targetNode,targetIdx,true);
        }

        targetIdxs.add(targetIdx);
    }

    private NodeTerminals obtainNodeTerminals(int node){
        NodeTerminals terminals = nodeTerminals.get(node);
        if(terminals == null){
            terminals = new NodeTerminals(node);
            nodeTerminals.put(node,terminals);
        }

        return terminals;
    }

    private void initializeTerminals(int terminal, int terminalIdx, boolean reverse){

        NodeTerminals terminals = obtainNodeTerminals(terminal);
        Terminal initial = new Terminal(0,0,0,terminal,terminalIdx);
        terminals.addInitialTerminal(initial,reverse,-1,-1,-1);
    }

    //NO REMOVE - Used for debugging routes errors purposes
    private void checkResult(VertexWithTerminals in,
                             int node, int adjNode, double weight, long time, double distance, boolean reverse){
        /*
        if(reverse){
            if(backwardComparatorOrder.containsKey(node)){
                if(backwardComparatorOrder.get(node) == adjNode){
                    double expected = backwardComparator.get(adjNode);
                    if(expected == time){
                        System.out.println("BACKWARD: " + node + " -> " + adjNode + " : " + time + " = " + expected);
                        in.getTerminalsList().forEach(new IntObjectProcedure<Terminal>() {
                            @Override
                            public void apply(int i, Terminal terminal) {
                                System.out.println(" Target : " + terminal.nodeIdx + " - " + terminal.time);
                            }
                        });
                    }
                }
            }
        }else{
            if(forwardComparatorOrder.containsKey(node)){
                if(forwardComparatorOrder.get(node) == adjNode){
                    double expected = forwardComparator.get(adjNode);
                    if(expected == time){
                        System.out.println("FORWARD: " + in.adjNode + " -> " + adjNode + " : " + time + " = " + expected);
                        in.getTerminalsList().forEach(new IntObjectProcedure<Terminal>() {
                            @Override
                            public void apply(int i, Terminal terminal) {
                                System.out.println(" Source : " + terminal.nodeIdx + " - " + terminal.time);
                            }
                        });
                    }
                }
            }

        }

         */
    }

    public void addInitialOutVertex(int node, int terminal, int terminalIdx, double weight, long time,
                                    double distance, boolean reverse, int origEdgeId, int origEdgeFirst, int origEdgeLast) {

        NodeTerminals terminals = obtainNodeTerminals(node);
        Terminal initial = new Terminal(weight,time,distance,terminal,terminalIdx);
        terminals.addInitialTerminal(initial,reverse,origEdgeId,origEdgeFirst,origEdgeLast);
        if(reverse) {
            findRoutes(node,null);
        }
    }

    private void addOutVertex(VertexWithTerminals in, Vertex out, boolean reverse){

        //NOT REMOVE!!! checkResult(in, out.baseNode, out.adjNode,out.weight,out.time,out.distance,reverse);

        NodeTerminals nodeTerminals = obtainNodeTerminals(out.adjNode);
        nodeTerminals.addTerminal(in,out,reverse);

    }

    private ObjectContainer<VertexWithTerminals> EMPTY_INS = new ObjectArrayList<>();


    private ObjectContainer<VertexWithTerminals> getIns(int node, boolean reverse){

        NodeTerminals nt = nodeTerminals.get(node);

        if(nt == null){
            return EMPTY_INS;
        }else{
            return nt.getTerminals(reverse);
        }
    }

    private double calculateTurnCost(int inVertexId, int outVertexId, int node, boolean reverse){

            return reverse
                    ? graph.getTurnWeight(outVertexId, node, inVertexId)
                    : graph.getTurnWeight(inVertexId, node, outVertexId);
    }

    private boolean isAccessible(double turnCost){
        return Double.isFinite(turnCost);
    }

    public void addOuts(int node, NodeOutVertices outs, boolean reverse){

        if(outs.hasOutValues()){

            getIns(node,reverse).forEach(new ObjectProcedure<VertexWithTerminals>() {
                @Override
                public void apply(VertexWithTerminals in) {
                    int inOrigEdgeId = reverse ? in.origEdgeFirst : in.origEdgeLast;

                    outs.outVertexes().forEach(new ObjectProcedure<Vertex>() {
                        @Override
                        public void apply(Vertex out) {
                            //If an in vertex origEdgeId is equal of this vertex origEdgeId we discard it
                            if(in.origEdgeId != out.origEdgeId){
                                int outOrigEdgeId = reverse ? out.origEdgeLast : out.origEdgeFirst;

                                double cost = calculateTurnCost(inOrigEdgeId,outOrigEdgeId,node,reverse);

                                //Out is accessible
                                if(isAccessible(cost)){

                                    Vertex outWithCost = out.withTurnCost(cost);
                                    addOutVertex(in,outWithCost,reverse);
                                }

                                //Check Selfs
                                if(outs.hasSelfLoops()){

                                    outs.selfVertexes().forEach(new ObjectProcedure<Vertex>() {
                                        @Override
                                        public void apply(Vertex self) {
                                            if(in.origEdgeId != self.origEdgeId){
                                                int selfOrigEdgeId = reverse ? self.origEdgeLast : self.origEdgeFirst;

                                                double costInToSelf = calculateTurnCost(inOrigEdgeId,selfOrigEdgeId,node,reverse);

                                                //In to self is accessible
                                                if(isAccessible(costInToSelf)){

                                                    double costSelfToOut = calculateTurnCost(selfOrigEdgeId,outOrigEdgeId,node,reverse);
                                                    //Self to out is accessible
                                                    if(isAccessible((costSelfToOut))){
                                                        Vertex outWithSelfCost = out.withSelf(self).withSelfTurnCost(costInToSelf,costSelfToOut);
                                                        addOutVertex(in,outWithSelfCost,reverse);
                                                    }
                                                }
                                            }
                                        }
                                    });


                                }
                            }
                        }
                    });
                }
            });
        }
    }

    private void saveToDistanceMatrix(int sourceIdx, int targetNode, long time, double distance){
        IntArrayList targets = targetIndexesNodes.get(targetNode);
            targets.forEach(new IntProcedure() {
                @Override
                public void apply(int targetIdx) {
                    dm.setCell(sourceIdx,targetIdx, distance, time);
                }
            });
    }

    private void saveShortRoutes(VertexWithTerminals forward, VertexWithTerminals backward, int node){

        int fEdgeId = forward.origEdgeLast;
        int bEdgeId = backward.origEdgeFirst;

        double turnCost = calculateTurnCost(fEdgeId,bEdgeId,node,false);


        if(isAccessible(turnCost)){

            forward.getTerminals().forEach(new IntObjectProcedure<Terminal>() {
                @Override
                public void apply(int i, Terminal ft) {

                    backward.getTerminals().forEach(new IntObjectProcedure<Terminal>() {
                        @Override
                        public void apply(int i, Terminal bt) {
                            int sourceNode = ft.node;
                            int sourceIdx = ft.nodeIdx;
                            int targetNode = bt.node;

                            long key = PairingUtils.pair(sourceNode, targetNode);
                            double possible = ft.weight + bt.weight + turnCost;
                            double current = shortestRoutes.get(key);

                            if(sourceNode == targetNode){
                                saveToDistanceMatrix(sourceIdx,targetNode,0,0);
                            }else if (current == 0.0 || current > possible){
                                shortestRoutes.put(key,possible);
                                long time = ft.time + bt.time;
                                double distance = ft.distance + bt.distance;
                                saveToDistanceMatrix(sourceIdx,targetNode,time,distance);
                            }
                        }
                    });

                }
            });
        }
    }

    private void saveShortRoutesWithSelfLoop(VertexWithTerminals forward, VertexWithTerminals backward, int node,
                                             NodeOutVertices outs){

        int fEdgeId = forward.origEdgeLast;
        int bEdgeId = backward.origEdgeFirst;

        outs.selfVertexes().forEach(new ObjectProcedure<Vertex>() {
            @Override
            public void apply(Vertex self) {
                double turnCostFromForwardToSelf =  calculateTurnCost(fEdgeId,self.origEdgeFirst,node,false);
                if(isAccessible(turnCostFromForwardToSelf)){
                    double turnCostFromSelfToBackward =  calculateTurnCost(self.origEdgeLast,bEdgeId,node,false);
                    if(isAccessible(turnCostFromSelfToBackward)){

                        //
                        forward.getTerminals().forEach(new IntObjectProcedure<Terminal>() {
                            @Override
                            public void apply(int i, Terminal ft) {

                                backward.getTerminals().forEach(new IntObjectProcedure<Terminal>() {
                                    @Override
                                    public void apply(int i, Terminal bt) {
                                        int sourceNode = ft.node;
                                        int sourceIdx = ft.nodeIdx;
                                        int targetNode = bt.node;

                                        long key = PairingUtils.pair(sourceNode, targetNode);
                                        double possible = ft.weight + bt.weight + + self.weight + turnCostFromForwardToSelf + turnCostFromSelfToBackward;
                                        double current = shortestRoutes.get(key);

                                        if (current == 0.0 || current > possible){
                                            shortestRoutes.put(key,possible);
                                            long time = ft.time + bt.time + self.time;
                                            double distance = ft.distance + bt.distance + self.distance;
                                            saveToDistanceMatrix(sourceIdx,targetNode,time,distance);
                                        }
                                    }
                                });

                            }
                        });
                    }
                }
            }
        });
    }

    public void findRoutes(int node,NodeOutVertices outs){

        NodeTerminals nodeTerminals = obtainNodeTerminals(node);

        if(nodeTerminals.hasPossibleShortRoutes()){

            nodeTerminals.getForwardTerminals().forEach(new ObjectProcedure<VertexWithTerminals>() {
                @Override
                public void apply(VertexWithTerminals f) {

                    nodeTerminals.getBackwardTerminals().forEach(new ObjectProcedure<VertexWithTerminals>() {
                        @Override
                        public void apply(VertexWithTerminals b) {
                            saveShortRoutes(f,b,node);
                            if(outs != null && outs.hasSelfLoops()){
                                saveShortRoutesWithSelfLoop(f,b,node,outs);
                            }
                        }
                    });
                }
            });

        }
    }


    private class NodeTerminals {

        private int node;

        private LongObjectMap<VertexWithTerminals> forwardTerminals;
        private LongObjectMap<VertexWithTerminals> backwardTerminals;

        public NodeTerminals(int node) {
            this.node = node;
            this.forwardTerminals = new LongObjectHashMap<>();
            this.backwardTerminals = new LongObjectHashMap<>();
        }


        public ObjectContainer<VertexWithTerminals> getTerminals(boolean reverse){
            if(reverse){
                return backwardTerminals.values();
            }else{
                return forwardTerminals.values();
            }
        }

        public boolean hasPossibleShortRoutes(){
            return !this.backwardTerminals.isEmpty() && !this.forwardTerminals.isEmpty();
        }

        public void addInitialTerminal(Terminal terminal, boolean reverse, int originalEdgeId, int origEdgeFirst, int origEdgeLast){

            VertexWithTerminals vertex = new VertexWithTerminals(originalEdgeId,origEdgeFirst,origEdgeLast);

            long key = PairingUtils.pair(terminal.node, vertex.origEdgeId);

            vertex.addTerminal(terminal);

            if(reverse){
                backwardTerminals.put(key,vertex);
            }else{
                forwardTerminals.put(key,vertex);
            }
        }

        public void addTerminal(VertexWithTerminals vt, Vertex out, boolean reverse){

            LongObjectMap<VertexWithTerminals> search = (reverse) ? backwardTerminals : forwardTerminals;

            vt.getTerminals().forEach(new IntObjectProcedure<Terminal>() {
                @Override
                public void apply(int i, Terminal t) {
                    int edgeId = (reverse) ? out.origEdgeFirst : out.origEdgeLast;

                    long key = PairingUtils.pair(t.node,edgeId);

                    VertexWithTerminals current = search.get(key);
                    if(current == null){
                        current = VertexWithTerminals.createFromOut(out);
                        search.put(key,current);
                    }

                    current.addTerminal(t.with(out.weight,out.time,out.distance));
                }
            });
        }


        public ObjectContainer<VertexWithTerminals> getBackwardTerminals() {
            return getTerminals(true);
        }

        public ObjectContainer<VertexWithTerminals> getForwardTerminals() {
            return getTerminals(false);
        }
    }
}
