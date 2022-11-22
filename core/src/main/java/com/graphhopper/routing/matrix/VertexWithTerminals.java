package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.IntObjectHashMap;
import com.carrotsearch.hppc.IntObjectMap;
import com.carrotsearch.hppc.procedures.IntObjectProcedure;
import com.graphhopper.storage.RoutingCHEdgeIterator;

public class VertexWithTerminals {

    public int baseNode;
    public int adjNode;
    public int edge;
    public int origEdgeId;
    public int origEdgeFirst;
    public int origEdgeLast;
    private IntObjectMap<Terminal> terminals;


    public VertexWithTerminals(RoutingCHEdgeIterator iter){

        this(iter.getBaseNode(), iter.getEdge(), iter.getAdjNode(),iter.getOrigEdge(),
                iter.getOrigEdgeFirst(), iter.getOrigEdgeLast());

    }

    public VertexWithTerminals(int base, int edge, int adj, int origEdgeId, int origEdgeFirst, int origEdgeLast) {
        this.baseNode = base;
        this.edge = edge;
        this.adjNode = adj;
        this.origEdgeId = origEdgeId;
        this.origEdgeFirst = origEdgeFirst;
        this.origEdgeLast = origEdgeLast;
        this.terminals = new IntObjectHashMap<>();
    }

    public void addTerminal(Terminal terminal){

        Terminal current = terminals.get(terminal.nodeIdx);
        if(current == null || current.weight >= terminal.weight){
            terminals.put(terminal.nodeIdx,terminal);
        }

    }

    public VertexWithTerminals combine(RoutingCHEdgeIterator iter, double turnCost, boolean reverse){

       VertexWithTerminals vt = new VertexWithTerminals(iter);

       long turnCostMillis = turnCostInMillis(turnCost);

       double weight = iter.getWeight(reverse);
       long time = iter.getTime(reverse);
       double distance = iter.getDistance();

       this.terminals.forEach(new IntObjectProcedure<Terminal>() {
           @Override
           public void apply(int i, Terminal terminal) {

               double w = terminal.weight + turnCost + weight;
               long t = terminal.time + turnCostMillis + time;
               double d = terminal.distance + distance;

               addTerminal(terminal.with(w, t, d));
           }
       });

       return vt;
    }



    private long turnCostInMillis(double turnCost){
        return (long) turnCost * 1000;
    }

    public Terminal[] getTerminals() {
        return terminals.values().toArray(Terminal.class);
    }


    public IntObjectMap<Terminal> getTerminalsList(){
        return terminals;
    }

    public static VertexWithTerminals createFromOut(Vertex out){

        VertexWithTerminals vt = new VertexWithTerminals(out.baseNode, out.edge, out.adjNode, out.origEdgeId,
               out.origEdgeFirst, out.origEdgeLast);

        return vt;

    }

    @Override
    public String toString() {
        return "VertexWithTerminals{" +
                "baseNode=" + baseNode +
                ", adjNode=" + adjNode +
                ", edge=" + edge +
                ", origEdgeId=" + origEdgeId +
                ", origEdgeFirst=" + origEdgeFirst +
                ", origEdgeLast=" + origEdgeLast +
                '}';
    }
}
