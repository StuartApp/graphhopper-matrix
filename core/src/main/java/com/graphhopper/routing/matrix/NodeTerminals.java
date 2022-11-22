package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.LongObjectHashMap;
import com.carrotsearch.hppc.LongObjectMap;
import com.graphhopper.util.PairingUtils;

public class NodeTerminals {

    private LongObjectMap<VertexWithTerminals> forwardTerminals;
    private LongObjectMap<VertexWithTerminals> backwardTerminals;

    public NodeTerminals() {
        this.forwardTerminals = new LongObjectHashMap<>();
        this.backwardTerminals = new LongObjectHashMap<>();
    }


    public VertexWithTerminals[] getTerminals(boolean reverse){
        if(reverse){
            return backwardTerminals.values().toArray(VertexWithTerminals.class);
        }else{
            return forwardTerminals.values().toArray(VertexWithTerminals.class);
        }
    }

    public boolean hasPossibleShortRoutes(){
        return !this.backwardTerminals.isEmpty() && !this.forwardTerminals.isEmpty();
    }

    public void addInitialTerminal(Terminal terminal, boolean reverse, int originalEdgeId, int origEdgeFirst, int origEdgeLast){

        VertexWithTerminals vertex = new VertexWithTerminals(terminal.node,-1,-1,originalEdgeId
                ,origEdgeFirst,origEdgeLast);

        long key = PairingUtils.pair(terminal.node, vertex.origEdgeId);

        vertex.addTerminal(terminal);

        if(reverse){
            backwardTerminals.put(key,vertex);
        }else{
            forwardTerminals.put(key,vertex);
        }
    }

    public void addTerminal(VertexWithTerminals vt, Vertex out, boolean reverse){

        Terminal[] terminals = vt.getTerminals();
        int size = terminals.length;

        LongObjectMap<VertexWithTerminals> search = (reverse) ? backwardTerminals : forwardTerminals;

        for(int i = 0; i < size; i++ ){
            Terminal t = terminals[i];

            int edgeId = (reverse) ? out.origEdgeFirst : out.origEdgeLast;

            long key = PairingUtils.pair(t.node,edgeId);

            VertexWithTerminals current = search.get(key);
            if(current == null){
                current = VertexWithTerminals.createFromOut(out);
                search.put(key,current);
            }

            current.addTerminal(t.with(out.weight,out.time,out.distance));

        }

    }

    public VertexWithTerminals[] getForwardTerminals() {
        return forwardTerminals.values().toArray(VertexWithTerminals.class);
    }

    public VertexWithTerminals[] getBackwardTerminals() {
        return backwardTerminals.values().toArray(VertexWithTerminals.class);
    }
}
