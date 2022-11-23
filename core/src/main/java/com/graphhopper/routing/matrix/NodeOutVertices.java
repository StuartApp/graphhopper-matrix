package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.ObjectArrayList;


public class NodeOutVertices {

    private ObjectArrayList<Vertex> outs;
    private ObjectArrayList<Vertex> selfs;

    public NodeOutVertices(){
        this.outs = new ObjectArrayList<>();
        this.selfs = new ObjectArrayList<>();
    }

    public void add(Vertex v){

        if(v.isSelfLoop()){
            selfs.add(v);
        }else{
            outs.add(v);
        }
    }

    public boolean hasSelfLoops(){
        return !this.selfs.isEmpty();
    }

    public boolean hasOutValues(){
        return !this.outs.isEmpty();
    }

    public ObjectArrayList<Vertex> outVertexes() {
        return this.outs;
    }

    public ObjectArrayList<Vertex> selfVertexes(){
        return this.selfs;
    }

}
