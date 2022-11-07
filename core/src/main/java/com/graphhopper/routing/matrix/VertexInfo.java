package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.IntHashSet;
import com.carrotsearch.hppc.IntSet;

import java.util.Objects;

public class VertexInfo {
    protected double turnCost;
    protected IntSet nodes;

    public VertexInfo() {
        this.turnCost = Double.POSITIVE_INFINITY;
        this.nodes = new IntHashSet();
    }

    public IntSet getNodes() {
        return nodes;
    }

    public void setTurnCost(double turnCost) {
        this.turnCost = turnCost;
    }

    public double getTurnCost() {
        return turnCost;
    }

    public void addNode(int node){
        this.nodes.add(node);
    }
}
