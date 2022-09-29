package com.graphhopper.routing.matrix;

public class SBIEntry {
    public int edge;
    public double weight;
    public long time;
    public double distance;
    public int adj;
    public int incEdge;
    public int traversalId;
    public int snapClosestNode;


    public SBIEntry(int traversalId,int edge, int adj, double weight, long time, double distance) {
        this.traversalId = traversalId;
        this.edge = edge;
        this.weight = weight;
        this.time = time;
        this.distance = distance;
        this.adj = adj;
    }

    @Override
    public String toString() {
        return "SBIEntry{" +
                "edge=" + edge +
                ", weight=" + weight +
                ", time=" + time +
                ", distance=" + distance +
                ", adj=" + adj +
                ", incEdge=" + incEdge +
                ", traversalId=" + traversalId +
                '}';
    }
}
