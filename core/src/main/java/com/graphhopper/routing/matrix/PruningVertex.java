package com.graphhopper.routing.matrix;

public class PruningVertex {
    public int vertex;
    public double weight;

    public PruningVertex(int vertex, double weight) {
        this.vertex = vertex;
        this.weight = weight;
    }
}
