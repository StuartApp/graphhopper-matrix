package com.graphhopper.routing.matrix;

public class Vertex {
    public int base;
    public double weight;
    public long time;
    public double distance;
    public int adj;


    public Vertex(int base, int adj,double weight, long time, double distance) {
        this.base = base;
        this.weight = weight;
        this.time = time;
        this.distance = distance;
        this.adj = adj;
    }
}
