package com.graphhopper.routing.matrix;

public class Vertex {
    public int base;
    public double weight;
    public long time;
    public double distance;


    public Vertex(int base,double weight, long time, double distance) {
        this.base = base;
        this.weight = weight;
        this.time = time;
        this.distance = distance;
    }

    @Override
    public String toString() {
        return "Vertex{" +
                "base=" + base +
                ", weight=" + weight +
                ", time=" + time +
                ", distance=" + distance +
                '}';
    }
}
