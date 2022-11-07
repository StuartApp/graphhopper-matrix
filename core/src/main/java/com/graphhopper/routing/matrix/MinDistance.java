package com.graphhopper.routing.matrix;

public class MinDistance {
    public double weight;
    public long time;
    public double distance;




   public MinDistance(double weight, long time, double distance){
       this.weight = weight;
       this.time = time;
       this.distance = distance;

   }

    @Override
    public String toString() {
        return "MinDistance{" +
                " weight=" + weight +
                ", time=" + time +
                ", distance=" + distance +
                '}';
    }
}
